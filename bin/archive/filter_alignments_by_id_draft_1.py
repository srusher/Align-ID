#!/usr/bin/env python3
"""Python port of the tax-ID alignment filtering bash script.

Usage: filter_by_tax_id.py PREFIX BAM TAX_IDS INCLUDE_CHILDREN NODES PYTHON_SCRIPT SAM_DB

  PREFIX			prefix for output files
  BAM			   original BAM (only its header is used)
  TAX_IDS		   text file with one tax ID per line
  INCLUDE_CHILDREN  'true'/'True' to also keep every descendant of each tax ID
  NODES			 NCBI nodes.dmp (only used when INCLUDE_CHILDREN is true)
  PYTHON_SCRIPT	 child-node lookup script, run as `python PYTHON_SCRIPT nodes_mod.dmp TAX_IDS`;
					it must write tax_id_child_nodes.tsv (only used when INCLUDE_CHILDREN is true)
  SAM_DB			database created by the import script (contains sam_complete)

Keeps every read (all of its alignments) that has at least one alignment whose tax_id is in
the wanted list, then writes PREFIX-primary-all-filtered_sorted.bam.

sqlite3 (stdlib) replaces the sqlite container and pysam (bundles samtools/htslib)
replaces the samtools container, so no Singularity/Apptainer is needed.
Requires: pysam.
"""

import os
import sqlite3
import subprocess
import sys
import time
from datetime import datetime

import pysam


THREADS = 16


# =============================================================================
# Helper functions
# =============================================================================

def now():

	return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg):

	print(f"{msg}: {now()}", flush=True)


def format_duration(seconds):

	seconds = int(round(seconds))

	h, rem = divmod(seconds, 3600)
	m, s = divmod(rem, 60)

	return f"{h}h {m}m {s}s"


def bam_header_text(bam):
	"""Equivalent of `samtools view -H BAM`."""

	with pysam.AlignmentFile(bam, "rb", check_sq=False) as af:

		text = str(af.header)

	return text if text.endswith("\n") or not text else text + "\n"


def read_tax_ids(path):
	"""One tax ID per line; blank lines are ignored."""

	with open(path) as fh:

		return [line.strip() for line in fh if line.strip()]


def trim_nodes_file(nodes_path, out_path):
	"""Equivalent of `sed 's/\\t|\\t/\\t/g' NODES | awk -F'\\t' -v OFS='\\t' '{print $1,$2}'`.

	Reduces nodes.dmp to two tab-separated columns: tax_id, parent_tax_id.
	"""

	with open(nodes_path) as fin, open(out_path, "w") as fout:

		for line in fin:

			fields = line.rstrip("\n").replace("\t|\t", "\t").split("\t")

			fields += [""] * (2 - len(fields))

			fout.write(f"{fields[0]}\t{fields[1]}\n")


def run_child_node_script(python_script, nodes_mod, my_tax_ids):
	"""Equivalent of `python PYTHON_SCRIPT nodes_mod.dmp MY_TAX_IDS`.

	The script is expected to write tax_id_child_nodes.tsv in the working directory.
	It runs with the same interpreter as this script (previously a python container).
	"""

	subprocess.run([sys.executable, python_script, nodes_mod, my_tax_ids], check=True)

	child_file = "tax_id_child_nodes.tsv"

	if not os.path.exists(child_file):

		sys.exit(f"ERROR: {python_script} did not create {child_file}")

	return child_file


# =============================================================================
# Main
# =============================================================================

def main():

	if len(sys.argv) < 8:

		sys.exit(__doc__)

	prefix, bam, my_tax_ids, include_children, nodes, python_script, sam_db = sys.argv[1:8]

	primary_all = f"{prefix}-primary_all.sam"
	read_pairs_count_file = f"{prefix}_read_pairs_count.tsv"

	# -------------------------------------------------------------------------
	# Tax IDs to keep
	# -------------------------------------------------------------------------
	if include_children in ("true", "True"):

		log("Collecting child tax IDs START")

		trim_nodes_file(nodes, "nodes_mod.dmp")

		tax_ids_i_want = run_child_node_script(python_script, "nodes_mod.dmp", my_tax_ids)

	else:

		tax_ids_i_want = my_tax_ids

	tax_ids = read_tax_ids(tax_ids_i_want)

	log("Collecting child tax IDs END")
	# -------------------------------------------------------------------------
	# BAM header
	# -------------------------------------------------------------------------
	log("Parsing BAM headers START")

	with open(primary_all, "w") as fh:

		fh.write(bam_header_text(bam))

	log("Parsing BAM headers END")
	# -------------------------------------------------------------------------
	# Database connection
	# -------------------------------------------------------------------------

	# autocommit mode so each statement is committed immediately
	conn = sqlite3.connect(sam_db, timeout=60, isolation_level=None)

	conn.execute("PRAGMA cache_size = -50000000")
	conn.execute("PRAGMA temp_store = MEMORY")

	# -------------------------------------------------------------------------
	# Load the wanted tax IDs (clearing any past table first)
	# -------------------------------------------------------------------------
	conn.executescript("""
		DROP TABLE IF EXISTS tax_ids_i_want;

		CREATE TABLE tax_ids_i_want (tax_id TEXT);
	""")

	conn.execute("BEGIN")
	conn.executemany("INSERT INTO tax_ids_i_want VALUES (?)", [(t,) for t in tax_ids])
	conn.execute("COMMIT")

	# -------------------------------------------------------------------------
	# Pre-filter read-pair count
	# -------------------------------------------------------------------------
	pre_count = conn.execute("SELECT COUNT(DISTINCT read_id) FROM sam_complete").fetchone()[0]

	# -------------------------------------------------------------------------
	# Add indexes to the tax_id column in the sam_complete and tax_ids_i_want tables
	# -------------------------------------------------------------------------

	log("Indexing tax_id columns START")

	conn.executescript("""
	CREATE INDEX idx_sc_tax_id ON sam_complete (tax_id);
	CREATE UNIQUE INDEX idx_want ON tax_ids_i_want (tax_id);
	VACUUM;
	""")

	log("Indexing tax_id columns END")

	# -------------------------------------------------------------------------
	# Remove all reads that don't match the tax ID filter
	# -------------------------------------------------------------------------
	log("Removing all SAM entries from DB that didn't match tax ID filter START")

	conn.executescript("""
		CREATE TABLE sam_complete_new AS
		SELECT *
		FROM sam_complete
		WHERE read_id IN (
			SELECT sc.read_id
			FROM sam_complete AS sc
			INNER JOIN tax_ids_i_want as t ON (sc.tax_id = t.tax_id)
		);

		DROP TABLE sam_complete;

		ALTER TABLE sam_complete_new RENAME TO sam_complete;
	""")

	log("Removing all SAM entries from DB that didn't match tax ID filter END")

	# -------------------------------------------------------------------------
	# Post-filter read-pair count
	# -------------------------------------------------------------------------
	post_count = conn.execute("SELECT COUNT(DISTINCT read_id) FROM sam_complete").fetchone()[0]

	with open(read_pairs_count_file, "w") as fh:

		fh.write(f"{pre_count}\n")
		fh.write(f"{post_count}\n")

	# -------------------------------------------------------------------------
	# Print DB entries to SAM file
	# -------------------------------------------------------------------------
	log("Printing DB entries to SAM file START")

	def fmt(v):

		return "" if v is None else str(v)

	with open(primary_all, "a") as fh:

		for row in conn.execute("SELECT * FROM sam_complete"):

			# first 15 columns are the SAM columns (tax_id / name come after them)
			cols = [fmt(v) for v in row[:15]]
			cols += [""] * (15 - len(cols))

			# blank sequence / quality are invalid SAM - use the '*' placeholder
			for i in (9, 10):

				if cols[i] == "":

					cols[i] = "*"

			# keep the 11 mandatory fields, drop empty optional fields
			cols = cols[:11] + [c for c in cols[11:] if c != ""]

			fh.write("\t".join(cols) + "\n")

	log("Printing DB entries to SAM file END")

	# clearing out intermediate tables and reclaiming unused space
	conn.execute("VACUUM;")

	conn.close()

	# -------------------------------------------------------------------------
	# Convert SAM back to a sorted BAM
	# -------------------------------------------------------------------------
	print("Converting SAM back to BAM and sorting", flush=True)

	# --verbosity 0 silences 'unrecognized reference name "unclassified"' warnings
	pysam.sort("--verbosity", "0", "-@", str(THREADS),
			   "-o", f"{prefix}-primary-all-filtered_sorted.bam", primary_all)

	# -------------------------------------------------------------------------
	# Cleanup
	# -------------------------------------------------------------------------
	os.remove(primary_all)


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":

	start = time.monotonic()

	try:

		main()

	finally:

		print(f"Total run time: {format_duration(time.monotonic() - start)}", flush=True)