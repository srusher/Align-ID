#!/usr/bin/env python3
"""Python port of the BAM/SAM -> SQLite import bash script.

Usage: import_sam_to_db.py PREFIX BAM_OR_SAM SEQID2TAXID_DB SINGLE_END INCLUDE_PRIMARY_SEQ

  PREFIX			  prefix for output files; the database is written to PREFIX-sam.db
  BAM_OR_SAM		  input alignment file (a .sam input is sorted into PREFIX.bam first)
  SEQID2TAXID_DB	  prebuilt SQLite database containing a tax_map(seq_id, tax_id, name)
					  table; it is copied to PREFIX-sam.db and left untouched
  SINGLE_END		  'false' for paired-end data, anything else means single-end
  INCLUDE_PRIMARY_SEQ 'false' to drop sequence / quality strings, anything else keeps them

Tables written to PREFIX-sam.db:
  sam_partial_1	everything except sequence and phred (plus alignment_id, tax_id, name)
  sam_partial_2	alignment_id, sequence, phred (linked to sam_partial_1 by alignment_id)

sqlite3 (stdlib) replaces the sqlite container and pysam (bundles samtools/htslib)
replaces the samtools container, so no Singularity/Apptainer is needed.
Requires: pysam.
"""

import itertools
import os
import shutil
import sqlite3
import sys
import time
import csv
from datetime import datetime
import pysam
import shutil

tmp_dir = os.path.join(os.getcwd(), "sqlite_tmp")
os.makedirs(tmp_dir, exist_ok=True)
os.environ["SQLITE_TMPDIR"] = tmp_dir

SORT_THREADS = 4
BATCH_SIZE = 200_000


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


def trimmed_sam_rows(bam, single_end, include_primary_seq):
	"""Equivalent of the `samtools view | awk | sed | awk` trimming pipeline.

	Yields one (partial_1_row, partial_2_row) pair per kept alignment:

	partial_1_row (14 values):
		alignment_id, read_id, flag, ref_id, pos, mapq, cigar, r_next, p_next,
		t_length, nm_tag, ms_tag, as_tag, tp_tag

	partial_2_row (3 values):
		alignment_id, sequence, phred

	The optional fields are positional, exactly as in the awk script
	($12, $13, $14 and $16). The double-quote escaping done by `sed` is not needed
	because rows are inserted with bound parameters instead of `.import`.
	"""

	alignment_id = 0

	with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=4) as alignment_file:

		for alignment in alignment_file:

			fields = alignment.to_string().split("\t")

			# awk returns an empty string for fields that don't exist
			if len(fields) < 16:

				fields += [""] * (16 - len(fields))

			flag = alignment.flag
			tp_tag = fields[15]

			# ---- filtering (awk '$2 < 2000' / the single-end tp tag filter)
			if not single_end:

				if not flag < 2000:

					continue

			else:

				if not (tp_tag == "tp:A:P" or tp_tag == "tp:A:S" or flag == 4):

					continue

				if flag == 2048 or flag == 2064:

					continue

			# ---- sequence / quality
			if include_primary_seq:

				sequence, phred = fields[9], fields[10]

			else:

				# '*' is the valid SAM placeholder (an empty field is rejected by htslib)
				sequence, phred = "*", "*"

			alignment_id += 1

			partial_1_row = (
				alignment_id,
				fields[0], fields[1], fields[2], fields[3], fields[4],
				fields[5], fields[6], fields[7], fields[8],
				fields[11], fields[12], fields[13], tp_tag,
			)

			partial_2_row = (alignment_id, sequence, phred)

			yield partial_1_row, partial_2_row


# =============================================================================
# Main
# =============================================================================

def main():

	if len(sys.argv) < 4:

		sys.exit(__doc__)

	prefix = sys.argv[1]
	bam = sys.argv[2]
	tax_map_tsv = sys.argv[3]

	# same semantics as the bash tests: only the literal string 'false' turns these off
	single_end_arg = sys.argv[4] if len(sys.argv) > 4 else ""
	include_seq_arg = sys.argv[5] if len(sys.argv) > 5 else ""

	paired_end = single_end_arg == "false"
	single_end = not paired_end
	include_primary_seq = include_seq_arg != "false"

	sam_db = f"{prefix}-sam.db"

	# deleting DB if it exists already (for reruns in work dir)
	try:
		os.remove(sam_db)
	except FileNotFoundError:
		pass

	# -------------------------------------------------------------------------
	# If a SAM file was used as input, sort it into BAM format first
	# -------------------------------------------------------------------------
	if bam.endswith(".sam"):

		sorted_bam = f"{prefix}.bam"

		pysam.sort("-@", str(SORT_THREADS), "-o", sorted_bam, bam)

		bam = sorted_bam

	# -------------------------------------------------------------------------
	# Database connection
	# -------------------------------------------------------------------------

	# autocommit mode; the bulk insert below opens its own explicit transaction
	conn = sqlite3.connect(sam_db, isolation_level=None)

	conn.execute("PRAGMA cache_size = -20000000")
	conn.execute("PRAGMA temp_store = FILE") 
	conn.execute("PRAGMA synchronous = OFF")
	# conn.execute("PRAGMA mmap_size = 50000000000")
	conn.execute("PRAGMA journal_mode = OFF") 

	# -------------------------------------------------------------------------
	# Create tax_map table and import tax map TSV file
	# -------------------------------------------------------------------------

	log("tax map import START")

	conn.execute("""
		DROP TABLE IF EXISTS tax_map;
	""")

	conn.execute("""
		CREATE TABLE tax_map (seq_id TEXT PRIMARY KEY, tax_id TEXT, name TEXT) WITHOUT ROWID;
	""")

	with open(tax_map_tsv, newline="", encoding="utf-8") as f:

		reader = csv.reader(f, delimiter="\t")

		conn.execute("BEGIN")
		conn.executemany("INSERT INTO tax_map VALUES (?, ?, ?)", reader)
		conn.execute("COMMIT")

	log("tax map import END")

	# -------------------------------------------------------------------------
	# Create sam_partial_1 / sam_partial_2 and import the trimmed alignments
	# -------------------------------------------------------------------------
	conn.execute("""
		DROP TABLE IF EXISTS sam_partial_1;
	""")

	conn.execute("""
		DROP TABLE IF EXISTS sam_partial_2;
	""")

	conn.execute("""
		CREATE TABLE sam_partial_1 (alignment_id INTEGER PRIMARY KEY, read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next INTEGER, t_length TEXT, nm_tag TEXT, ms_tag TEXT, as_tag INTEGER, tp_tag TEXT);
	""")

	conn.execute("""
		CREATE TABLE sam_partial_2 (alignment_id INTEGER PRIMARY KEY, sequence TEXT, phred BLOB);
	""")

	print("Converting BAM to SAM with no headers", flush=True)

	log("SAM import START")

	rows = trimmed_sam_rows(bam, single_end, include_primary_seq)
	total_rows = 0

	conn.execute("BEGIN")

	while True:

		batch = list(itertools.islice(rows, BATCH_SIZE))

		if not batch:

			break

		partial_1_rows, partial_2_rows = zip(*batch)

		conn.executemany("INSERT INTO sam_partial_1 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", partial_1_rows)
		conn.executemany("INSERT INTO sam_partial_2 VALUES (?,?,?)", partial_2_rows)

		total_rows += len(batch)

	conn.execute("COMMIT")

	log(f"SAM import FINISH ({total_rows} alignments)")

	# -------------------------------------------------------------------------
	# Remove accession trailing version number (e.g. NC_000001.1 -> NC_000001)
	# -------------------------------------------------------------------------
	# log("Remove accession trailing version num START")

	# conn.executescript("""
	# 	ALTER TABLE sam_partial_1 ADD COLUMN ref_id_alt TEXT;

	# 	UPDATE sam_partial_1
	# 	SET ref_id_alt = CASE
	# 		WHEN INSTR(ref_id, '.') > 0 THEN SUBSTR(ref_id, 1, INSTR(ref_id, '.') - 1)
	# 		ELSE ref_id
	# 	END;
	# """)

	# log("Remove accession trailing version num END")

	# -------------------------------------------------------------------------
	# Add tax_id / name from tax_map to sam_partial_1
	# -------------------------------------------------------------------------
	log("Adding tax_id to sam_partial_1 table START")

	conn.execute("""
		DROP TABLE IF EXISTS sam_partial_1_new;
	""")

	conn.executescript("""
		CREATE TABLE sam_partial_1_new AS
		SELECT sp.alignment_id, sp.read_id, sp.flag, sp.ref_id, sp.left_most_position, sp.mapq, sp.cigar, sp.r_next, sp.p_next, sp.t_length, sp.nm_tag, sp.ms_tag, sp.as_tag, sp.tp_tag, tax_map.tax_id, tax_map.name
		FROM sam_partial_1 AS sp
		LEFT JOIN tax_map ON (sp.ref_id = tax_map.seq_id);

		DROP TABLE sam_partial_1;

		ALTER TABLE sam_partial_1_new RENAME TO sam_partial_1;

		DROP TABLE tax_map;
	""")

	log("Adding tax_id to sam_partial_1 table END")

	# -------------------------------------------------------------------------
	# Label unmapped reads as 'unclassified'
	# -------------------------------------------------------------------------
	if paired_end:

		# paired reads which are both unmapped have their ref_id = *
		# a mate that is unmapped when the other is mapped will have r_next = '=' and cigar = *
		conn.executescript("""
			UPDATE sam_partial_1
			SET ref_id = 'unclassified', tax_id = 'unclassified', name = 'unclassified'
			WHERE ref_id = '*' or (r_next = '=' and cigar = '*');
		""")

	else:

		conn.executescript("""
			UPDATE sam_partial_1
			SET ref_id = 'unclassified', tax_id = 'unclassified', name = 'unclassified'
			WHERE ref_id = '*' and flag = 4;
		""")

	# -------------------------------------------------------------------------
	# Fall back to ref_id when no tax_id was found
	# -------------------------------------------------------------------------
	log("Assigning ref_id to tax_id if tax_id value is NULL START")

	conn.executescript("""
		UPDATE sam_partial_1
		SET tax_id = ref_id, name = ref_id
		WHERE tax_id IS NULL;
	""")

	log("Assigning ref_id to tax_id if tax_id value is NULL END")

	log("Cleaning up databse START")
	# clearing out intermediate tables and reclaiming unused space
	conn.execute("PRAGMA temp_store = MEMORY")
	conn.execute("VACUUM;")

	log("Cleaning up databse END")

	conn.close()


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":

	start = time.monotonic()

	try:

		main()

	finally:

		print(f"Total run time: {format_duration(time.monotonic() - start)}", flush=True)