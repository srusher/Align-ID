#!/usr/bin/env python3
"""Python port of the BAM/SAM -> SQLite import bash script.

Usage: import_sam_to_db.py PREFIX BAM_OR_SAM SEQID2TAXID_DB SINGLE_END INCLUDE_PRIMARY_SEQ

  PREFIX			  prefix for output files; the database is written to PREFIX-sam.db
  BAM_OR_SAM		  input alignment file (a .sam input is sorted into PREFIX.bam first)
  SEQID2TAXID_DB	  prebuilt SQLite database containing a tax_map(seq_id, tax_id, name)
					  table; it is copied to PREFIX-sam.db and left untouched
  SINGLE_END		  'false' for paired-end data, anything else means single-end
  INCLUDE_PRIMARY_SEQ 'false' to drop sequence / quality strings, anything else keeps them

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
import itertools, time

import pysam


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

	Yields one 15-tuple per kept alignment:
		read_id, flag, ref_id, pos, mapq, cigar, r_next, p_next, t_length,
		sequence, phred, nm_tag, ms_tag, as_tag, tp_tag

	The optional fields are positional, exactly as in the awk script
	($12, $13, $14 and $16). The double-quote escaping done by `sed` is not needed
	because rows are inserted with bound parameters instead of `.import`.
	"""

	with pysam.AlignmentFile(bam, "rb", check_sq=False) as alignment_file:

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

			yield (
				fields[0], fields[1], fields[2], fields[3], fields[4],
				fields[5], fields[6], fields[7], fields[8],
				sequence, phred,
				fields[11], fields[12], fields[13], tp_tag,
			)


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

	conn.execute("PRAGMA cache_size = -5000000")
	conn.execute("PRAGMA temp_store = MEMORY")

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
	# Create sam_complete and import the trimmed alignments
	# -------------------------------------------------------------------------
	conn.execute("""
		DROP TABLE IF EXISTS sam_complete;
	""")
			  
	conn.execute("""
		CREATE TABLE sam_complete (read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next INTEGER, t_length TEXT, sequence TEXT, phred BLOB, nm_tag TEXT, ms_tag TEXT, as_tag INTEGER, tp_tag TEXT);
	""")

	print("Converting BAM to SAM with no headers", flush=True)

	log("SAM import START")

	rows = trimmed_sam_rows(bam, single_end, include_primary_seq)
	total_rows = 0

	# t = time.monotonic()
	# n = sum(1 for _ in itertools.islice(trimmed_sam_rows(bam, single_end, include_primary_seq), 1_000_000))
	# print(n, time.monotonic() - t)

	conn.execute("BEGIN")

	while True:

		batch = list(itertools.islice(rows, BATCH_SIZE))

		if not batch:

			break

		conn.executemany("INSERT INTO sam_complete VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", batch)

	conn.execute("COMMIT")

	log(f"SAM import FINISH ({total_rows} alignments)")

	# -------------------------------------------------------------------------
	# Remove accession trailing version number (e.g. NC_000001.1 -> NC_000001)
	# -------------------------------------------------------------------------
	# log("Remove accession trailing version num START")

	# conn.executescript("""
	# 	ALTER TABLE sam_complete ADD COLUMN ref_id_alt TEXT;

	# 	UPDATE sam_complete
	# 	SET ref_id_alt = CASE
	# 		WHEN INSTR(ref_id, '.') > 0 THEN SUBSTR(ref_id, 1, INSTR(ref_id, '.') - 1)
	# 		ELSE ref_id
	# 	END;
	# """)

	# log("Remove accession trailing version num END")

	# -------------------------------------------------------------------------
	# Add tax_id / name from tax_map to sam_complete
	# -------------------------------------------------------------------------
	log("Adding tax_id to sam_complete table START")

	conn.execute("""
		DROP TABLE IF EXISTS sam_complete_new;
	""")

	conn.executescript("""
		CREATE TABLE sam_complete_new AS
		SELECT sc.read_id, sc.flag, sc.ref_id, sc.left_most_position, sc.mapq, sc.cigar, sc.r_next, sc.p_next, sc.t_length, sc.sequence, sc.phred, sc.nm_tag, sc.ms_tag, sc.as_tag, sc.tp_tag, tax_map.tax_id, tax_map.name
		FROM sam_complete AS sc
		LEFT JOIN tax_map ON (sc.ref_id = tax_map.seq_id);

		DROP TABLE sam_complete;

		ALTER TABLE sam_complete_new RENAME TO sam_complete;

		DROP TABLE tax_map;
	""")

	log("Adding tax_id to sam_complete table END")

	# -------------------------------------------------------------------------
	# Label unmapped reads as 'unclassified'
	# -------------------------------------------------------------------------
	if paired_end:

		# paired reads which are both unmapped have their ref_id = *
		# a mate that is unmapped when the other is mapped will have r_next = '=' and cigar = *
		conn.executescript("""
			UPDATE sam_complete
			SET ref_id = 'unclassified', tax_id = 'unclassified', name = 'unclassified'
			WHERE ref_id = '*' or (r_next = '=' and cigar = '*');
		""")

	else:

		conn.executescript("""
			UPDATE sam_complete
			SET ref_id = 'unclassified', tax_id = 'unclassified', name = 'unclassified'
			WHERE ref_id = '*' and flag = 4;
		""")

	# -------------------------------------------------------------------------
	# Fall back to ref_id when no tax_id was found
	# -------------------------------------------------------------------------
	log("Assigning ref_id to tax_id if tax_id value is NULL START")

	conn.executescript("""
		UPDATE sam_complete
		SET tax_id = ref_id, name = ref_id
		WHERE tax_id IS NULL;
	""")

	log("Assigning ref_id to tax_id if tax_id value is NULL END")

	# clearing out intermediate tables and reclaiming unused space
	conn.execute("VACUUM;")

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