#!/usr/bin/env python3
"""Python port of the BAM/SAM -> SQLite import bash script.

Usage: import_sam_to_db.py PREFIX BAM_OR_SAM SEQID2TAXID_TSV SINGLE_END INCLUDE_PRIMARY_SEQ

  PREFIX			  prefix for output files; the database is written to PREFIX-sam.db
  BAM_OR_SAM		  input alignment file (a .sam input is sorted into PREFIX.bam first)
  SEQID2TAXID_TSV	  TSV with three columns: seq_id, tax_id, name
  SINGLE_END		  'false' for paired-end data, anything else means single-end
  INCLUDE_PRIMARY_SEQ 'false' to drop sequence / quality strings, anything else keeps them

Tables written to PREFIX-sam.db:
  sam_partial_1	everything except sequence and phred (plus alignment_id, tax_id, name)
  sam_partial_2	alignment_id, sequence, phred (linked to sam_partial_1 by alignment_id)

Requires: pysam.
"""

import csv
import itertools
import os
import shutil
import sqlite3
import sys
import time
from datetime import datetime

import pysam

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


def load_tax_map(path):
	"""seq_id -> (tax_id, name). The first occurrence of a seq_id wins."""

	tax_map = {}

	with open(path, newline="", encoding="utf-8") as f:

		for seq_id, tax_id, name in csv.reader(f, delimiter="\t"):

			tax_map.setdefault(seq_id, (tax_id, name))

	return tax_map


def trimmed_sam_rows(bam, single_end, include_primary_seq, tax_map):
	"""Equivalent of the `samtools view | awk | sed | awk` trimming pipeline, plus the
	tax_map lookup and 'unclassified' labeling (formerly done with SQL after the import).

	Yields one (partial_1_row, partial_2_row) pair per kept alignment:

	partial_1_row (16 values):
		alignment_id, read_id, flag, ref_id, pos, mapq, cigar, r_next, p_next,
		t_length, nm_tag, ms_tag, as_tag, tp_tag, tax_id, name

	partial_2_row (3 values):
		alignment_id, sequence, phred
	"""

	paired_end = not single_end
	alignment_id = 0

	with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=SORT_THREADS) as alignment_file:

		for alignment in alignment_file:

			flag = alignment.flag

			# ---- cheap flag-only filters first, before building any strings
			if paired_end:

				if flag >= 2000:

					continue

			elif flag == 2048 or flag == 2064:

				continue

			fields = alignment.to_string().split("\t")

			# awk returns an empty string for fields that don't exist
			if len(fields) < 16:

				fields += [""] * (16 - len(fields))

			tp_tag = fields[15]

			# ---- single-end tp tag filter
			if single_end and not (tp_tag == "tp:A:P" or tp_tag == "tp:A:S" or flag == 4):

				continue

			# ---- tax_id / name lookup and 'unclassified' labeling
			ref_id = fields[2]

			if paired_end:

				# both mates unmapped: ref_id = '*'; a mate that is unmapped while the
				# other is mapped has r_next = '=' and cigar = '*'
				unclassified = ref_id == "*" or (fields[6] == "=" and fields[5] == "*")

			else:

				unclassified = ref_id == "*" and flag == 4

			if unclassified:

				ref_id = tax_id = name = "unclassified"

			else:

				hit = tax_map.get(ref_id)

				# no tax_map entry: fall back to the ref_id itself
				tax_id, name = hit if hit else (ref_id, ref_id)

			# ---- sequence / quality
			if include_primary_seq:

				sequence, phred = fields[9], fields[10]

			else:

				# '*' is the valid SAM placeholder (an empty field is rejected by htslib)
				sequence, phred = "*", "*"

			alignment_id += 1

			partial_1_row = (
				alignment_id,
				fields[0], fields[1], ref_id, fields[3], fields[4],
				fields[5], fields[6], fields[7], fields[8],
				fields[11], fields[12], fields[13], tp_tag,
				tax_id, name,
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

	# deleting DB (and any leftover side files) if it exists already (for reruns in work dir)
	for suffix in ("", "-journal", "-wal", "-shm"):

		try:

			os.remove(sam_db + suffix)

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
	# Load the tax map into memory (about 50k entries, so tiny)
	# -------------------------------------------------------------------------
	log("tax map load START")

	tax_map = load_tax_map(tax_map_tsv)

	log(f"tax map load END ({len(tax_map)} entries)")

	# -------------------------------------------------------------------------
	# Database connection
	# -------------------------------------------------------------------------

	# autocommit mode; the bulk insert below opens its own explicit transaction
	conn = sqlite3.connect(sam_db, isolation_level=None)

	# page_size only takes effect on a new, empty database, so set it first.
	# Larger pages suit wide rows and sequential writes.
	conn.execute("PRAGMA page_size = 32768")

	# appends in alignment_id order need little cache; there are no joins or sorts anymore
	conn.execute("PRAGMA cache_size = -2000000")
	conn.execute("PRAGMA temp_store = FILE")
	conn.execute("PRAGMA synchronous = OFF")
	conn.execute("PRAGMA journal_mode = OFF")

	# -------------------------------------------------------------------------
	# Create sam_partial_1 / sam_partial_2 and import the labeled alignments
	# -------------------------------------------------------------------------
	conn.execute("""
		CREATE TABLE sam_partial_1 (alignment_id INTEGER PRIMARY KEY, read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next INTEGER, t_length TEXT, nm_tag TEXT, ms_tag TEXT, as_tag INTEGER, tp_tag TEXT, tax_id TEXT, name TEXT);
	""")

	conn.execute("""
		CREATE TABLE sam_partial_2 (alignment_id INTEGER PRIMARY KEY, sequence TEXT, phred BLOB);
	""")

	print("Converting BAM to SAM with no headers", flush=True)

	log("SAM import START")

	rows = trimmed_sam_rows(bam, single_end, include_primary_seq, tax_map)
	total_rows = 0

	conn.execute("BEGIN")

	while True:

		batch = list(itertools.islice(rows, BATCH_SIZE))

		if not batch:

			break

		partial_1_rows, partial_2_rows = zip(*batch)

		conn.executemany("INSERT INTO sam_partial_1 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", partial_1_rows)
		conn.executemany("INSERT INTO sam_partial_2 VALUES (?,?,?)", partial_2_rows)

		total_rows += len(batch)

		# log(f"{total_rows} alignments imported")

	conn.execute("COMMIT")

	log(f"SAM import FINISH ({total_rows} alignments)")

	# no join copy, no UPDATE passes, and nothing was dropped, so there are no free
	# pages to reclaim and no VACUUM is needed

	conn.close()


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":

	start = time.monotonic()

	try:

		main()

	finally:

		shutil.rmtree(tmp_dir, ignore_errors=True)

		print(f"Total run time: {format_duration(time.monotonic() - start)}", flush=True)