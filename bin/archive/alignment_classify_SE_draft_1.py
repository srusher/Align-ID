#!/usr/bin/env python3
"""Python port of the alignment_classify_SE.sh bash script (primary + secondary
alignment-score version).

Usage: alignment_classify_SE.py PREFIX BAM SAM_DB [TAXA_NAMES] [MAPQ]

TAXA_NAMES and MAPQ are accepted for compatibility but unused (as in the original).
sqlite3 (stdlib) replaces the sqlite container and pysam (bundles samtools/htslib)
replaces the samtools container, so no Singularity/Apptainer is needed.
Requires: pysam.
"""

import glob
import os
import sqlite3
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


def log(msg, stamp=True):

	print(f"{msg}: {now()}" if stamp else msg, flush=True)


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


def sam_to_sorted_bam(sam, sorted_bam):
	"""Equivalent of `samtools view -b SAM | samtools sort` (sort reads SAM directly)."""

	if "unmapped" in os.path.basename(sam):

		# silences 'unrecognized reference name "unclassified"' warnings
		pysam.sort("--verbosity", "0", "-@", str(THREADS), "-o", sorted_bam, sam)

	else:

		pysam.sort("-@", str(THREADS), "-o", sorted_bam, sam)


# =============================================================================
# Main
# =============================================================================

def main():

	if len(sys.argv) < 4:

		sys.exit(__doc__)

	prefix, bam, sam_db = sys.argv[1:4]

	# defining all alignment output files
	summary_file = f"{prefix}-taxonomic-summary.tsv"
	primary_all = f"{prefix}-primary_all.sam"
	primary_unambiguous = f"{prefix}-primary_unambiguous.sam"
	primary_amb_single = f"{prefix}-primary_ambiguous_single_genome.sam"
	primary_amb_multi = f"{prefix}-primary_ambiguous_multi_genome.sam"
	unmapped = f"{prefix}-unmapped.sam"

	# -------------------------------------------------------------------------
	# BAM headers
	# -------------------------------------------------------------------------
	log("Parsing BAM headers", stamp=False)

	header = bam_header_text(bam)

	for sam in (primary_all, primary_unambiguous, primary_amb_single, primary_amb_multi, unmapped):

		with open(sam, "w") as fh:

			fh.write(header)

	# -------------------------------------------------------------------------
	# Database connection
	# -------------------------------------------------------------------------

	# autocommit mode so each statement is committed immediately
	conn = sqlite3.connect(sam_db, timeout=5, isolation_level=None)

	conn.execute("PRAGMA cache_size = -50000000")
	conn.execute("PRAGMA temp_store = MEMORY")

	def run(sql, fetch=False):
		"""Run SQL, retrying while the database is locked. Returns rows if fetch."""

		while True:

			try:

				if fetch:

					return conn.execute(sql).fetchall()

				conn.executescript(sql)

				return None

			except sqlite3.OperationalError as e:

				if "locked" in str(e):

					print(f"{sam_db} is locked - trying again in 10 seconds", flush=True)
					time.sleep(10)

				else:

					raise

	def stream(sql):
		"""Like run(fetch=True) but yields rows lazily (for large alignment dumps)."""

		while True:

			try:

				cur = conn.execute(sql)
				break

			except sqlite3.OperationalError as e:

				if "locked" in str(e):

					print(f"{sam_db} is locked - trying again in 10 seconds", flush=True)
					time.sleep(10)

				else:

					raise

		yield from cur

	def fmt(v):

		return "" if v is None else str(v)

	# #########################################################################
	# SAM parsing section
	# #########################################################################

	# -------------------------------------------------------------------------
	# total_reads
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS total_reads;

		CREATE TABLE total_reads (num_reads FLOAT);
	""")

	log("Calculating total reads")

	run("""
		INSERT INTO total_reads
		SELECT COUNT(read_id)
		FROM sam_complete
		WHERE flag = 4 or flag = 16 or flag = 0;
	""")

	# -------------------------------------------------------------------------
	# unmapped
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS unmapped;

		CREATE TABLE unmapped (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER);
	""")

	log("Populating unmapped table")

	run("""
		INSERT INTO unmapped
		SELECT read_id, flag, ref_id, mapq
		FROM sam_complete
		WHERE ref_id = 'unclassified';
	""")

	# -------------------------------------------------------------------------
	# primary_multimap
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS primary_multimap;

		CREATE TABLE primary_multimap (read_id TEXT, flag INTEGER, ref_id TEXT, tax_id TEXT, name TEXT, mapq INTEGER, as_tag TEXT);
	""")

	log("Populating primary_multimap table")

	run("""
		INSERT INTO primary_multimap
		SELECT read_id, flag, ref_id, tax_id, name, mapq, as_tag
		FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and mapq = 0;
	""")

	# -------------------------------------------------------------------------
	# secondary_mappings
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS secondary_mappings;

		CREATE TABLE secondary_mappings (sec_read_id TEXT, sec_flag INTEGER, sec_ref_id TEXT, sec_tax_id TEXT, sec_name TEXT, sec_mapq INTEGER, sec_as_tag TEXT);
	""")

	log("Populating secondary_mappings table")

	run("""
		INSERT INTO secondary_mappings
		SELECT read_id, flag, ref_id, tax_id, name, mapq, as_tag
		FROM sam_complete
		WHERE tp_tag = 'tp:A:S' and read_id IN (SELECT read_id FROM primary_multimap);
	""")

	# -------------------------------------------------------------------------
	# primary_secondary_aggregate
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS primary_secondary_aggregate;

		CREATE TABLE primary_secondary_aggregate (read_id TEXT, sec_tax_id INTEGER, sec_name TEXT, sec_tax_id_count INTEGER, prim_tax_id INTEGER, prim_name TEXT);
	""")

	log("Populating primary_secondary_aggregate table")

	run("""
		INSERT INTO primary_secondary_aggregate
		SELECT secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id, secondary_mappings.sec_name, COUNT(secondary_mappings.sec_tax_id), primary_multimap.tax_id, primary_multimap.name
		FROM secondary_mappings
		JOIN primary_multimap ON (secondary_mappings.sec_read_id = primary_multimap.read_id)
		WHERE secondary_mappings.sec_as_tag = primary_multimap.as_tag or cast(substr(secondary_mappings.sec_as_tag, 6) AS INTEGER) > cast(substr(primary_multimap.as_tag, 6) AS INTEGER)
		GROUP BY secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id;
	""")

	# -------------------------------------------------------------------------
	# primary_unambiguous
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS primary_unambiguous;

		CREATE TABLE primary_unambiguous (read_id TEXT, flag INTEGER, ref_id TEXT, tax_id TEXT, name TEXT, mapq INTEGER, as_tag TEXT);
	""")

	log("Populating primary_unambiguous table")

	# a MAPQ of 0 doesn't guarantee a secondary alignment with a >= alignment score (AS) as the
	# primary alignment, so to grab all unambiguous primary alignments we use the tp_tag and make
	# sure the read_id does not appear in the primary_secondary_aggregate table
	run("""
		INSERT INTO primary_unambiguous
		SELECT read_id, flag, ref_id, tax_id, name, mapq, as_tag
		FROM sam_complete
		WHERE tp_tag = 'tp:A:P' AND read_id NOT IN (
			SELECT read_id
			FROM primary_secondary_aggregate
		);
	""")

	# -------------------------------------------------------------------------
	# single genome hits
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS single_genome_hits_by_read;

		CREATE TABLE single_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_name TEXT, sec_tax_id_count INTEGER, prim_tax_id INTEGER, prim_name TEXT);
	""")

	log("Populating single_genome_hits_by_read table")

	run("""
		INSERT INTO single_genome_hits_by_read
		SELECT read_id, sec_tax_id, sec_name, sec_tax_id_count, prim_tax_id, prim_name
		FROM primary_secondary_aggregate
		WHERE sec_tax_id = prim_tax_id
		GROUP BY read_id, prim_tax_id
		HAVING COUNT(read_id) = 1;
	""")

	run("""
		DROP TABLE IF EXISTS single_genome_hits_by_name;

		CREATE TABLE single_genome_hits_by_name (prim_name TEXT, total_hits);
	""")

	log("Populating single_genome_hits_by_name table")

	run("""
		INSERT INTO single_genome_hits_by_name
		SELECT prim_name, COUNT(prim_name)
		FROM single_genome_hits_by_read
		GROUP BY prim_name;
	""")

	# -------------------------------------------------------------------------
	# multi genome hits
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS multi_genome_hits_by_read;

		CREATE TABLE multi_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_name TEXT, sec_tax_id_count INTEGER, prim_tax_id INTEGER, prim_name TEXT);
	""")

	log("Populating multi_genome_hits_by_read table")

	run("""
		INSERT INTO multi_genome_hits_by_read
		SELECT read_id, sec_tax_id, sec_name, sec_tax_id_count, prim_tax_id, prim_name
		FROM primary_secondary_aggregate
		WHERE read_id NOT IN (
			SELECT read_id
			FROM single_genome_hits_by_read
		);
	""")

	run("""
		DROP TABLE IF EXISTS multi_genome_hits_by_name;

		CREATE TABLE multi_genome_hits_by_name (prim_name TEXT, total_hits);
	""")

	log("Populating multi_genome_hits_by_name table")

	run("""
		INSERT INTO multi_genome_hits_by_name
		SELECT prim_name, COUNT(DISTINCT read_id)
		FROM multi_genome_hits_by_read
		GROUP BY prim_name;
	""")

	# -------------------------------------------------------------------------
	# unambiguous hits
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS unambiguous_hits_by_name;

		CREATE TABLE unambiguous_hits_by_name (prim_name TEXT, total_hits);
	""")

	log("Populating unambiguous_hits_by_name table")

	run("""
		INSERT INTO unambiguous_hits_by_name
		SELECT name, COUNT(DISTINCT read_id)
		FROM primary_unambiguous
		GROUP BY name;
	""")

	# #########################################################################
	# Summary file
	# #########################################################################
	log("Populating Summary File")

	summary_queries = [

		("none", """
			SELECT unmapped.ref_id, COUNT(unmapped.read_id),
				   ROUND(COUNT(unmapped.read_id) / total_reads.num_reads * 100, 6) as percentage
			FROM unmapped
			LEFT JOIN total_reads
			GROUP BY ref_id
		"""),

		("none", """
			SELECT unambiguous_hits_by_name.prim_name, unambiguous_hits_by_name.total_hits,
				   ROUND(unambiguous_hits_by_name.total_hits / total_reads.num_reads * 100, 6) as percentage
			FROM unambiguous_hits_by_name
			LEFT JOIN total_reads
			GROUP BY prim_name
		"""),

		("single_genome", """
			SELECT single_genome_hits_by_name.prim_name, single_genome_hits_by_name.total_hits,
				   ROUND(single_genome_hits_by_name.total_hits / total_reads.num_reads * 100, 6) as percentage
			FROM single_genome_hits_by_name
			LEFT JOIN total_reads
			GROUP BY prim_name
		"""),

		("multi_genome", """
			SELECT multi_genome_hits_by_name.prim_name, multi_genome_hits_by_name.total_hits,
				   ROUND(multi_genome_hits_by_name.total_hits / total_reads.num_reads * 100, 6) as percentage
			FROM multi_genome_hits_by_name
			LEFT JOIN total_reads
			GROUP BY prim_name
		"""),

	]

	rows = []  # (species, num_reads, percent, ambiguity)

	for ambiguity, sql in summary_queries:

		for r in run(sql, fetch=True):

			rows.append([fmt(r[0]), fmt(r[1]), fmt(r[2]), ambiguity])

	# order by species total (desc), then per-row read count (desc), then full line
	totals = {}

	for sp, n, _, _ in rows:

		totals[sp] = totals.get(sp, 0) + float(n or 0)

	rows.sort(key=lambda r: (-totals[r[0]], -float(r[1] or 0), "\t".join(r)))

	with open(summary_file, "w") as fh:

		fh.write("species\tnum_reads\tpercent_reads\tambiguity\n")

		for r in rows:

			fh.write("\t".join(r) + "\n")

	# #########################################################################
	# Alignment files
	# #########################################################################

	def dump_alignments(sql, out_path):
		"""Append first 15 columns of each row; drop empty optional fields (col 12+)."""

		with open(out_path, "a") as fh:

			for row in stream(sql):

				cols = [fmt(v) for v in row[:15]]
				cols += [""] * (15 - len(cols))

				cols = cols[:11] + [c for c in cols[11:] if c != ""]

				fh.write("\t".join(cols) + "\n")

	log("Populating primary_all alignment file")

	dump_alignments("SELECT * FROM sam_complete WHERE tp_tag = 'tp:A:P'", primary_all)

	log("Populating primary_unambiguous alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and read_id IN (SELECT read_id FROM primary_unambiguous)
	""", primary_unambiguous)

	log("Populating primary_ambig_single alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id IN (SELECT read_id FROM single_genome_hits_by_read)
	""", primary_amb_single)

	log("Populating primary_ambig_multi alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id NOT IN (SELECT read_id FROM single_genome_hits_by_read)
	""", primary_amb_multi)

	log("Populating unmapped alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE read_id IN (SELECT read_id FROM unmapped)
	""", unmapped)

	# clearing out intermediate tables and reclaiming unused space
	run("VACUUM;")

	conn.close()

	# #########################################################################
	# BAM output
	# #########################################################################
	log("Formatting all primary output BAM file")

	sam_to_sorted_bam(primary_all, f"{prefix}-primary_all_sorted.bam")

	log("Formatting primary unambiguous single genome output BAM file")

	sam_to_sorted_bam(primary_unambiguous, f"{prefix}-primary_unambiguous_sorted.bam")

	log("Formatting primary ambiguous single genome output BAM file")

	sam_to_sorted_bam(primary_amb_single, f"{prefix}-primary_ambiguous_single_genome_sorted.bam")

	log("Formatting primary ambiguous multi genome output BAM file")

	sam_to_sorted_bam(primary_amb_multi, f"{prefix}-primary_ambiguous_multi_genome_sorted.bam")

	log("Formatting Unmapped output BAM file")

	sam_to_sorted_bam(unmapped, f"{prefix}-unmapped-sorted.bam")

	# #########################################################################
	# Cleanup
	# #########################################################################
	for pattern in (f"{prefix}-chunk_*", f"{prefix}-part-*", "sam_chunk*", "*.sam", f"{prefix}-ambig.*"):

		for path in glob.glob(pattern):

			try:

				os.remove(path)

			except OSError:

				pass


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":

	start = time.monotonic()

	try:

		main()

	finally:

		print(f"Total run time: {format_duration(time.monotonic() - start)}", flush=True)