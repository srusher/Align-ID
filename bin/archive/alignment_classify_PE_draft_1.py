#!/usr/bin/env python3
"""Python port of alignment_classify_PE.sh bash script.

Usage: alignment_classify_PE.py PREFIX BAM SAM_DB [TAXA_NAMES] [MAPQ]

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
	# total_read_pairs
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS total_read_pairs;

		CREATE TABLE total_read_pairs (num_reads FLOAT);
	""")

	log("Calculating total read pairs")

	run("""
		INSERT INTO total_read_pairs
		SELECT COUNT(DISTINCT read_id)
		FROM sam_complete;
	""")

	# -------------------------------------------------------------------------
	# unmapped_pairs
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS unmapped_pairs;

		CREATE TABLE unmapped_pairs (read_id TEXT, r1_flag INTEGER, r2_flag INTEGER, ref_id TEXT, mapq INTEGER);
	""")

	log("Populating unmapped pairs table")

	run("""
		INSERT INTO unmapped_pairs
		SELECT read_id, MIN(flag), MAX(flag), ref_id, mapq
		FROM sam_complete
		WHERE ref_id = 'unclassified' and (flag = 77 or flag = 141)
		GROUP BY read_id;
	""")

	# -------------------------------------------------------------------------
	# r1_reads
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS r1_reads;

		CREATE TABLE r1_reads (read_id TEXT, r1_flag INTEGER, r1_ref_id TEXT, r1_tax_id TEXT, r1_name TEXT, r1_mapq INTEGER, r1_p_next INTEGER, r1_as_tag INTEGER);
	""")

	log("Populating r1-reads table")

	run("""
		INSERT INTO r1_reads
		SELECT read_id, MIN(flag), ref_id, tax_id, name, mapq, p_next, as_tag
		FROM sam_complete
		WHERE tp_tag = 'tp:A:P' or (r_next = '=' and cigar = '*')
		GROUP BY read_id;
	""")

	# -------------------------------------------------------------------------
	# r2_reads
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS r2_reads;

		CREATE TABLE r2_reads (read_id TEXT, r2_flag INTEGER, r2_ref_id TEXT, r2_tax_id TEXT, r2_name TEXT, r2_mapq INTEGER, r2_p_next INTEGER, r2_as_tag INTEGER);
	""")

	log("Populating r2-reads table")

	run("""
		INSERT INTO r2_reads
		SELECT read_id, MAX(flag), ref_id, tax_id, name, mapq, p_next, as_tag
		FROM sam_complete
		WHERE tp_tag = 'tp:A:P' or (r_next = '=' and cigar = '*')
		GROUP BY read_id;
	""")

	# -------------------------------------------------------------------------
	# mapped_pairs
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS mapped_pairs;

		CREATE TABLE mapped_pairs (read_id TEXT, r1_flag INTEGER, r2_flag INTEGER, r1_ref_id TEXT, r2_ref_id TEXT, r1_tax_id TEXT, r1_name TEXT, r2_tax_id TEXT, r2_name TEXT, r1_mapq INTEGER, r2_mapq INTEGER, r1_p_next INTEGER, r2_p_next INTEGER, r1_as_tag INTEGER, r2_as_tag INTEGER);
	""")

	log("Populating mapped_pairs table")

	run("""
		INSERT INTO mapped_pairs
		SELECT r1.read_id, r1.r1_flag, r2.r2_flag, r1.r1_ref_id, r2.r2_ref_id, r1.r1_tax_id, r1.r1_name, r2.r2_tax_id, r2.r2_name, r1.r1_mapq, r2.r2_mapq, r1.r1_p_next, r2.r2_p_next, r1.r1_as_tag, r2.r2_as_tag
		FROM r1_reads as r1
		JOIN r2_reads as r2 ON (r1.read_id = r2.read_id)
		GROUP BY r1.read_id;
	""")

	log("Updating mapped_pairs table unclassified tax IDs")

	# one mate unmapped while the other is mapped -> convert to unclassified
	run("""
		UPDATE mapped_pairs
		SET r1_tax_id = 'unclassified', r1_name = 'unclassified'
		WHERE r1_ref_id = 'unclassified' and r1_tax_id IS NULL;

		UPDATE mapped_pairs
		SET r2_tax_id = 'unclassified', r2_name = 'unclassified'
		WHERE r2_ref_id = 'unclassified' and r2_tax_id IS NULL;
	""")

	# -------------------------------------------------------------------------
	# secondary mappings
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS top_secondary_mappings;

		CREATE TABLE top_secondary_mappings (read_id TEXT, sec_tax_id TEXT, sec_name TEXT, sec_p_next INTEGER, sec_as_tag INTEGER);
	""")

	log("Populating top_secondary_mappings table")

	run("""
		INSERT INTO top_secondary_mappings
		SELECT read_id, tax_id as tax_id, name, p_next, MAX(as_tag) as as_tag
		FROM sam_complete
		WHERE tp_tag = 'tp:A:S'
		GROUP BY read_id, p_next, tax_id;
	""")

	run("""
		DROP TABLE IF EXISTS top_secondary_mappings_taxa_count;

		CREATE TABLE top_secondary_mappings_taxa_count (read_id TEXT, taxa_count INTEGER);
	""")

	log("Populating top_secondary_mappings_taxa_count table")

	run("""
		INSERT INTO top_secondary_mappings_taxa_count
		SELECT read_id, COUNT(DISTINCT sec_tax_id) as taxa_count
		FROM top_secondary_mappings
		GROUP BY read_id;
	""")

	run("""
		DROP TABLE IF EXISTS secondary_mapping_aggregate;

		CREATE TABLE secondary_mapping_aggregate (read_id TEXT, sec_tax_id TEXT, sec_name TEXT, sec_p_next INTEGER, sec_as_tag INTEGER, sec_taxa_count INTEGER);
	""")

	log("Populating secondary_mapping_aggregate table")

	run("""
		INSERT INTO secondary_mapping_aggregate
		SELECT top_secondary_mappings.read_id, top_secondary_mappings.sec_tax_id, top_secondary_mappings.sec_name, top_secondary_mappings.sec_p_next, top_secondary_mappings.sec_as_tag, top_secondary_mappings_taxa_count.taxa_count
		FROM top_secondary_mappings
		LEFT JOIN top_secondary_mappings_taxa_count ON (top_secondary_mappings.read_id = top_secondary_mappings_taxa_count.read_id);
	""")

	# -------------------------------------------------------------------------
	# unambiguous_pairs
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS unambiguous_pairs;

		CREATE TABLE unambiguous_pairs (read_id TEXT, r1_flag INTEGER, r2_flag INTEGER, r1_ref_id TEXT, r2_ref_id TEXT, r1_tax_id TEXT, r1_name TEXT, r2_tax_id TEXT, r2_name TEXT, r1_mapq INTEGER, r2_mapq INTEGER, r1_p_next INTEGER, r2_p_next INTEGER, r1_as_tag INTEGER, r2_as_tag INTEGER);
	""")

	log("Populating unambiguous_pairs table")

	run("""
		INSERT INTO unambiguous_pairs
		SELECT *
		FROM mapped_pairs
		WHERE r1_tax_id = r2_tax_id and r1_mapq > 0 and r2_mapq > 0;
	""")

	# -------------------------------------------------------------------------
	# ambiguous_single_genome_pairs
	# -------------------------------------------------------------------------
	log("Populating ambiguous_single_genome_pairs table")

	run("""
		DROP TABLE IF EXISTS ambiguous_single_genome_pairs;

		CREATE TABLE ambiguous_single_genome_pairs (read_id TEXT, r1_tax_id TEXT, r1_name TEXT, r2_tax_id, r2_name TEXT, r1_mapq INTEGER, r2_mapq INTEGER, sec_tax_id TEXT, sec_taxa_count INTEGER);
	""")

	run("""
		INSERT INTO ambiguous_single_genome_pairs
		SELECT mp.read_id, mp.r1_tax_id, mp.r1_name, mp.r2_tax_id, mp.r2_name, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
		FROM mapped_pairs as mp
		JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r1_tax_id = sma.sec_tax_id
		WHERE mp.r1_tax_id = mp.r2_tax_id and (mp.r1_mapq = 0 or mp.r2_mapq = 0) and sma.sec_taxa_count = 1
		GROUP BY mp.read_id

		UNION ALL

		SELECT mp.read_id, mp.r1_tax_id, mp.r1_name, mp.r2_tax_id, mp.r2_name, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
		FROM mapped_pairs as mp
		JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r1_tax_id = sma.sec_tax_id
		WHERE mp.r1_tax_id != 'unclassified' and mp.r2_tax_id = 'unclassified' and mp.r1_mapq = 0 and sma.sec_taxa_count = 1
		GROUP BY mp.read_id

		UNION ALL

		SELECT mp.read_id, mp.r1_tax_id, mp.r1_name, mp.r2_tax_id, mp.r2_name, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
		FROM mapped_pairs as mp
		JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r2_tax_id = sma.sec_tax_id
		WHERE mp.r1_tax_id = 'unclassified' and mp.r2_tax_id != 'unclassified' and mp.r2_mapq = 0 and sma.sec_taxa_count = 1
		GROUP BY mp.read_id

		UNION ALL

		SELECT mp.read_id, mp.r1_tax_id, mp.r1_name, mp.r2_tax_id, mp.r2_name, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
		FROM mapped_pairs as mp
		JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r2_tax_id = sma.sec_tax_id
		WHERE mp.r1_tax_id != 'unclassified' and mp.r2_tax_id = 'unclassified' and mp.r1_mapq > 0
		GROUP BY mp.read_id

		UNION ALL

		SELECT mp.read_id, mp.r1_tax_id, mp.r1_name, mp.r2_tax_id, mp.r2_name, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
		FROM mapped_pairs as mp
		JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r2_tax_id = sma.sec_tax_id
		WHERE mp.r1_tax_id = 'unclassified' and mp.r2_tax_id != 'unclassified' and mp.r2_mapq > 0
		GROUP BY mp.read_id;
	""")

	# -------------------------------------------------------------------------
	# ambiguous_multi_genome_pairs
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS ambiguous_multi_genome_pairs;

		CREATE TABLE ambiguous_multi_genome_pairs (read_id TEXT, r1_tax_id TEXT, r1_name TEXT, r2_tax_id TEXT, r2_name TEXT, r1_mapq INTEGER, r2_mapq INTEGER, r1_as_tag INTEGER, r2_as_tag INTEGER);
	""")

	log("Populating ambiguous_multi_genome_pairs table")

	run("""
		INSERT INTO ambiguous_multi_genome_pairs
		SELECT mp.read_id, mp.r1_tax_id, mp.r1_name, mp.r2_tax_id, mp.r2_name, mp.r1_mapq, mp.r2_mapq, mp.r1_as_tag, mp.r2_as_tag
		FROM mapped_pairs as mp
		WHERE mp.read_id NOT IN (
			SELECT read_id FROM unambiguous_pairs

			UNION ALL

			SELECT read_id FROM ambiguous_single_genome_pairs
		);
	""")

	# -------------------------------------------------------------------------
	# hits by name
	# -------------------------------------------------------------------------
	run("""
		DROP TABLE IF EXISTS unambiguous_hits_by_name;

		CREATE TABLE unambiguous_hits_by_name (name TEXT, total_hits);
	""")

	log("Populating unambiguous_hits_by_name table")

	run("""
		INSERT INTO unambiguous_hits_by_name
		SELECT r1_name, COUNT(DISTINCT read_id)
		FROM unambiguous_pairs
		GROUP BY r1_name;
	""")

	run("""
		DROP TABLE IF EXISTS ambiguous_single_genome_hits_by_name;

		CREATE TABLE ambiguous_single_genome_hits_by_name (name TEXT, total_hits);
	""")

	log("Populating ambiguous_single_genome_hits_by_name table")

	# primary tax ID = the mate that is NOT 'unclassified' when applicable
	run("""
		INSERT INTO ambiguous_single_genome_hits_by_name
		SELECT q1.name, COUNT(q1.name)
		FROM (
			SELECT CASE
				WHEN r1_name = 'unclassified' THEN r2_name
				WHEN r2_name = 'unclassified' THEN r1_name
				ELSE r1_name
			END name
			FROM ambiguous_single_genome_pairs
		) as q1
		GROUP BY q1.name;
	""")

	run("""
		DROP TABLE IF EXISTS ambiguous_multi_genome_hits_by_name;

		CREATE TABLE ambiguous_multi_genome_hits_by_name (name TEXT, total_hits);
	""")

	log("Populating ambiguous_multi_genome_hits_by_name table")

	run("""
		INSERT INTO ambiguous_multi_genome_hits_by_name
		SELECT q1.name, COUNT(q1.name)
		FROM (
			SELECT CASE
				WHEN r1_name = 'unclassified' THEN r2_name
				WHEN r2_name = 'unclassified' THEN r1_name
				WHEN r1_as_tag > r2_as_tag THEN r1_name
				WHEN r2_as_tag > r1_as_tag THEN r2_name
				ELSE r1_name
			END name
			FROM ambiguous_multi_genome_pairs
		) as q1
		GROUP BY q1.name;
	""")

	# #########################################################################
	# Summary file
	# #########################################################################
	log("Populating Summary File")

	summary_queries = [

		("none", """
			SELECT unmapped_pairs.ref_id, COUNT(unmapped_pairs.read_id),
				   ROUND(COUNT(unmapped_pairs.read_id) / total_read_pairs.num_reads * 100, 6) as percentage
			FROM unmapped_pairs
			LEFT JOIN total_read_pairs
			GROUP BY ref_id
		"""),

		("none", """
			SELECT unambiguous_hits_by_name.name, unambiguous_hits_by_name.total_hits,
				   ROUND(unambiguous_hits_by_name.total_hits / total_read_pairs.num_reads * 100, 6) as percentage
			FROM unambiguous_hits_by_name
			LEFT JOIN total_read_pairs
			GROUP BY name
		"""),

		("single_genome", """
			SELECT ambiguous_single_genome_hits_by_name.name, ambiguous_single_genome_hits_by_name.total_hits,
				   ROUND(ambiguous_single_genome_hits_by_name.total_hits / total_read_pairs.num_reads * 100, 6) as percentage
			FROM ambiguous_single_genome_hits_by_name
			LEFT JOIN total_read_pairs
			GROUP BY name
		"""),

		("multi_genome", """
			SELECT ambiguous_multi_genome_hits_by_name.name, ambiguous_multi_genome_hits_by_name.total_hits,
				   ROUND(ambiguous_multi_genome_hits_by_name.total_hits / total_read_pairs.num_reads * 100, 6) as percentage
			FROM ambiguous_multi_genome_hits_by_name
			LEFT JOIN total_read_pairs
			GROUP BY name
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

	def dump_alignments(sql, out_path, strip_empty_aux=True):
		"""Append first 15 columns of each row; drop empty optional fields (col 12+)."""

		with open(out_path, "a") as fh:

			for row in stream(sql):

				cols = [fmt(v) for v in row[:15]]
				cols += [""] * (15 - len(cols))

				if strip_empty_aux:

					cols = cols[:11] + [c for c in cols[11:] if c != ""]

				fh.write("\t".join(cols) + "\n")

	log("Populating primary_all alignment file")

	dump_alignments("SELECT * FROM sam_complete WHERE tp_tag = 'tp:A:P'", primary_all)

	log("Populating primary_unambiguous alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and read_id IN (SELECT read_id FROM unambiguous_pairs)
	""", primary_unambiguous)

	log("Populating primary_ambig_single alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and read_id IN (SELECT read_id FROM ambiguous_single_genome_pairs)
	""", primary_amb_single)

	log("Populating primary_ambig_multi alignment file")

	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE tp_tag = 'tp:A:P' and read_id IN (SELECT read_id FROM ambiguous_multi_genome_pairs)
	""", primary_amb_multi)

	log("Populating unmapped alignment file")

	# empty optional fields are dropped (blank fields are invalid SAM and htslib rejects them)
	dump_alignments("""
		SELECT * FROM sam_complete
		WHERE read_id IN (SELECT read_id FROM unmapped_pairs)
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