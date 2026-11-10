#!/bin/bash
SAMTOOLS_CONTAINER="https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
seqid2taxid=$3
single_end=$4
include_primary_seq=$5

cur_dir="$(pwd)"

# sqlite_pragmas="
# PRAGMA cache_size = -5000000000;
# PRAGMA temp_store_directory = '$cur_dir';
# PRAGMA journal_mode = OFF;
# PRAGMA synchronous = OFF;
# "

sqlite_pragmas="
PRAGMA cache_size = -50000000;
PRAGMA temp_store = MEMORY;
"

# if a sam file was used as input, we need to convert it to bam format first
if [[ $bam == *".sam" ]]; then

    singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort -@ 4 -o "$prefix".bam $bam
    bam="$prefix".bam

fi

# singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam > "$prefix-complete.sam" #converting bam to sam for easier parsing in the loop below

# complete_sam=$prefix-complete.sam

sam_db="$prefix-sam.db"
trimmed_sam="$prefix-trimmed.sam"

echo "Tax DB copy START: $(date '+%Y-%m-%d %H:%M:%S')"
# creating a copy of the original tax_map DB file to preserve structure for furture runs - prevents us from having to import, index, and dedup tax_map table each run
cp $seqid2taxid $cur_dir/tax_map_copy.db

mv $cur_dir/tax_map_copy.db $sam_db

echo "Tax DB copy END: $(date '+%Y-%m-%d %H:%M:%S')"

echo "Converting BAM to SAM with no headers"

if [[ $include_primary_seq == 'false' ]]; then

	awk_script='{print $1, $2, $3, $4, $5, $6, $7, $8, $9, "", "", $12, $13, $14, $16}'

else

	awk_script='{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}'

fi


if [[ $single_end == "false" ]]; then

	# trimming SAM file for entries we want, the 'sed' line will escape any rogue double quotes ["] in the phred score - if we don't do this, SQLite won't be able to import the SAM file
	singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam | awk '$2 < 2000' | sed 's|"|\\"|g' | awk -v OFS='\t' "$awk_script" > $trimmed_sam

else

	# trimming SAM file for entries we want, the 'sed' line will escape any rogue double quotes ["] in the phred score - if we don't do this, SQLite won't be able to import the SAM file
	singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam | awk '$16 == "tp:A:P" || $16 == "tp:A:S" || $2 == 4' | awk '$2 != 2048 && $2 != 2064' | sed 's|"|\\"|g' | awk -v OFS='\t' "$awk_script" > $trimmed_sam

fi


singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE sam_complete (read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next INTEGER, t_length TEXT, sequence TEXT, phred BLOB, nm_tag TEXT, ms_tag TEXT, as_tag INTEGER, tp_tag TEXT);
	"

echo "SAM import START: $(date '+%Y-%m-%d %H:%M:%S')"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
$sqlite_pragmas
.mode tabs
.import $trimmed_sam sam_complete
EOF

echo "SAM import FINISH: $(date '+%Y-%m-%d %H:%M:%S')"

echo "Remove escape characters START: $(date '+%Y-%m-%d %H:%M:%S')"

# Getting rid of extra escape character that was added to double quotes to allow SQLite import
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
$sql_pragmas
UPDATE sam_complete
SET phred = replace(phred, '\\\"', '\"');
"

echo "Remove escape characters END: $(date '+%Y-%m-%d %H:%M:%S')"

echo "Remove accession trailing version num START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	$sqlite_pragmas
# 	BEGIN TRANSACTION;
# 	UPDATE sam_complete
# 	SET ref_id = SUBSTR(ref_id, 1, INSTR(ref_id, '.') - 1)
# 	WHERE INSTR(ref_id, '.') > 0;
# 	COMMIT;
# 	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
ALTER TABLE sam_complete ADD COLUMN ref_id_alt TEXT;
UPDATE sam_complete
SET ref_id_alt = CASE
    WHEN INSTR(ref_id, '.') > 0 THEN SUBSTR(ref_id, 1, INSTR(ref_id, '.') - 1)
    ELSE ref_id
END;
"
echo "Remove accession trailing version num END: $(date '+%Y-%m-%d %H:%M:%S')"


# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	CREATE TABLE tax_map (seq_id TEXT, tax_id TEXT);
# 	"

# echo "Taxonomy import START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
# $sqlite_pragmas
# .mode tabs
# .import $seqid2taxid tax_map
# EOF

# echo "Taxonomy import END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "DB indexing START: $(date '+%Y-%m-%d %H:%M:%S')"
# # creating indicies for faster queries
# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	$sqlite_pragmas
# 	CREATE INDEX IF NOT EXISTS idx_sam_read_id   ON sam_complete(read_id);
# 	CREATE INDEX IF NOT EXISTS idx_sam_ref_id    ON sam_complete(ref_id);
# 	CREATE INDEX IF NOT EXISTS idx_sam_tp_tag    ON sam_complete(tp_tag);
# 	CREATE INDEX IF NOT EXISTS idx_sam_flag      ON sam_complete(flag);
# 	CREATE INDEX IF NOT EXISTS idx_sam_complete_filter ON sam_complete(r_next, cigar);
# 	"

# echo "DB indexing FINISH: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "Manual ref seq ID add START: $(date '+%Y-%m-%d %H:%M:%S')"
# # manually adding reference sequence IDs to taxonomy map if they had not already been added
# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	$sqlite_pragmas
# 	INSERT INTO tax_map
# 	SELECT ref_id, ref_id
# 	FROM sam_complete
# 	WHERE ref_id NOT IN (
# 		SELECT seq_id
# 		FROM tax_map
# 	);

# 	"
# echo "Manual ref seq ID add END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "Creating index on tax_map(seq_id) START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	$sqlite_pragmas
# 	CREATE INDEX IF NOT EXISTS idx_tax_map_seq_id ON tax_map(seq_id);
# 	"

# echo "Creating index on tax_map(seq_id) END: $(date '+%Y-%m-%d %H:%M:%S')"

echo "Adding tax_id to sam_complete table START: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	$sqlite_pragmas
	CREATE TABLE sam_complete_new AS
	SELECT sc.read_id, sc.flag, sc.ref_id, sc.left_most_position, sc.mapq, sc.cigar, sc.r_next, sc.p_next, sc.t_length, sc.sequence, sc.phred, sc.nm_tag, sc.ms_tag, sc.as_tag, sc.tp_tag, tax_map.tax_id, tax_map.name
	FROM sam_complete AS sc
	LEFT JOIN tax_map ON (sc.ref_id_alt = tax_map.seq_id);

	DROP TABLE sam_complete;

	ALTER TABLE sam_complete_new RENAME TO sam_complete;

	DROP TABLE tax_map;
	"

echo "Adding tax_id to sam_complete table END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "DB indexing round 2 START: $(date '+%Y-%m-%d %H:%M:%S')"
# # creating indicies for any sequences that didn't appearin the original tax_map
# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	$sqlite_pragmas
# 	CREATE INDEX IF NOT EXISTS idx_taxmap_seqid  ON tax_map(seq_id);
# 	CREATE INDEX IF NOT EXISTS idx_taxmap_taxid  ON tax_map(tax_id);
# 	"

# echo "DB indexing round 2 END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "seq ID dedup START: $(date '+%Y-%m-%d %H:%M:%S')"
# # removing duplicate rows from seq id to tax id map file
# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
# 	"
# 	$sqlite_pragmas
# 	DELETE FROM tax_map
# 	WHERE EXISTS (
# 		SELECT 1
# 		FROM tax_map t2
# 		WHERE t2.seq_id = tax_map.seq_id
# 		AND t2.rowid < tax_map.rowid
# 	);
# 	"
# echo "seq ID dedup END: $(date '+%Y-%m-%d %H:%M:%S')"

# labeling unmapped reads as 'unclassified'
if [[ $single_end == "false" ]]; then
	
	# paired reads which are both unmapped have their ref_id = *
	# a mate that is unmapped when the other is mapped will have it's r_next = '=' and cigar = *
	singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
		"
		$sqlite_pragmas
		UPDATE sam_complete
		SET ref_id = 'unclassified', tax_id = 'unclassified', name = 'unclassified'
		WHERE ref_id = '*' or (r_next = '=' and cigar = '*');
		"

else

	singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
		"
		$sqlite_pragmas
		UPDATE sam_complete
		SET ref_id = 'unclassified', tax_id = 'unclassified', name = 'unclassified'
		WHERE ref_id = '*' and flag = 4;
		"

fi

echo "Assigning ref_id to tax_id if tax_id value is NULL START: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
	$sql_pragmas
	UPDATE sam_complete
	SET tax_id = ref_id, name = ref_id
	WHERE tax_id IS NULL;

"

echo "Assigning ref_id to tax_id if tax_id value is NULL END: $(date '+%Y-%m-%d %H:%M:%S')"
# cleaning up temporary files
rm -f *.sam