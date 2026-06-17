#!/bin/bash
SAMTOOLS_CONTAINER="https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
seqid2taxid=$3
single_end=$4

# if a sam file was used as input, we need to convert it to bam format first
if [[ $bam == *".sam" ]]; then

    singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort -@ 4 -o "$prefix".bam $bam
    bam="$prefix".bam

fi

echo "Converting BAM to SAM with no headers"

singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam > "$prefix-complete.sam" #converting bam to sam for easier parsing in the loop below

complete_sam=$prefix-complete.sam

sam_db="$prefix-sam.db"
trimmed_sam="$prefix-trimmed.sam"

rm -f $sam_db

if [[ $single_end == "false" ]]; then

	# trimming SAM file for entries we want, the 'sed' line will escape any rogue double quotes ["] in the phred score - if we don't do this, SQLite won't be able to import the SAM file
	singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam | awk '$2 < 2000' | sed 's|"|\\"|g' | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' > $trimmed_sam

else

	# trimming SAM file for entries we want, the 'sed' line will escape any rogue double quotes ["] in the phred score - if we don't do this, SQLite won't be able to import the SAM file
	singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam | awk '$16 == "tp:A:P" || $16 == "tp:A:S" || $2 == 4' | awk '$2 != 2048 && $2 != 2064' | sed 's|"|\\"|g' | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' > $trimmed_sam

fi


singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE sam_complete (read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next INTEGER, t_length TEXT, sequence TEXT, phred BLOB, nm_tag TEXT, ms_tag TEXT, as_tag TEXT, tp_tag TEXT);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.mode tabs
.import $trimmed_sam sam_complete
EOF

# Getting rid of extra escape character that was added to double quotes to allow SQLite import
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
UPDATE sam_complete
SET phred = replace(phred, '\\\"', '\"')
"

if [[ $single_end == "false" ]]; then
	
	# paired reads which are both unmapped have their ref_id = *
	# a mate that is unmapped when the other is mapped will have it's r_next = '=' and cigar = *
	singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
		"
		UPDATE sam_complete
		SET ref_id = 'unclassified'
		WHERE ref_id = '*' or (r_next = '=' and cigar = '*')
		"

else

	singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
		"
		UPDATE sam_complete
		SET ref_id = 'unclassified'
		WHERE ref_id = '*' and flag = 4
		"

fi


singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE tax_map (seq_id TEXT, tax_id TEXT);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.mode tabs
.import $seqid2taxid tax_map
EOF

# manually adding reference sequence IDs to taxonomy map if they had not already been added
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO tax_map
	SELECT ref_id, ref_id
	FROM sam_complete
	WHERE ref_id NOT IN (
		SELECT seq_id
		FROM tax_map
	)
	"

# removing duplicate rows from seq id to tax id map file
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	DELETE FROM tax_map
	WHERE rowid NOT IN (
		SELECT MIN(rowid)
		FROM tax_map
		GROUP BY seq_id	
	)
	"

# cleaning up temporary files
rm -f *.sam