#!/bin/bash
SAMTOOLS_CONTAINER="https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"
PYTHON_CONTAINER="https://depot.galaxyproject.org/singularity/python%3A3.12.12"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
my_tax_ids=$3
include_children=$4
nodes=$5
python_script=$6
sam_db=$7

if [[ $include_children == 'true' || $include_children == 'True' ]]; then

    sed 's/\t|\t/\t/g' $nodes | awk -F'\t' -v OFS='\t' '{print $1,$2}' > nodes_mod.dmp

    singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python $python_script nodes_mod.dmp $my_tax_ids

    tax_ids_i_want="tax_id_child_nodes.tsv"

else

    tax_ids_i_want=$my_tax_ids

fi

sqlite_pragmas="
PRAGMA cache_size = -50000000;
PRAGMA temp_store = MEMORY;
"

# # Determining child nodes based on parent tax ID provided
# if [[ "$filter_alignment_by_id" == "true" && "$include_children" == "true" ]]; then

#     echo "Collecting all children tax ids from SQL database"

#     while IFS= read -r line || [[ -n "$line" ]]; do

#         data=$(singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $nodes_sqlite "SELECT * FROM TAX_IDS WHERE parent_id = "$line"")

#         parent_id="$line"        
#         echo $parent_id >> $tax_ids_i_want

#         if [[ -n $data ]]; then

#             child_ids=$(echo $data | cut -d '|' -f3 | sed 's/,/ /g')

#             for x in $child_ids; do

#                 echo $x >> $tax_ids_i_want
            
#             done

#         fi

#     done < $my_tax_ids

# else
#     # if we do not want to include child taxonomy nodes then we simply use the tax ids provided
# 	while IFS= read -r line || [[ -n "$line" ]]; do
    
# 		echo $line >> $tax_ids_i_want

# 	done < $my_tax_ids

# fi

# defining all alignment output files
primary_all="$prefix-primary_all.sam"

echo "Parsing BAM headers"

# printing bam headers to all primary alignment output sam files
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_all

# clearing past SQL table entries
singularity -q exec --bind /scicomp $SQLITE3_CONTAINER sqlite3 $sam_db \
"
DROP TABLE IF EXISTS tax_ids_i_want;
CREATE TABLE tax_ids_i_want (tax_id TEXT);
"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.import $tax_ids_i_want tax_ids_i_want
EOF

# singularity -q exec --bind /scicomp $SQLITE3_CONTAINER sqlite3 $sam_db \
# "
# CREATE INDEX tax_id_idx ON tax_ids_i_want (tax_id);
# "



# count pre-filter read-pairs
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
$sql_pragmas
SELECT COUNT(DISTINCT read_id)
FROM sam_complete;
" > $prefix'_read_pairs_count.tsv'

echo "Removing all SAM entires from DB that didn't match tax ID filter START: $(date '+%Y-%m-%d %H:%M:%S')"
# grabbing all alignments whos tax_id are found within the list of specified tax ids
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".mode tabs" \
"
$sql_pragmas
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
" 

echo "Removing all SAM entires from DB that didn't match tax ID filter END: $(date '+%Y-%m-%d %H:%M:%S')"

# count post-filter read-pairs
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
$sql_pragmas
SELECT COUNT(DISTINCT read_id)
FROM sam_complete;
" >> $prefix'_read_pairs_count.tsv'

echo "Printing DB entries to SAM file START: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".mode tabs" \
"
$sql_pragmas
SELECT *
FROM sam_complete
" | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' | awk -F'\t' -v OFS='\t' '
{
    line = $1
    for (i = 2; i <= 11; i++) {
        line = line OFS $i
    }
    for (i = 12; i <= NF; i++) {
        if ($i != "") {
            line = line OFS $i
        }
    }
    print line
}' >> $primary_all

echo "Printing DB entries to SAM file END: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".mode tabs" \
# "
# $sql_pragmas
# SELECT *
# FROM sam_complete
# " | head -2

# awk -F'\t' -v OFS='\t' '
# {
#     line = $1
#     for (i = 2; i <= 11; i++) {
#         line = line OFS $i
#     }
#     for (i = 12; i <= NF; i++) {
#         if ($i != "") {
#             line = line OFS $i
#         }
#     }
#     print line
# }' $primary_all > $prefix-primary-all-fixed.sam

echo Converting SAM back to BAM

singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_all > $prefix-primary-all.bam

echo Sorting BAM

singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-all.bam > $prefix-primary-all-filtered_sorted.bam

rm -f $prefix-primary-all.bam
rm -f *.sam