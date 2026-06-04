#!/bin/bash
SAMTOOLS_CONTAINER="/scicomp/groups-pure/WDPB/EMEL/singularity/samtools/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="/scicomp/groups-pure/WDPB/EMEL/singularity/sqlite3/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
seqid2taxid=$3
my_tax_ids=$4
include_children=$5
nodes_sqlite=$6
sam_db=$7

tax_ids_i_want="$prefix-parent_and_child_ids.txt"
>$tax_ids_i_want

# Determining child nodes based on parent tax ID provided
if [[ "$filter_alignment_by_id" == "true" && "$include_children" == "true" ]]; then

    echo "Collecting all children tax ids from SQL database"

    while IFS= read -r line || [[ -n "$line" ]]; do

        data=$(singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $nodes_sqlite "SELECT * FROM TAX_IDS WHERE parent_id = "$line"")

        parent_id="$line"        
        echo $parent_id >> $tax_ids_i_want

        if [[ -n $data ]]; then

            child_ids=$(echo $data | cut -d '|' -f3 | sed 's/,/ /g')

            for x in $child_ids; do

                echo $x >> $tax_ids_i_want
            
            done

        fi

    done < $my_tax_ids

else
    # if we do not want to include child taxonomy nodes then we simply use the tax ids provided
	while IFS= read -r line || [[ -n "$line" ]]; do
    
		echo $line >> $tax_ids_i_want

	done < $my_tax_ids

fi

# defining all alignment output files
primary_all="$prefix-primary_all.sam"

echo "Parsing BAM headers"

# printing bam headers to all primary alignment output sam files
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_all

echo "Converting BAM to SAM with no headers"

singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam > "$prefix-complete.sam" #converting bam to sam for easier parsing in the loop below

complete_sam=$prefix-complete.sam

# clearing past SQL table entries
singularity -q exec --bind /scicomp $SQLITE3_CONTAINER sqlite3 $sam_db \
"
DROP TABLE IF EXISTS tax_ids_i_want;
"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
"
CREATE TABLE tax_ids_i_want (tax_id TEXT);
"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.import $tax_ids_i_want tax_ids_i_want
EOF

# grabbing all alignments whos tax_id are found within the list of specified tax ids
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
"
SELECT *
FROM sam_complete
WHERE tp_tag = 'tp:A:P' and read_id IN (
    SELECT read_id
    FROM sam_complete
    JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
    WHERE tp_tag = 'tp:A:P' and tax_id IN (SELECT tax_id FROM tax_ids_i_want)
)
" | sed 's/!##########!/\t/g' >> $primary_all

echo Converting SAM back to BAM

singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_all > $prefix-primary-all.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-all.bam > $prefix-primary-all-filtered_sorted.bam
rm -f $prefix-primary-all.bam