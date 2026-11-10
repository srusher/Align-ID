#!/bin/bash
SAMTOOLS_CONTAINER="https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
sam_db=$3
taxa_names=$4
mapq=$5

# defining all alignment output files
summary_file="$prefix-taxonomic-summary.tsv"
primary_all="$prefix-primary_all.sam"
primary_unambiguous="$prefix-primary_unambiguous.sam"
primary_ambiguous_single_genome="$prefix-primary_ambiguous_single_genome.sam"
primary_ambiguous_multi_genome="$prefix-primary_ambiguous_multi_genome.sam"

cur_dir="$(pwd)"

sqlite_pragmas="
PRAGMA cache_size = -50000000;
PRAGMA temp_store = MEMORY;
"

echo "Parsing BAM headers"

# printing bam headers to all primary alignment output sam files
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_all 
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_unambiguous
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_single_genome
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_multi_genome
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $unmapped

# this loop will ensure that the SAM database isn't locked from the previous the import module
run_with_retry() {

    local sql="$1"
    local output
    local loop_again='true'

    while [[ "$loop_again" == 'true' ]]; do

        output=$(singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".mode tabs" "$sql" 2>&1)
        
        if [[ "$output" == *'database is locked'* ]]; then

            echo "$sam_db is locked - trying again in 10 seconds"
            sleep 10

        else

            loop_again='false'

        fi
    done

	if [[ $output != '' ]]; then

		echo "$output"
	
	fi

}

#############################################################################################
###                                 SAM Parsing Section                                  ####
#############################################################################################


run_with_retry \
"
DROP TABLE IF EXISTS total_read_pairs; 
CREATE TABLE total_read_pairs (num_reads FLOAT);
"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO total_reads
	SELECT COUNT(read_id)
	FROM sam_complete
	WHERE flag = 4 or flag = 16 or flag = 0;
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS unmapped;
	CREATE TABLE unmapped (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO unmapped
	SELECT read_id, flag, ref_id, mapq
	FROM sam_complete
	WHERE ref_id = 'unclassified'
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS primary_multimap;
	CREATE TABLE primary_multimap (read_id TEXT, flag INTEGER, ref_id TEXT, tax_id TEXT, name TEXT, mapq INTEGER, as_tag TEXT);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO primary_multimap
	SELECT read_id, flag, ref_id, tax_id, name, mapq, as_tag
	FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS secondary_mappings;
	CREATE TABLE secondary_mappings (sec_read_id TEXT, sec_flag INTEGER, sec_ref_id TEXT, sec_tax_id TEXT, sec_name TEXT, sec_mapq INTEGER, sec_as_tag TEXT);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO secondary_mappings
	SELECT read_id, flag, ref_id, tax_id, name, mapq, as_tag
	FROM sam_complete
	WHERE tp_tag = 'tp:A:S' and read_id IN (SELECT read_id FROM primary_multimap)
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS primary_secondary_aggregate;
	CREATE TABLE primary_secondary_aggregate (read_id TEXT, sec_tax_id INTEGER, sec_name TEXT, sec_tax_id_count INTEGER, prim_tax_id INTEGER, prim_name TEXT);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO primary_secondary_aggregate
	SELECT secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id, secondary_mappings.sec_name, COUNT(secondary_mappings.sec_tax_id), primary_multimap.tax_id, primary_multimap.name
	FROM secondary_mappings
	JOIN primary_multimap ON (secondary_mappings.sec_read_id = primary_multimap.read_id) 
	WHERE secondary_mappings.sec_as_tag = primary_multimap.as_tag or cast(substr(secondary_mappings.sec_as_tag, 6) AS INTEGER) > cast(substr(primary_multimap.as_tag, 6) AS INTEGER)
	GROUP BY secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id
	" 

run_with_retry \
	"
	DROP TABLE IF EXISTS primary_unambiguous;
	CREATE TABLE primary_unambiguous (read_id TEXT, flag INTEGER, ref_id TEXT, tax_id TEXT, name TEXT, mapq INTEGER, as_tag TEXT);
	"

# had to change the condition for WHERE clause - a MAPQ of 0 doesn't guarantee a secondary alignment with a >= alignment score (AS) as the primary alignment meaning that the read may not be accounted for in the ambiguous (single or multi genome) queries; so in order to grab all unambiguous primary alignments we will use the tp_tag and make sure the read_id does not appear in the primary_secondary_aggregate table
run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO primary_unambiguous
	SELECT read_id, flag, ref_id, tax_id, name, mapq, as_tag
	FROM sam_complete
	WHERE tp_tag = 'tp:A:P' AND read_id NOT IN (
		SELECT read_id
		FROM primary_secondary_aggregate
	)
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS single_genome_hits_by_read;
	CREATE TABLE single_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_name TEXT, sec_tax_id_count INTEGER, prim_tax_id INTEGER, prim_name TEXT);
	"	

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO single_genome_hits_by_read
	SELECT read_id, sec_tax_id, sec_name, sec_tax_id_count, prim_tax_id, prim_name
	FROM primary_secondary_aggregate
	WHERE sec_tax_id = prim_tax_id 
	GROUP BY read_id, prim_tax_id
	HAVING COUNT(read_id) = 1
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS single_genome_hits_by_name;
	CREATE TABLE single_genome_hits_by_name (prim_name TEXT, total_hits);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO single_genome_hits_by_name
	SELECT prim_name, COUNT(prim_name)
	FROM single_genome_hits_by_read
	GROUP BY prim_name
	" 

run_with_retry \
	"
	DROP TABLE IF EXISTS multi_genome_hits_by_read;
	CREATE TABLE multi_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_name TEXT, sec_tax_id_count INTEGER, prim_tax_id INTEGER, prim_name TEXT);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO multi_genome_hits_by_read
	SELECT read_id, sec_tax_id, sec_name, sec_tax_id_count, prim_tax_id, prim_name
	FROM primary_secondary_aggregate
	WHERE read_id NOT IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	) 
	" 

run_with_retry \
	"
	DROP TABLE IF EXISTS multi_genome_hits_by_name;
	CREATE TABLE multi_genome_hits_by_name (prim_name TEXT, total_hits);
	"

run_with_retry \
	"
	$sqlite_pragmas
	INSERT INTO multi_genome_hits_by_name
	SELECT prim_name, COUNT(DISTINCT read_id)
	FROM multi_genome_hits_by_read
	GROUP BY prim_name
	"

run_with_retry \
	"
	DROP TABLE IF EXISTS unambiguous_hits_by_name;
	CREATE TABLE unambiguous_hits_by_name (prim_name TEXT, total_hits);
	"

run_with_retry \
	"
	$sqlite_pragmas
    INSERT INTO unambiguous_hits_by_name
	SELECT name, COUNT(DISTINCT read_id)
	FROM primary_unambiguous
	GROUP BY name;
	"

echo -e "species\tnum_reads\tpercent_reads\tambiguity" > "$summary_file"

run_with_retry \
	"
	$sqlite_pragmas
    SELECT unmapped.ref_id, COUNT(unmapped.read_id), ROUND(COUNT(unmapped.read_id) / total_reads.num_reads * 100,6) as percentage
    FROM unmapped
	LEFT JOIN total_reads
	GROUP BY ref_id
	" | awk -F'\t' -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

run_with_retry \
	"
	$sqlite_pragmas
    SELECT unambiguous_hits_by_name.prim_name, unambiguous_hits_by_name.total_hits, ROUND(unambiguous_hits_by_name.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM unambiguous_hits_by_name
	LEFT JOIN total_reads
	GROUP BY prim_name
	" | awk -F'\t' -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

run_with_retry \
	"
	$sqlite_pragmas
    SELECT single_genome_hits_by_name.name, single_genome_hits_by_name.total_hits, ROUND(single_genome_hits_by_name.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM single_genome_hits_by_name
	LEFT JOIN total_reads
	GROUP BY prim_name
	" | awk -F'\t' -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

run_with_retry \
	"
	$sqlite_pragmas
    SELECT multi_genome_hits_by_name.name, multi_genome_hits_by_name.total_hits, ROUND(multi_genome_hits_by_name.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM multi_genome_hits_by_name
	LEFT JOIN total_reads
	GROUP BY name
	" | awk -F'\t' -v OFS='\t' '{print $0, "none"}' >> "$summary_file"


# reordering summary file to list classifications from highest read count to lowest and to group by species
{

  head -n 1 $summary_file
  tail -n +2 $summary_file | awk -F'\t' '
    { total[$1]+=$2 }
    { rows[NR]=$0; sp[NR]=$1 }
    END {
      for (i=1;i<=NR;i++) print total[sp[i]]"\t"rows[i]
    }' | sort -t$'\t' -k1,1nr -k3,3nr | cut -f2-

} > $summary_file".sorted"

rm $summary_file
mv $summary_file".sorted" $summary_file

# grab all primary alignments
run_with_retry \
	"
	$sqlite_pragmas
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P'
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

# grab all nonambiguous primary alignments
run_with_retry \
	"
	$sqlite_pragmas
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM primary_unambiguous
	)
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
	}' >> $primary_unambiguous

# grab all primary ambiguous single genome hits
run_with_retry \
	"
	$sqlite_pragmas
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	)
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
	}' >> $primary_ambiguous_single_genome

# grab all primary ambiguous multi-genome hits
run_with_retry \
	"
	$sqlite_pragmas
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id NOT IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	) 
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
	}' >> $primary_ambiguous_multi_genome

# grab all unmapped reads - unmapped reads will often contain empty aux fields so we need to remove trailing tab characters from the SQL query
run_with_retry \
	"
	$sqlite_pragmas
    SELECT *
    FROM sam_complete
	WHERE read_id IN (
		SELECT read_id
		FROM unmapped
	) 
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
	}' >> $unmapped

#############################################################################################
#############################################################################################
#############################################################################################


echo "Formatting all primary output BAM file: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view --verbosity 0 -@ 16 -b $primary_all > $prefix-primary-all.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-all.bam > $prefix-primary_all_sorted.bam
rm -f $prefix-primary-all.bam

echo "Formatting primary unambiguous single genome output BAM file"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view --verbosity 0 -@ 16 -b $primary_unambiguous > $prefix-primary-unambiguous.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-unambiguous.bam > $prefix-primary_unambiguous_sorted.bam
rm -f $prefix-primary-unambiguous.bam

echo "Formatting primary ambiguous single genome output BAM file: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view --verbosity 0 -@ 16 -b $primary_ambiguous_single_genome > $prefix-primary_ambiguous_single_genome.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary_ambiguous_single_genome.bam > $prefix-primary_ambiguous_single_genome_sorted.bam
rm -f $prefix-primary_ambiguous_single_genome.bam

echo "Formatting primary ambiguous multi genome output BAM file: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view --verbosity 0 -@ 16 -b $primary_ambiguous_multi_genome > $prefix-primary_ambiguous_multi_genome.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary_ambiguous_multi_genome.bam > $prefix-primary_ambiguous_multi_genome_sorted.bam
rm -f $prefix-primary_ambiguous_multi_genome.bam

echo "Formatting Unmapped output BAM file: $(date '+%Y-%m-%d %H:%M:%S')"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view --verbosity 0 -@ 16 -b $unmapped > $prefix-unmapped.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-unmapped.bam > $prefix-unmapped-sorted.bam
rm -f $prefix-unmapped.bam


# cleaning up temporary files
rm -f $prefix-chunk_*
rm -f $prefix-part-* 
rm -f sam_chunk*
rm -f *.sam
rm -f $prefix-ambig.*