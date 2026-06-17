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

echo "Parsing BAM headers"

# printing bam headers to all primary alignment output sam files
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_all 
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_unambiguous
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_single_genome
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_multi_genome
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $unmapped


#############################################################################################
###                                 SAM Parsing Section                                  ####
#############################################################################################

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE total_reads (num_reads FLOAT);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO total_reads
	SELECT COUNT(read_id)
	FROM sam_complete
	WHERE flag = 4 or flag = 16 or flag = 0;
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unmapped (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO unmapped
	SELECT read_id, flag, ref_id, mapq
	FROM sam_complete
	WHERE ref_id = 'unclassified'
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE primary_multimap (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER, as_tag TEXT, tax_id INTEGER, FOREIGN KEY (ref_id) REFERENCES tax_map (seq_id));
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO primary_multimap
	SELECT sam_complete.read_id, sam_complete.flag, sam_complete.ref_id, sam_complete.mapq, sam_complete.as_tag, tax_map.tax_id
	FROM sam_complete
	JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:P' and mapq = 0
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE secondary_mappings (sec_read_id TEXT, sec_flag INTEGER, sec_ref_id TEXT, sec_mapq INTEGER, sec_as_tag TEXT, sec_tax_id INTEGER, FOREIGN KEY (sec_read_id) REFERENCES primary_multimap (read_id), FOREIGN KEY (sec_ref_id) REFERENCES tax_map (seq_id));
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO secondary_mappings
	SELECT sam_complete.read_id, sam_complete.flag, sam_complete.ref_id, sam_complete.mapq, sam_complete.as_tag, tax_map.tax_id
	FROM sam_complete
	JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:S' and read_id IN (SELECT read_id FROM primary_multimap)
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE primary_secondary_aggregate (read_id TEXT, sec_tax_id INTEGER, sec_tax_id_count INTEGER, prim_tax_id INTEGER);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO primary_secondary_aggregate
	SELECT secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id, COUNT(secondary_mappings.sec_tax_id), primary_multimap.tax_id
	FROM secondary_mappings
	JOIN primary_multimap ON (secondary_mappings.sec_read_id = primary_multimap.read_id) 
	WHERE secondary_mappings.sec_as_tag = primary_multimap.as_tag or cast(substr(secondary_mappings.sec_as_tag, 6) AS INTEGER) > cast(substr(primary_multimap.as_tag, 6) AS INTEGER)
	GROUP BY secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id
	" 

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE primary_unambiguous (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER, as_tag TEXT, tax_id INTEGER, FOREIGN KEY (ref_id) REFERENCES tax_map (seq_id));
	"

# had to change the condition for WHERE clause - a MAPQ of 0 doesn't guarantee a secondary alignment with a >= alignment score (AS) as the primary alignment meaning that the read may not be accounted for in the ambiguous (single or multi genome) queries; so in order to grab all unambiguous primary alignments we will use the tp_tag and make sure the read_id does not appear in the primary_secondary_aggregate table
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO primary_unambiguous
	SELECT sam_complete.read_id, sam_complete.flag, sam_complete.ref_id, sam_complete.mapq, sam_complete.as_tag, tax_map.tax_id
	FROM sam_complete
	JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE sam_complete.tp_tag = 'tp:A:P' AND sam_complete.read_id NOT IN (
		SELECT read_id
		FROM primary_secondary_aggregate
	)
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE single_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_tax_id_count INTEGER, prim_tax_id INTEGER);
	"	

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO single_genome_hits_by_read
	SELECT read_id, sec_tax_id, sec_tax_id_count, prim_tax_id
	FROM primary_secondary_aggregate
	WHERE sec_tax_id = prim_tax_id 
	GROUP BY read_id, prim_tax_id
	HAVING COUNT(read_id) = 1
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE single_genome_hits_by_taxid (prim_tax_id INTEGER, total_hits);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO single_genome_hits_by_taxid
	SELECT prim_tax_id, COUNT(prim_tax_id)
	FROM single_genome_hits_by_read
	GROUP BY prim_tax_id
	" 

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE multi_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_tax_id_count INTEGER, prim_tax_id INTEGER);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO multi_genome_hits_by_read
	SELECT read_id, sec_tax_id, sec_tax_id_count, prim_tax_id
	FROM primary_secondary_aggregate
	WHERE read_id NOT IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	) 
	" 

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE multi_genome_hits_by_taxid (prim_tax_id INTEGER, total_hits);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO multi_genome_hits_by_taxid
	SELECT prim_tax_id, COUNT(DISTINCT read_id)
	FROM multi_genome_hits_by_read
	GROUP BY prim_tax_id
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unambiguous_hits_by_taxid (prim_tax_id INTEGER, total_hits);
	"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    INSERT INTO unambiguous_hits_by_taxid
	SELECT tax_id, COUNT(DISTINCT read_id)
	FROM primary_unambiguous
	GROUP BY tax_id;
	"

echo -e "species\tnum_reads\tpercent_reads\tambiguity" > "$summary_file"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT unmapped.ref_id, COUNT(unmapped.read_id), ROUND(COUNT(unmapped.read_id) / total_reads.num_reads * 100,6) as percentage
    FROM unmapped
	LEFT JOIN total_reads
	GROUP BY ref_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT unambiguous_hits_by_taxid.prim_tax_id, unambiguous_hits_by_taxid.total_hits, ROUND(unambiguous_hits_by_taxid.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM unambiguous_hits_by_taxid
	LEFT JOIN total_reads
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT single_genome_hits_by_taxid.prim_tax_id, single_genome_hits_by_taxid.total_hits, ROUND(single_genome_hits_by_taxid.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM single_genome_hits_by_taxid
	LEFT JOIN total_reads
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "single_genome"}' >> "$summary_file"

singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT multi_genome_hits_by_taxid.prim_tax_id, multi_genome_hits_by_taxid.total_hits, ROUND(multi_genome_hits_by_taxid.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM multi_genome_hits_by_taxid
	LEFT JOIN total_reads
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "multi_genome"}' >> "$summary_file"

count=0

while IFS= read -r line || [[ -n "$line" ]]; do

	#skipping header line
	if [ $count -lt 1 ]; then

		count=$((count+1))
		continue
	
	fi

	tax_id=$(echo "$line" | cut -f1)

	if [[ $tax_id == "unclassified" ]]; then

		continue
	
	elif [[ "$tax_id" =~ [a-zA-Z] ]]; then # if the tax_id var contains letters then it is most likely a nonstandard reference

		continue

	fi

	# this line is grabbing the scientific name from the names.dmp taxonomy file
	species_name=$(grep "^$tax_id.[|].*scientific name" $taxa_names | cut -d "|" -f 2 | tr -d '\t')

	sed -i "s|^\<$tax_id\>|$species_name|g" "$summary_file"

done < "$summary_file"

# grab all primary alignments
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P'
	" | sed 's/!##########!/\t/g' >> $primary_all

# grab all nonambiguous primary alignments
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM primary_unambiguous
	)
	" | sed 's/!##########!/\t/g' >> $primary_unambiguous

# grab all primary ambiguous single genome hits
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	)
	" | sed 's/!##########!/\t/g' >> $primary_ambiguous_single_genome

# grab all primary ambiguous multi-genome hits
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id NOT IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	) 
	" | sed 's/!##########!/\t/g' >> $primary_ambiguous_multi_genome

# grab all unmapped reads - unmapped reads will often contain empty aux fields so we need to remove trailing tab characters from the SQL query
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE read_id IN (
		SELECT read_id
		FROM unmapped
	) 
	" | sed 's/!##########!/\t/g' | sed 's/\t*$//' >> $unmapped

#############################################################################################
#############################################################################################
#############################################################################################


echo "Formatting all primary output BAM file"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_all > $prefix-primary-all.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-all.bam > $prefix-primary-all-sorted.bam
rm -f $prefix-primary-all.bam

echo "Formatting primary unambiguous single genome output BAM file"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_unambiguous > $prefix-primary-unambiguous.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-unambiguous.bam > $prefix-primary-unambiguous-sorted.bam
rm -f $prefix-primary-unambiguous.bam

echo "Formatting primary ambiguous single genome output BAM file"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_ambiguous_single_genome > $prefix-primary_ambiguous_single_genome.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary_ambiguous_single_genome.bam > $prefix-primary_ambiguous_single_genome-sorted.bam
rm -f $prefix-primary_ambiguous_single_genome.bam

echo "Formatting primary ambiguous multi genome output BAM file"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_ambiguous_multi_genome > $prefix-primary_ambiguous_multi_genome.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary_ambiguous_multi_genome.bam > $prefix-primary_ambiguous_multi_genome-sorted.bam
rm -f $prefix-primary_ambiguous_multi_genome.bam

echo "Formatting Unmapped output BAM file"
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $unmapped > $prefix-unmapped.bam
singularity -q exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-unmapped.bam > $prefix-unmapped-sorted.bam
rm -f $prefix-unmapped.bam

# cleaning up temporary files
rm -f $prefix-chunk_*
rm -f $prefix-part-* 
rm -f sam_chunk*
rm -f *.sam
rm -f $prefix-ambig.*