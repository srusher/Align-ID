#!/bin/bash

PYTHON_CONTAINER='https://depot.galaxyproject.org/singularity/python%3A3.10'
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

cur_dir=$(pwd)
script_dir=$(dirname "$(realpath "$0")")
tax_dir="$script_dir/../assets/taxonomy"

strain_arg=$1

cd $tax_dir

# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

# gzip -d nucl_wgs.accession2taxid.gz

# awk '{print $1 "\t" $3}' nucl_wgs.accession2taxid > seqid2taxid.map

# rm nucl_wgs.accession2taxid

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

gzip -d nucl_gb.accession2taxid.gz

awk '{print $1 "\t" $3}' nucl_gb.accession2taxid > nucl_gb_seqid2taxid.map

rm nucl_gb.accession2taxid

# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

# gzip -d nucl_wgs.accession2taxid.gz

# awk '{print $1 "\t" $3}' nucl_wgs.accession2taxid > nucl_wgs_seqid2taxid.map

# rm nucl_wgs.accession2taxid

# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

tar -xvf taxdump.tar.gz

rm taxdump.tar.gz

rm citations.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp images.dmp readme.txt gc.prt

seqid2taxid_map=$tax_dir/seqid2taxid.map
tax_db="$tax_dir/tax_map.db"

# removing header row from file
sed -i '1d' $seqid2taxid_map

if [[ $strain_arg == 'no-strain' ]]; then

	singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python3 $script_dir/remove_strains_from_seqid2taxid.py $tax_dir/nodes.dmp $tax_dir/seqid2taxid.map $tax_dir/seqid2taxid_no-strain.map

	seqid2taxid_map=$tax_dir/seqid2taxid_no-strain.map
	# tax_db="$tax_dir/tax_map_no_strain.db"

fi

########################

# tax_map DB generation phase

cur_dir="$(pwd)"

# sqlite_pragmas="
# PRAGMA cache_size = -50000000;
# PRAGMA temp_store_directory = '$cur_dir';
# "

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
# 	"
# 	DROP TABLE IF EXISTS tax_map;
# 	CREATE TABLE tax_map (seq_id TEXT, tax_id TEXT);
# 	"

# echo "Taxonomy import START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db <<EOF
# $sqlite_pragmas
# .mode tabs
# .import $seqid2taxid_map tax_map
# EOF

# echo "Taxonomy import END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "seq ID dedup START: $(date '+%Y-%m-%d %H:%M:%S')"
# # removing duplicate rows from seq id to tax id map file
# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
# 	"
# 	$sqlite_pragmas
# 	DELETE FROM tax_map
# 	WHERE rowid NOT IN (
# 		SELECT MIN(rowid)
# 		FROM tax_map
# 		GROUP BY seq_id	
# 	)
# 	"

# echo "seq ID dedup END: $(date '+%Y-%m-%d %H:%M:%S')"


echo "Parsing names.dmp taxonomy file START: $(date '+%Y-%m-%d %H:%M:%S')"

sed -i 's/\t|\t/\t/g' "$tax_dir/names.dmp"

sed -i $'s/\t|$//' "$tax_dir/names.dmp"

grep -P "\tscientific name" "$tax_dir/names.dmp" > "$tax_dir/names_scientific.dmp"

grep -vP "\tauthority" "$tax_dir/names.dmp" | grep -vP "\tscientific name" > "$tax_dir/names_other.dmp"

echo "Parsing names.dmp taxonomy file END: $(date '+%Y-%m-%d %H:%M:%S')"

# awk -F'\t' -v OFS='\t' '{print $1, $2}' "$tax_dir/names_new_1.dmp" > "$tax_dir/names_new_2.dmp"

# echo "Formatting names_dmp for DB Import END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "Importting names_dmp into DB START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
# 	"
# 	DROP TABLE IF EXISTS names_dmp;
# 	CREATE TABLE names_dmp (tax_id TEXT, name TEXT);
# 	"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db <<EOF
# $sqlite_pragmas
# .mode tabs
# .import "$tax_dir/names_new_2.dmp" names_dmp
# EOF

# echo "Importting names_dmp into DB END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "Adding tax_name column to tax_map table START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
# 	"
# 	$sqlite_pragmas
# 	CREATE TABLE tax_map_new AS
# 	SELECT tm.seq_id, tm.tax_id, names_dmp.name
# 	FROM tax_map AS tm
# 	LEFT JOIN names_dmp ON (tm.tax_id = names_dmp.tax_id);

# 	DROP TABLE tax_map;

# 	ALTER TABLE tax_map_new RENAME TO tax_map;

# 	"

# echo "Adding tax_name column to tax_map table END: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "Creating index on tax_map(seq_id) START: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
# 	"
# 	$sqlite_pragmas
# 	CREATE INDEX IF NOT EXISTS idx_tax_map_seq_id ON tax_map(seq_id);
# 	"

# echo "Creating index on tax_map(seq_id) END: $(date '+%Y-%m-%d %H:%M:%S')"

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
# 	"
# 	$sqlite_pragmas
# 	VACUUM;
# 	"

# rm $tax_dir/names_new_1.dmp
# rm $tax_dir/names_new_2.dmp