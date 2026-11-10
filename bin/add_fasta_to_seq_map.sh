#!/bin/bash
PYTHON_CONTAINER='https://depot.galaxyproject.org/singularity/python%3A3.10'
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)
cur_dir=$(pwd)
script_dir=$(dirname "$(realpath "$0")")

fasta_tax_file=$1

tax_dir="$script_dir/../assets/taxonomy"
map_file="$tax_dir/custom_seqid2taxid.map"
# tax_db="$tax_dir/tax_map.db"
# tax_db_no_strain="$tax_dir/tax_map_no_strain.db"


while IFS= read -r line || [[ -n "$line" ]]; do

    fasta=$(echo "$line" | cut -d '|' -f1)
    tax_id=$(echo "$line" | cut -d '|' -f2)

    if [ -f "$fasta" ]; then

        continue

    else

        echo "FASTA: '$fasta' does not exist"
        exit

    fi

    for id in $(grep '>' $fasta | cut -d '>' -f2 | cut -d ' ' -f1); do

        # removing any trailing version numbers from the accession number
        id=$(sed -E 's/\.[0-9]+$//' <<< "$id")

        dup_check=$(cut -f1 $map_file | grep "$id")

        if [[ -z $dup_check ]]; then

            echo -e "$id\t$tax_id" >> $map_file
        
        else

            echo "Sequence ID [$id] already in mapping file [$map_file] - skipping contig"

        fi

    done

done < "$fasta_tax_file"

# sqlite_pragmas="
# PRAGMA cache_size = -500000;
# PRAGMA temp_store_directory = '$cur_dir';
# "

# if [ -f "$tax_db" ]; then

#     singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
#         "
#         DROP TABLE IF EXISTS tax_map_custom;
#         CREATE TABLE tax_map_custom (seq_id TEXT, tax_id TEXT);
#         "   

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db <<EOF
# $sqlite_pragmas
# .mode tabs
# .import $map_file tax_map_custom
# EOF

#     singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
#         "
#         $sqlite_pragmas
#         CREATE TABLE tax_map_custom_new AS
#         SELECT tax_map_custom.seq_id, tax_map_custom.tax_id, names_dmp.name
#         FROM tax_map_custom
#         LEFT JOIN names_dmp ON (tax_map_custom.tax_id = names_dmp.tax_id);

#         DROP TABLE tax_map_custom;

# 	    ALTER TABLE tax_map_custom_new RENAME TO tax_map;
#         "

#     singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
#         "
#         $sqlite_pragmas
#         INSERT INTO tax_map
#         SELECT *
#         FROM tax_map_custom
#         "

# else

#     echo "No tax_map.db file found in taxonomy directory, checking for tax_map_no_strain.db"

# fi

# if [ -f "$tax_db_no_strain" ]; then

#     no_strain_tax_map=$tax_dir/seqid2taxid_no-strain_custom.map

# 	singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python3 $script_dir/remove_strains_from_seqid2taxid.py $script_dir/../assets/taxonomy/nodes.dmp $map_file $no_strain_tax_map

# 	tax_db="$tax_dir/tax_map_no_strain.db"

#     singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db \
#         "
#         DROP TABLE IF EXISTS tax_map_custom_no_strain;
#         CREATE TABLE tax_map_custom_no_strain (seq_id TEXT, tax_id TEXT);
#         "   

# singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $tax_db <<EOF
# $sqlite_pragmas
# .mode tabs
# .import $no_strain_tax_map tax_map_custom_no_strain
# EOF

#     singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
#         "
#         $sqlite_pragmas
#         CREATE TABLE tax_map_custom_no_strain_new AS
#         SELECT tax_map_custom_no_strain.seq_id, tax_map_custom_no_strain.tax_id, names_dmp.name
#         FROM tax_map_custom_no_strain
#         LEFT JOIN names_dmp ON (tax_map_custom_no_strain.tax_id = names_dmp.tax_id);

#         DROP TABLE tax_map_custom_no_strain;

# 	    ALTER TABLE tax_map_custom_no_strain_new RENAME TO tax_map_custom_no_strain;
#         "

#     singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
#         "
#         $sqlite_pragmas
#         INSERT INTO tax_map
#         SELECT *
#         FROM tax_map_custom_no_strain
#         "

# else

#     echo "No tax_map_no_strain.db file found in taxonomy directory -- skipping"

# fi

