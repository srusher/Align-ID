#!/bin/bash
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

local_nodes_db=$1
nodes=$2
my_tax_ids=$3

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

# if the DB does not exist, have sqlite create the DB and the TAX_IDS table
singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db \
"
CREATE TABLE IF NOT EXISTS TAX_IDS (
    parent_id  TEXT,
    child_ids  TEXT
);
"

parent_id_column=$(singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "SELECT parent_id FROM TAX_IDS")

while IFS= read -r id || [[ -n "$id" ]]; do #iterating through tax ids listed in the tax id input file

    data=$(echo $parent_id_column | grep -w "$id" )

    # if this parent ID entry already exists in the DB then go to the next one
    if [[ -n $data ]]; then

        echo "Tax ID: $id already in database"
        continue

    fi

    singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "INSERT INTO TAX_IDS (parent_id, child_ids) VALUES ("$id", '');"

    parent_id="$id"

    tax_array=("$parent_id") #creating array for tax ids
    loop_again=true
    count_other=0
    all_child_ids=""

    while $loop_again; do

        count=0

        for i in "${tax_array[@]}"; do

            if [[ $count -eq 0 ]]; then

                tax_array=() #clearing out array to rebuild it with child tax ids for the next for loop iteration
                ((count++))

            fi

            child_ids="$(awk -v id="$i" '$3 == id { print $1 }' $nodes)" #finding child nodes that have the "parent tax id" column set to current tax id

            if [[ ! -z $child_ids ]]; then #checking to see if we found any child nodes

                for x in $child_ids; do

                    all_child_ids=$all_child_ids"$x,"
                    tax_array+=("$x")
                    # if [[ $count_other -eq 0 ]]; then
                        
                    #     all_child_ids=$all_child_ids"$x,"
                    #     # this should only execute once per specified parent tax id
                    #     singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "UPDATE TAX_IDS SET child_ids = "$x" WHERE parent_id = "$parent_id";"
                    #     ((count_other++))

                    #     tax_array+=("$x") #appending child tax ids to array to use in
                    
                    # else

                    #     singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "UPDATE TAX_IDS SET child_ids = child_ids || ',' || "$x" WHERE parent_id = "$parent_id";"

                    #     tax_array+=("$x") #appending child tax ids to array to use in

                    # fi

                done                        
            
            # if there are no child IDs from the original parent ID then add the parent ID as the only child ID value
            elif [[ $count -eq 1 ]]; then

                singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "UPDATE TAX_IDS SET child_ids = '$i' WHERE parent_id = '$parent_id';"
                continue

            else

                continue

            fi

        done

        if [ ${#tax_array[@]} -eq 0 ]; then

            loop_again=false

        fi

    done

    if [[ -n $all_child_ids ]]; then

        all_child_ids="${all_child_ids%,}"
        singularity -q exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "UPDATE TAX_IDS SET child_ids = '$all_child_ids' WHERE parent_id = '$parent_id';"

    fi        

done < "$my_tax_ids"

# create placeholder output text file that tells nextflow to resume the workflow
touch complete.txt
