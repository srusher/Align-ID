#!/usr/bin/env python3
"""
Replace all strain-level tax IDs and below with the parent species tax ID
in seqid2taxid.map file - useful for custom taxonomic classification.
"""

from datetime import datetime
import sys

nodes_path = sys.argv[1]
seqid2taxid = sys.argv[2]
outfile   = sys.argv[3]

node_dict= {}
strain_species_dict = {}

print("Populating taxonomy node dictionary")
print(f"{datetime.now()}\n")

with open(nodes_path, 'r') as nodes:

    for line in nodes:

        parts = line.replace("\t", "").split("|")

        id = parts[0]
        parent_id = parts[1]
        clasification = parts[2]

        node_dict[id] = [parent_id,clasification]

print("Taxonomy node dictionary populated")
print(f"{datetime.now()}\n")

print("Extracting strain ID / species ID pairs")
print(f"{datetime.now()}\n")

for id in node_dict:

    clasification = node_dict[id][1]

    if clasification == 'strain' or clasification == 'no rank':

        strain_id = id
        continue_loop = True
        node = [id, node_dict[id][0], node_dict[id][1]]

        while continue_loop:

            node_id = node[0]
            parent_id = node[1]
            clasification = node[2]

            parent_node = [parent_id, node_dict[parent_id][0], node_dict[parent_id][1]]

            if node_id == parent_id:
                break

            if parent_node[2] == 'species':

                species_id = parent_id

                if node_id not in strain_species_dict:

                    strain_species_dict[strain_id] = species_id

                continue_loop = False
            
            elif parent_node is not None:

                node = parent_node

            else:

                continue_loop = False

print("Strain species pairs extracted")
print(f"{datetime.now()}\n")

print(f"Replacing all strain IDs with parent species IDs and writing to new mapping file: {outfile}")
print(f"{datetime.now()}\n")

count=0

with open(seqid2taxid, "r") as f_in, open(outfile, "w") as f_out:

    for line in f_in:
        count += 1
        parts = line.replace(" ", "").strip().split("\t")

        seqid = parts[0]
        taxid = parts[1]

        if taxid in strain_species_dict:

            new_taxid = strain_species_dict[taxid]
            f_out.write(f"{seqid}\t{new_taxid}\n")

        else:

            f_out.write(line)


print("All strain IDs replaced")
print(f"{datetime.now()}\n")