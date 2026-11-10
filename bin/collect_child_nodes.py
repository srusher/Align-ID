import sys

input_nodes = sys.argv[1]
tax_id_filter = sys.argv[2]
outfile = "tax_id_child_nodes.tsv"
outfile2 = "tax_parent_child_map.tsv"

open(outfile, "w").close()

nodes_dict = {}

with open(input_nodes, "r") as fin:

	for line in fin:
		
		value, key = line.rstrip("\n").split("\t")

		if key in nodes_dict:

			nodes_dict[key] += "," + value
			   
		else:
			   
			nodes_dict[key] = value
		

all_nodes = []

with open(tax_id_filter, "r") as fin:
	
	for line in fin:
		
		all_child_nodes = []
		id = line.strip()
		loop_again = True

		if id in nodes_dict:

			all_child_nodes.append(id)
			tmp_child_ids = nodes_dict[id].split(',')


			while loop_again:

				count = 0

				for x in tmp_child_ids:

					count = count + 1

					# clear out tmp_child_ids list on first iteration
					if count == 1:

						tmp_child_ids = []

					if x in nodes_dict:

						all_child_nodes.append(x)

						tmp_child_ids.extend(nodes_dict[x].split(','))

					else:

						all_child_nodes.append(x)
				
				if not tmp_child_ids:

					loop_again = False

		else:

			all_child_nodes.append(id)
		
		# with open(outfile2, "a") as fin:

		# 	fin.write(id + '\t' + ",".join(list(dict.fromkeys(all_child_nodes))) + '\n')

		all_nodes.extend(list(dict.fromkeys(all_child_nodes)))
		
all_nodes = list(dict.fromkeys(all_nodes))

with open(outfile, "a") as fout:

	for i in all_nodes:

		fout.write(i + '\n')

	




		
		
