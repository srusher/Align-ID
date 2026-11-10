#!/usr/bin/env python3
"""
Extract accession IDs (without descriptions) from a FASTA file and
write them to an output file.

Usage:
	python extract_accessions.py input.fasta
	python extract_accessions.py input.fasta output.txt
"""

import sys
import os
import re
from pathlib import Path

def extract_scientific_names(names):

	name_dict = {}
	tax_id_dict = {}

	tax_id_dict['3136119'] = 'Bovine respiratory syncytial virus'
	tax_id_dict['11247'] = 'Bovine respiratory syncytial virus'
	tax_id_dict['488538'] = 'Candidatus Puniceispirillum marinum IMCC1322'

	with open(names, 'r') as f:

		for line in f:

			name = line.split('\t')[1].strip()
			tax_id = line.split('\t')[0].strip()

			if name not in name_dict:

				name_dict[name] = tax_id

			if tax_id not in tax_id_dict:

				tax_id_dict[tax_id] = name

	return name_dict, tax_id_dict

def extract_other_names(names, name_dict, tax_id_dict):

	names_list_lower = set([k.lower() for k in name_dict])

	with open(names, 'r') as f:

		for line in f:

			name = line.split('\t')[1].strip()
			tax_id = line.split('\t')[0].strip()

			if name not in name_dict and name not in names_list_lower:

				name_dict[name] = tax_id

			if tax_id not in tax_id_dict:

				tax_id_dict[tax_id] = name
	
	return name_dict, tax_id_dict

def build_seqid2taxid_dict(nucl_gb_seqid2taxid):

	seqid2taxid_dict = {}

	with open(nucl_gb_seqid2taxid, 'r') as f:

		for line in f:

			seq_id = line.split('\t')[0].strip()
			tax_id = line.split('\t')[1].strip()

			if seq_id in seqid2taxid_dict:

				continue
			
			else:

				seqid2taxid_dict[seq_id] = tax_id
	
	return seqid2taxid_dict

def extract_accessions(fasta_path, name_dict, tax_id_dict, seqid2taxid_dict):

	names_list_lower = set([k.lower() for k in name_dict])

	accessions = []
	unknown_accessions = []

	with open(fasta_path, 'r') as f:

		for line in f:

			if line.startswith('>'):

				header = line[1:].strip()
				# Accession is everything before the first whitespace
				tokens = header.split()

				accession = tokens[0] if tokens else header
				accession_trimmed = re.sub(r"\.\d+$", "", accession)

				tax_id = seqid2taxid_dict.get(accession_trimmed)
				name = tax_id_dict.get(tax_id) if tax_id else None

				if name is not None and tax_id is not None:
					accessions.append(f"{accession}\t{tax_id}\t{name}")

				else:

					description = tokens[1:]

					if description and description[0].lower() == 'mag:':

						description.pop(0)

					found = False
					brackets = False

					if " ".join(description).startswith('['):

						brackets = True

					for i in range(len(description), 0, -1):
						
						name = " ".join(description[:i]).strip() 

						if name.endswith(","):
							name = name[:-1]

						if name in name_dict or name.lower() in names_list_lower:

							accessions.append(f"{accession}\t{name_dict[name]}\t{name}")
							found = True
							break

					if brackets and not found:

						for i in range(len(description), 0, -1):
							
							name = " ".join(description[:i]).strip()  
							name = name[1:]
							name = name.replace(']', '', 1)

							if name.endswith(","):
								name = name[:-1]

							if name in name_dict or name.lower() in names_list_lower:

								accessions.append(f"{accession}\t{name_dict[name]}\t{name.capitalize()}")
								found = True
								break					
					
					if found:

						continue

					else:

						unknown_accessions.append(accession)

	# Default output filename if none provided

	output_path = Path(fasta_path).stem + '_tax_map.tsv'

	with open(output_path, 'w') as out:
		for acc in accessions:
			out.write(acc + '\n')

	print(f"Found {len(accessions)} accession IDs")
	print(f"Written to: {output_path}")

	if len(unknown_accessions) > 0:

		unknown_output_path = Path(fasta_path).stem + '_unknown.tsv'

		with open(unknown_output_path, 'w') as out:
			for acc in unknown_accessions:
				out.write(acc + '\n')		

		print(f"{len(unknown_accessions)} accession IDs could not be mapped to a name.")
		print(f"Written to: {unknown_output_path}")

	else:

		print(f"All accession IDs successfully mapped to a name.")

	return accessions


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Usage: python extract_accessions.py <input.fasta> [output.txt]")
		sys.exit(1)

	fasta_file = sys.argv[1]
	scientific_names_file = sys.argv[2]
	other_names_file = sys.argv[3]
	nucl_gb_seqid2taxid = sys.argv[4] if len(sys.argv) > 4 else None

	science_name_dict, tax_id_dict = extract_scientific_names(scientific_names_file)
	other_name_dict, other_tax_id_dict = extract_other_names(other_names_file, science_name_dict, tax_id_dict)
	seqid2taxid_dict = build_seqid2taxid_dict(nucl_gb_seqid2taxid)
	extract_accessions(fasta_file, other_name_dict, other_tax_id_dict, seqid2taxid_dict)