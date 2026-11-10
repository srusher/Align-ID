#!/usr/bin/env python3
import sys
from pathlib import Path

def name_to_taxid_dict(tax_map_file):
	"""Read the 3-column mapping file (accession, tax_id, name).

	Returns a dict: {name: taxid}
	    {'Escherichia coli': '562', 'Foo virus': '10'}

	Many accessions share one name/taxid, so duplicates collapse to one entry.
	"""

	name_to_taxid = {}

	with open(tax_map_file) as fh:
		for line in fh:
			columns = line.rstrip('\r\n').split('\t', 2)

			if len(columns) < 3 or not columns[1].strip().isdigit():
				continue  # blank line, header, or missing name

			taxid = columns[1].strip()
			name = columns[2].strip()

			name_to_taxid[name] = taxid

	return name_to_taxid

def accession_to_tax_id_dict(tax_map_file):
	"""Read the 3-column mapping file (accession, tax_id, name).

	Returns a dict: {accession: taxid}
	The accession is also stored without its version (NC_000913.3 -> NC_000913).
	"""

	contig_to_taxid = {}

	with open(tax_map_file) as fh:
		for line in fh:
			columns = line.rstrip('\r\n').split('\t', 2)

			if len(columns) < 2 or not columns[1].strip().isdigit():
				continue  # blank line or header

			accession = columns[0].strip()
			taxid = columns[1].strip()

			contig_to_taxid[accession] = taxid
			contig_to_taxid.setdefault(accession.split('.')[0], taxid)

	return contig_to_taxid


def parse_coverage(coverage_file, tax_map_file):
	"""Combine samtools coverage rows into whole-genome coverage per tax ID.

	Returns a dict of dicts:

	    {
	        '562': {
	            'n_contigs': 2,
	            'genome_length': 4000,
	            'covbases': 3500,
	            'coverage_pct': 87.5,
	            'mean_depth': 3.5,
	            'mean_mapq': 60.0,
	        },
	        '10': {...},
	    }

	mean_mapq is weighted by numreads per contig (0.0 if no reads map).
	Contigs with no taxid in the mapping file are skipped.
	Samtools columns used: startpos, endpos, numreads, covbases, meandepth, meanmapq.
	"""

	contig_to_taxid = accession_to_tax_id_dict(tax_map_file)

	# running totals per taxid
	totals = {}

	with open(coverage_file) as fh:
		for line in fh:
			if not line.strip() or line.startswith('#'):
				continue

			columns = line.rstrip('\r\n').split('\t')

			contig = columns[0]
			start = int(columns[1])
			end = int(columns[2])
			numreads = int(columns[3])
			covbases = int(columns[4])
			mean_depth = float(columns[6])
			mean_mapq = float(columns[8])

			taxid = contig_to_taxid.get(contig) or contig_to_taxid.get(contig.split('.')[0])

			if taxid is None:
				continue

			contig_length = end - start + 1

			if taxid not in totals:
				totals[taxid] = {'n_contigs': 0, 'genome_length': 0, 'covbases': 0,
				                 'depth_x_length': 0.0, 'numreads': 0, 'mapq_x_reads': 0.0}

			totals[taxid]['n_contigs'] += 1
			totals[taxid]['genome_length'] += contig_length
			totals[taxid]['covbases'] += covbases
			totals[taxid]['depth_x_length'] += mean_depth * contig_length
			totals[taxid]['numreads'] += numreads
			totals[taxid]['mapq_x_reads'] += mean_mapq * numreads

	# turn running totals into final numbers
	results = {}

	for taxid, t in totals.items():
		length = t['genome_length']
		reads = t['numreads']

		results[taxid] = {
			'n_contigs': t['n_contigs'],
			'genome_length': length,
			'covbases': t['covbases'],
			'coverage_pct': 100 * t['covbases'] / length if length else 0.0,
			'mean_depth': t['depth_x_length'] / length if length else 0.0,
			'mean_mapq': t['mapq_x_reads'] / reads if reads else 0.0,
		}

	return results

def main():

	prefix = sys.argv[1]
	tax_sum = sys.argv[2]
	coverage_file = sys.argv[3]
	tax_map = sys.argv[4]
	read_count_file = sys.argv[5] if len(sys.argv) > 5 else None

	out_file = f"{prefix}_taxonomy_report.tsv"

	name_dict = name_to_taxid_dict(tax_map)

	coverage_stats = parse_coverage(coverage_file, tax_map)

	total_reads = 0

	if read_count_file:

		with open(read_count_file) as fh:

			top_line = fh.readline().strip()

		total_reads = int(top_line)

	total_organism_reads = 0
	count = 0

	with open(tax_sum, 'r') as f:
		
		for line in f:

			count = count + 1

			if count == 1:

				continue
		
		num_reads = line.split('\t')[1].strip()

		total_organism_reads = total_organism_reads + int(num_reads)



	tax_dict = {}
	count = 0

	with open(tax_sum, 'r') as f:
		
		for line in f:

			count = count + 1

			#skip header line
			if count == 1:

				continue

			if not line.strip(): 
				
				continue

			species = line.split('\t')[0].strip()
			num_reads = int(line.split('\t')[1].strip())
			percent = float(line.split('\t')[2].strip())
			ambiguity = line.split('\t')[3].strip()

			if species not in tax_dict:

				tax_dict[species] = [0,0,0,0,0,0,0,0,0,0,0,0,0]

				tax_id = name_dict.get(species)
				stats = coverage_stats.get(tax_id, {})

				tax_dict[species][8] = stats.get('coverage_pct', 'NA')
				tax_dict[species][9] = stats.get('mean_depth', 'NA')
				tax_dict[species][10] = stats.get('genome_length', 'NA')
				tax_dict[species][11] = stats.get('n_contigs', 'NA')
				tax_dict[species][12] = stats.get('mean_mapq', 'NA')
			
			if ambiguity == 'none':

				values = tax_dict[species]

				values[0] = values[0] + num_reads
				values[1] = round(float(values[1]) + 100 * (num_reads/total_reads) if (read_count_file) else float(values[1]) + percent, 3)
				values[2] = num_reads
				values[3] = percent

				tax_dict[species] = values
			
			elif ambiguity == 'single_genome':
		
				values = tax_dict[species]

				values[0] = values[0] + num_reads
				values[1] = round(float(values[1]) + 100 * (num_reads/total_reads) if (read_count_file) else float(values[1]) + percent, 3)
				values[4] = num_reads
				values[5] = percent

				tax_dict[species] = values

			elif ambiguity == 'multi_genome':
		
				values = tax_dict[species]

				values[0] = values[0] + num_reads
				values[1] = round(float(values[1]) + 100 * (num_reads/total_reads) if (read_count_file) else float(values[1]) + percent, 3)
				values[6] = num_reads
				values[7] = percent

				tax_dict[species] = values

	with open(out_file, 'w') as out:

		out.write('species\tspecies_read_count\tspecies_read_perc\tambigutiy_none_read_count\tambiguity_none_percent\tambigutiy_single_read_count\tambiguity_single_percent\tambigutiy_multi_read_count\tambiguity_multi_percent\tgenome_cov(%)\tmean_depth\tgenome_len\tn_contigs\tmean_mapq\n')

		for key in tax_dict:

			out.write(f"{key}\t{tax_dict[key][0]}\t{tax_dict[key][1]}%\t{tax_dict[key][2]}\t{tax_dict[key][3]}%\t{tax_dict[key][4]}\t{tax_dict[key][5]}%\t{tax_dict[key][6]}\t{tax_dict[key][7]}%\t{tax_dict[key][8]}%\t{tax_dict[key][9]}\t{tax_dict[key][10]}\t{tax_dict[key][11]}\t{tax_dict[key][12]}\n")

if __name__ == "__main__":
	main()