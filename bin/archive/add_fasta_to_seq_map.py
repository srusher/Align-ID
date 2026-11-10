#!/usr/bin/env python3
import re
import sys
from pathlib import Path

if len(sys.argv) != 2:
    sys.exit(f"Usage: {sys.argv[0]} <fasta_tax_file>")

fasta_tax_file = Path(sys.argv[1])
script_dir = Path(__file__).resolve().parent

map_file = script_dir / ".." / "assets" / "taxonomy" / "seqid2taxid.map"

# Load existing accession -> taxid mapping
seqid2taxid = {}

if map_file.exists():
    
    for line in map_file.read_text().splitlines():
        
        if line:
            
            acc, tax_id = line.split("\t")
            seqid2taxid[acc] = tax_id[0] if rest else ""

# Process each "fasta|taxid" line
with open(fasta_tax_file) as infile, open(map_file, "a") as outfile:
    
    for line in infile:
        
        line = line.strip()
        
        if not line:
            continue

        fasta, tax_id = line.split("|")[:2]

        if not Path(fasta).is_file():
            sys.exit(f"FASTA: '{fasta}' does not exist")

        with open(fasta) as fh:
            
            for fasta_line in fh:
                
                if not fasta_line.startswith(">"):
                    continue

                # first token of the header, minus any trailing version (.1, .2, ...)
                acc = fasta_line[1:].split()[0]
                acc = re.sub(r"\.[0-9]+$", "", acc)

                if acc in seqid2taxid:
                    
                    print(f"Sequence ID [{acc}] already in mapping file [{map_file}] - skipping contig")
                    continue

                outfile.write(f"{acc}\t{tax_id}\n")
                seqid2taxid[acc] = tax_id