#!/bin/bash

PYTHON_CONTAINER='https://depot.galaxyproject.org/singularity/python%3A3.10'

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

cur_dir=$(pwd)
script_dir=$(dirname "$(realpath "$0")")
tax_dir="$script_dir/../assets/taxonomy"

strain_arg=$1

cd $tax_dir

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

gzip -d nucl_wgs.accession2taxid.gz

awk '{print $1 "\t" $3}' nucl_wgs.accession2taxid > seqid2taxid.map

rm nucl_wgs.accession2taxid

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

tar -xvf taxdump.tar.gz

rm taxdump.tar.gz

rm citations.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp images.dmp readme.txt gc.prt

if [[ $strain_arg == 'no-strain' ]]; then

	singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python3 $script_dir/remove_strains_from_seqid2taxid.py $script_dir/../assets/taxonomy/nodes.dmp $script_dir/../assets/seqid2taxid.map $script_dir/../assets/seqid2taxid_no-strain.map

fi

if [[ $strain_arg == 'no-strain' ]]; then

	singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python3 $script_dir/remove_strains_from_seqid2taxid.py nodes.dmp seqid2taxid.map seqid2taxid_no-strain.map

fi