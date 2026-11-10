#!/bin/bash
PYTHON_CONTAINER="https://depot.galaxyproject.org/singularity/python%3A3.12.12"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
single_end=$2
python_script=$3
mode=$4
headcrop=$5
tailcrop=$6
which_reads=$7
template_seqs=$8
reads=$9

echo $reads
r1=$(echo $reads | cut -d ' ' -f1)
r2=$(echo $reads | cut -d ' ' -f2)

if [[ $mode == 'fixed-length' ]]; then

	if [[ $which_reads == 'r1' || $which_reads == 'both' ]]; then

		singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python $python_script --fastq $r1 --read-type forward --headcrop $headcrop --tailcrop $tailcrop --output $prefix'_manual_trim_r1.fastq'

	elif [[ $which_reads == 'r2' || $which_reads == 'both' ]]; then

		singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python $python_script --fastq $r2 --read-type reverse --headcrop $headcrop --tailcrop $tailcrop --output $prefix'_manual_trim_r2.fastq'

	fi

elif [[ $mode == 'seq-match' ]]; then

	if [[ $which_reads == 'r1' || $which_reads == 'both' ]]; then

		singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python $python_script --fastq $r1 --read-type forward --sequences $template_seqs --output $prefix'_manual_trim_r1.fastq'	

	elif [[ $which_reads == 'r2' || $which_reads == 'both' ]]; then

		singularity -q exec --bind $bind_dir $PYTHON_CONTAINER python $python_script --fastq $r2 --read-type reverse --sequences $template_seqs --output $prefix'_manual_trim_r2.fastq'

	fi

fi

if [[ $single_end == 'false' && $which_reads == 'r1' ]]; then

	gzip -f $prefix'_manual_trim_r1.fastq'
	cp $r2 ./$prefix'_manual_trim_r2.fastq.gz'

elif [[ $single_end == 'false' && $which_reads == 'r2' ]]; then

	gzip -f $prefix'_manual_trim_r2.fastq'
	cp $r1 ./$prefix'_manual_trim_r1.fastq.gz'

else

	gzip -f *_manual_trim*.fastq

fi