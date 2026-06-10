#!/bin/bash

fastq_dir=$(realpath $1)
delimiter=$2
read_type=$3
single_end=$4
samplesheet=$(realpath $5)


if [[ $read_type == 'short' ]]; then

	echo 'sample,fastq_1,fastq_2' > $samplesheet

else

	echo 'sample,fastq_long' > $samplesheet

fi

declare -A prefix_count

read_1=''
read_2=''

for i in $(ls $fastq_dir); do

	prefix=$(echo $i | cut -d "$delimiter" -f1)

	if [[ $single_end == 'true' ]]; then

		if [[ $read_type == 'short' ]]; then

			echo "$prefix,$fastq_dir/$i," >> $samplesheet
		
		else

			echo "$prefix,$fastq_dir/$i" >> $samplesheet
		
		fi

	elif [[ -v prefix_count["$prefix"] ]]; then

		read_2="$fastq_dir/$i"

		echo "$prefix,$read_1,$read_2" >> $samplesheet
	
	else

		prefix_count["$prefix"]=1
		read_1="$fastq_dir/$i"

	fi

done


