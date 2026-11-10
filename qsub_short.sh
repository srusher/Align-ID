#!/bin/bash
source /etc/profile

module load nextflow
module load python/3.11.3

config=$1

# primary nextflow execution command: alter as needed
nextflow run main.nf -profile singularity,cluster_submit,short
