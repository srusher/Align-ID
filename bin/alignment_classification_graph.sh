#!/bin/bash
R_CONTAINER="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/R/ggplot_dplyr/ggplot_dplyr.sif"

prefix=$1
alignment_report=$2
rscript=$3
data_type=$4

input_tsv="$prefix-rscript-input.csv"

singularity exec --bind /scicomp $R_CONTAINER Rscript $rscript $alignment_report "$(pwd)" $prefix $data_type

mv plot.png $prefix-taxonomic-distribution-plot_mqc.png