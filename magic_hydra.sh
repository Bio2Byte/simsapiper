#!/bin/bash
module load Nextflow/23.04.2
house=$VSC_SCRATCH/simsapiper
data=toy_example
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name
mkdir -p $house/results
mkdir -p $output_folder
nextflow run simsapiper.nf \
    -profile hydra,withsingularity \
    --data $house/$data/data \
    --magic \
    --outFolder $output_folder \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog