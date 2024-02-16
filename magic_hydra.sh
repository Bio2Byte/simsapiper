#!/bin/bash
module load Nextflow/23.04.2
house=/scratch/brussel/106/vsc10611/nf_tcoffee_bitbucket/
data=Daan_beta_barrels
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${now}_no_subset
output_folder=$house/results/$output_name
mkdir -p $house/$data/results
mkdir -p $output_folder
nextflow run simsapiper.nf \
    -profile hpc,withsingularity \
    --data $house/$data/data \
    --seqQC 5 \
    --createSubsets true \
    --minSubsetID "min" \
    --strucQC 5 \
    --reorder true \
    --outFolder $output_folder \
    -with-trace $output_folder/timing_$output_name.txt \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog
