#!/bin/bash
module load Nextflow/23.10.0
house=$VSC_SCRATCH_VO_USER/simsapiper
data=toy_example
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name
mkdir -p $house/results
mkdir -p $output_folder
nextflow run simsapiper.nf \
    -profile hydra \
    --data $house/$data/data \
    --magic \
    --outFolder $output_folder \
    -with-trace $output_folder/timing_$output_name.txt \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog

#screen -S nextflowalign bash -c ./magic_hydra.sh