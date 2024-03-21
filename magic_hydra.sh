#!/bin/bash
#screen -S nextflowalign bash -c ./magic_hydra.sh

data=toy_example

module load Nextflow/23.04.2
house=$VSC_SCRATCH_VO_USER/simsapiper
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
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog
