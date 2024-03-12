#!/bin/bash
house=$(pwd)
data=toy_example
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name
mkdir -p $house/results/
mkdir -p $output_folder
echo 'Starting nextflow'
nextflow run simsapiper.nf \
    -profile server,withsingularity \
    --data $house/$data/data \
    --magic \
    --outFolder $output_folder \
    -with-trace $output_folder/timings_$output_name.txt \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog

#screen -S nextflowalign bash -c ./magic_align.sh
