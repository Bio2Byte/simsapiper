#!/bin/bash


#Copy your_data_folder in simsapiper folder and edit this line to let magic happen!
data=toy_example

house=$(pwd)
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/$data/results/$output_name

mkdir -p $output_folder
echo 'Starting nextflow'

nextflow run simsapiper.nf \
    -profile standard,withdocker \
    --data $house/$data/data \
    --magic \
    --outFolder $output_folder \
    | tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog


#To run in background: screen -S nextflowalign bash -c ./magic_hydra.sh