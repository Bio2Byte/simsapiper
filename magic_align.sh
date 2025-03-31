#!/bin/bash
#Copy your_data_folder in simsapiper folder and edit this line to let magic happen!
data=toy_example

house=$(pwd)
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name

mkdir -p $house/results/
mkdir -p $output_folder
echo 'Starting nextflow'

nextflow run simsapiper.nf \
    -profile server,withsingularity \
    --data $house/$data/data \
    --localModel 0.5 \
    --outFolder $output_folder \
    --apptainerPath  "$VSC_SCRATCH_VO_USER/.apptainer" \
    | tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog

