#!/bin/bash
module load Nextflow/22.10.4
house=$VSC_SCRATCH/
data=toy_expample
now=`date +"%Y_%m_%d_%H_%M_%S"`
output=${data}_${now}_test
mkdir -p $house/results/
mkdir -p $house/results/$output
nextflow run simsapiper.nf \
    -profile hydra,withsingularity \
    --data $house/$data \
    --magic \
    --outFolder $house/results/$output \
    |& tee  $house/results/$output/run_report_$output.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $house/results/$output/run_report_$output.nflog)
nextflow log | grep $sessionName >> $house/results/$output/run_report_$output.nflog