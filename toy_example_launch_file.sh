#!/bin/bash
house=$(pwd)
data=toy_example
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=toy_example_${now}
output_folder=$house/results/$output_name
mkdir $house/results
mkdir $house/results/$output_name
nextflow run simsapiper.nf \
    -profile hpc,withsingularity \
    --data $house/$data/data \
    --seqQC 5 \
    --dropSimilar 90 \
    --createSubsets 30 \
    --minSubsetID 20 \
    --retrieve true \
    --model true \
    --strucQC 5 \
    --dssp true \
    --squeeze "H,G,E" \
    --reorder true \
    --outFolder $output_folder \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog
