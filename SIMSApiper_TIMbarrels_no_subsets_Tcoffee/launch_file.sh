#!/bin/bash
house=$(pwd)
data=SIMSApiper_TIMbarrels_no_subsets_Tcoffee
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name
mkdir -p $house/results
mkdir -p $output_folder
nextflow run simsapiper.nf \
    -profile hpc,withsingularity\
    --data $house/$data/data \
    --seqQC 5 \
    --dropSimilar 90 \
    --createSubsets false \
    --retrieve true \
    --model true \
    --strucQC 5 \
    --dssp true \
    --squeeze "H,G,E" \
    --squeezePerc 70 \
    --reorder true \
    --outFolder $output_folder \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog
