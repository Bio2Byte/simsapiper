#!/bin/bash
module load Nextflow/23.04.2
house=$(pwd)
data=toy_example
now=`date +"%Y_%m_%d_%H_%M_%S"`
output=toy_example_${now}
mkdir $house/results
mkdir $house/results/$output
nextflow run simsapiper.nf \
    -profile hpc,withsingularity \
    --data $house/$data/data \
    --seqQC 5 \
    --dropSimilar 90 \
    --createSubsets 30 \
    --retrieve true \
    --model true \
    --strucQC 5 \
    --dssp true \
    --squeeze "H,G,E" \
    --reorder true \
    --outFolder $house/results/$output \
    |& tee  $house/results/$output/$output.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $house/results/$output/$output.nflog)
nextflow log | grep $sessionName >> $house/results/$output/$output.nflog
