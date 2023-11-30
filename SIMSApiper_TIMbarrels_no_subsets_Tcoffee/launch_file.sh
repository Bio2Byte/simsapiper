#!/bin/bash
module load Nextflow/22.10.4
house=$VSC_SCRATCH/nf_tcoffee/SIMSApiper_TIMbarrels
data=data
now=`date +"%Y_%m_%d_%H_%M_%S"`
output=${now}
mkdir $house/results
mkdir $house/results/$output
nextflow run simsapiper.nf \
    -profile hpc,withsingularity \
    --data $house/$data \
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
    --outFolder $house/results/$output \
    |& tee  $house/results/$output/$output.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $house/results/$output/$output.nflog)
nextflow log | grep $sessionName >> $house/results/$output/$output.nflog