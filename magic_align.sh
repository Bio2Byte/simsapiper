#!/bin/bash
module load Nextflow/23.04.2
house=$VSC_SCRATCH_VO_USER/
data=/scratch/brussel/vo/000/bvo00023/vsc10611/simsapiper/tests
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${now}
output_folder=$data/results/$output_name
mkdir -p $data/results/
mkdir -p $output_folder
echo 'Starting nextflow'
nextflow run /scratch/brussel/vo/000/bvo00023/vsc10611/simsapiper/simsapiper.nf \
    -profile hydra \
    --data $data/data \
    --seqQC 5 \
    --dropSimilar 90 \
    --createSubsets 30 \
    --minSubsetID "min" \
    --retrieve true \
    --model true \
    --strucQC 5 \
    --reorder true \
    --outFolder $output_folder \
    -with-trace $output_folder/timings_$output_name.txt \
    |& tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog

#screen -S nextflowalign bash -c ./magic_align.sh
