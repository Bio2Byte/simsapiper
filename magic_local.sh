#!/bin/bash

# Copy this file in the same folder as your_data_folder to let magic happen!

# It should look like this:
## toy_example
## | data/
## |   | seqs/
## |   |   | mysequences.fasta
## | magic_local.sh

# Edit this line to the complete path of yor copy SIMSApiper directory. Do not add any spaces!
simsa_dir=/Users/someone/workspace/simsapiper

# Edit this line to give your alignment a name! 
name=toy_example


# Remember to start Docker!
# You may have to run 
##  chmod +x magic_local.sh
# Now run
## ./magic_local.sh

#####################magic_zone##############################


data=$(pwd)
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${name}_${now}
output_folder=$data/results/$output_name

mkdir -p $data/results/
mkdir -p $output_folder
echo 'Starting nextflow'

nextflow run $simsa_dir/simsapiper.nf \
    -profile standard,withdocker \
    --data $data/data \
    --magic \
    --outName $name \
    --outFolder $output_folder \
    | tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog


#To run in background: screen -S nextflowalign bash -c ./magic_hydra.sh
