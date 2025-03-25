#!/bin/bash
#SBATCH --job-name=simsa
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G


data=toy_example

module load Nextflow/23.10.0
house=$VSC_SCRATCH_VO_USER/simsapiper
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name
mkdir -p $house/results
mkdir -p $output_folder

nextflow run simsapiper.nf \
    -profile hydra \
    --data $house/$data/data \
    --magic \
    --outFolder $output_folder