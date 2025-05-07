#!/bin/bash
#SBATCH --job-name=simsa
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G

#Copy your_data_folder in simsapiper folder and edit this line to let magic happen!
data=toy_example

#If nextflow/apptainer is not on per default, load it here
module load Nextflow/23.10.0

#set apptainer cache to one place with lots of storage avoid keeping multiple copies
apptainercache=$VSC_SCRATCH_VO_USER/.apptainer
export APPTAINER_CACHEDIR=$apptainercache

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