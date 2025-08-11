#!/bin/bash
#SBATCH --job-name=simsa
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3G

# Copy this file in the same folder as your_data_folder to let magic happen!

# It should look like this:
## toy_example
## | data/
## |   | seqs/
## |   |   | mysequences.fasta
## | magic_apptainer_slurm.sh

# Edit this line to the complete path of yor copy SIMSApiper directory. Do not add any spaces!
simsa_dir=/Users/someone/workspace/simsapiper

# Edit this line to give your alignment a name! 
name=toy_example

#If nextflow/apptainer is not on per default, load it here
module load Nextflow/23.10.0

#set apptainer cache to one place with lots of storage avoid keeping multiple copies
apptainercache=$VSC_SCRATCH_VO_USER/.apptainer
export APPTAINER_CACHEDIR=$apptainercache


#####################magic_zone##############################


data=$(pwd)
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${name}_${now}_test
output_folder=$data/results/$output_name

mkdir -p $data/results/
mkdir -p $output_folder
echo 'Starting nextflow'

nextflow run simsapiper.nf \
    -profile server \
    --data $house/$data/data \
    --minimagic \
    --outFolder $output_folder \
    --apptainerPath  $apptainercache \
    | tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog
