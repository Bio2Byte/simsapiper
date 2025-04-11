#!/bin/bash
#Copy your_data_folder in simsapiper folder and edit this line to let magic happen!
#!/bin/bash
#SBATCH --job-name=simsa
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3G


data=toy_example

module load Nextflow/23.10.0

export APPTAINER_CACHEDIR=$VSC_SCRATCH_VO_USER/.apptainer


house=$(pwd)
now=`date +"%Y_%m_%d_%H_%M_%S"`
output_name=${data}_${now}_test
output_folder=$house/results/$output_name

mkdir -p $house/results/
mkdir -p $output_folder
echo 'Starting nextflow'

nextflow run simsapiper.nf \
    -profile server,withapptainer \
    --data $house/$data/data \
    --localModel 0.5 \
    --outFolder $output_folder \
    --apptainerPath  $VSC_SCRATCH_VO_USER/.apptainer \
    | tee  $output_folder/run_report_$output_name.nflog
sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $output_folder/run_report_$output_name.nflog)
nextflow log | grep $sessionName >> $output_folder/run_report_$output_name.nflog

