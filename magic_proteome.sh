#!/bin/bash
#SBATCH --job-name=simsaht
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G

#use this if you want to align many ortholog proteins containing on average less then 100 proteins per set


house=$VSC_SCRATCH_VO_USER/proteome


## create a file containing all orthogroup identifiers you want to align individually, each in a new line

## each ortholog group fasta file must be in a folder OGN/seqs/OGN.fasta
## if all fastas are in one location, use this script

# Extract all identifiers from filenames in the folder
#cd ogfolder
#temp_file="ogs_to_align.txt"
#touch $temp_file
#for fastaf in *.fasta; do
#    ogn=$(basename "$fastaf" .fast)
#    ogn >> "$temp_file"
#    mkdir -p $house/$ogn
#    mkdir -p $house/$ogn/seqs
#    cp $fastaf $house/$ogn/seqs/
#done



ogfile=$VSC_SCRATCH_VO_USER/ogs_to_align/ogs_to_align.txt


while IFS= read -r line; do
    echo "Processing line: $line"
    ogn=$line
    sdir=$house/$ogn
    cd $sdir

    nextflow run $VSC_SCRATCH_VO_USER/simsapiper/simsapiper.nf -resume\
        -profile minihydra \
        --data $sdir \
        --outName $ogn \
        --minimagic \
        --stopHyperconserved \
        --localModel 1 --retrieve false --model false \
        --outFolder $sdir 
    cd $house
done < $ogfile
