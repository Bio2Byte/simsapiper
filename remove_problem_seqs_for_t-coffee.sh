#!/bin/bash

log_filepath="$1"
log_file="${1}*.nflog"
simsapath="$2"

echo logfilepath: $1
echo pipelinepath: $2

tail -n 30 log_file


echo 'Find failed T-coffee runs'

# Extract information from the log file using regex
while read -r line; do
    if [[ $line == *"terminated with an error exit status (1) -- Execution is retried (1)"* ]]; then
        # Extract folder path and file name
        partial_folder_path=$(echo "$line" | grep -oE '\[.*\]' | tr -d '[]')
        file_name=$(echo "$line" |  grep -oE '\(.*fasta\)' | tr -d '()') 

        echo "$file_name at $partial_folder_path"

        IFS='/' read -ra folders <<< "$partial_folder_path"
       

        # Search for the tcoffee_*.log file in the ../../work/ directory
        folder_path=$(find ../../work/${folders[0]} -type d -name "${folders[1]}*" -print -quit)       
        echo $folder_path


        # Copy tcoffee_*.log file to the current location
        if [ -n "$folder_path" ]; then
            cp $folder_path/.command.log ./${file_name}_${folders[0]}_${folders[1]}_tcoffee.fail
            echo "Copied tcoffee log from $folder_path to current location."
        else
            echo "Error: Folder with name $folder_path not found in ../../work/."
        fi
 

    fi
done < "$log_file"


echo "Collect problematic sequence identifiers"

for file in *fail; do
    while read -r line; do
            if [[ $line == *".pl structures/"* ]]; then
            result=$(echo "$line" | grep -oE 'structures/[^.]+\.pdb')

            #match structures/A0A7S9PST3_E__Epichloe_festucae.pdb : /*_
            echo "$result" | grep -oE '\/([^_]+)_' | sed 's/\// /; s/_$//' >> collected_fail_ids.txt
        fi
    done < $file
done

echo "Create filtered sequence file, please rerun SIMSApiper with this"

#use bin/parse_failids.py to filter them out of the sequence file 
python3 ${2}bin/parse_failids.py collected_fail_ids.txt seqs/seqs_to_align.fasta 


