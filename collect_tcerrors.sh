#!/bin/bash

log_file="$1"
# Check if the log file exists
if [ ! -f "$log_file" ]; then
    echo "Log file $log_file not found."
    exit 1
fi


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
            cp $folder_path/tcoffee_* ./${file_name}_${folders[0]}_${folders[1]}_tcoffee.fail
            echo "Copied tcoffee log from $folder_path to current location."
        else
            echo "Error: Folder with name $folder_path not found in ../../work/."
        fi
 

    fi
done < "$log_file"

