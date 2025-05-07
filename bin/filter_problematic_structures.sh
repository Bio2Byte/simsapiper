#!/bin/bash

set -euo pipefail

echo 'Starting filtering problematic seqs'

# Inputs
seqsToAlign="$1"
outFolder="$2"

echo $seqsToAlign
echo $outFolder
# Derived filenames
original="original_${seqsToAlign}"
filtered="fil_${seqsToAlign}"
removed="removed_${seqsToAlign}"
entries="entries_to_remove.txt"

# Step 1: Backup original
cp "$seqsToAlign" "$original"

# Step 2: Collect problematic sequence identifiers
> "$entries"

echo "$outFolder"/resources_*.txt

while IFS= read -r line; do
    if [[ "$line" == *"runTcoffee"* && "$line" == *"FAILED"* && "$line" == *"$seqsToAlign"* ]]; then
        wdir=$(echo "$line" | awk '{print $2}')
        log_file=$(find ../../"${wdir}"* -name ".command.log" -print -quit)
        echo "Processing $log_file"

        
        if grep -q "attempt 3" "$log_file"; then
            echo "Detected structure pair error"
            grep "attempt 3" "$log_file" | \
                    awk -F'\\[' '{print $2}' | awk -F']' '{print $1}' | \
                    tr ' ' '\n' | sort | uniq | \
                    sed 's/\.pdb$//' >> "$entries" 
    
        elif grep -q "^[0-9]\+ -- ERROR: UNSPECIFIED UNSPECIFIED" "$log_file"; then
            echo "Detected structure error"
            awk '
                /-- ERROR: UNSPECIFIED UNSPECIFIED/ { in_block = 1; next }
                /\*{5,}/ { in_block = 0 }
                in_block && /structures\/[^ ]+\.pdb/ {
                    for (i = 1; i <= NF; i++) {
                        if ($i ~ /structures\/[^ ]+\.pdb/) {
                            split($i, a, "/")
                            gsub(".pdb$", "", a[length(a)])
                            print a[length(a)]
                        }
                    }
                }
            ' "$log_file" >> "$entries"
        else
            echo "No known error pattern found in $log_file"
        fi
    fi
done < "$outFolder"/resources_*.txt


echo "problematic seqs found:"
cat "$entries"

echo "Remove duplicates"
sort -u "$entries" -o "$entries"

cat "$entries"
# Step 3: Filter FASTA file
awk 'BEGIN {
    while ((getline line < "'"$entries"'") > 0) remove[line] = 1
}
/^>/ {
    header = substr($0, 2)
    write = !(header in remove)
    write_removed = (header in remove)
}
{
    if (write) print > "'"$filtered"'"
    if (write_removed) print > "'"$removed"'"
}' "$seqsToAlign"

# Step 4: Replace original with filtered
cp "$filtered" "$seqsToAlign"

echo "Filtering complete."
echo "Original saved as: $original"
echo "Filtered written to: $seqsToAlign"
echo "Removed sequences written to: $removed"
