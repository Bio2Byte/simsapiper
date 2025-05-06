#!/bin/bash

set -euo pipefail

# Inputs
seqsToAlign="$1"
outFolder="$2"

# Derived filenames
original="original_${seqsToAlign}"
filtered="fil_${seqsToAlign}"
removed="removed_${seqsToAlign}"
entries="entries_to_remove.txt"

# Step 1: Backup original
cp "$seqsToAlign" "$original"

# Step 2: Collect problematic sequence identifiers
> "$entries"

for resource in "$outFolder"/resources_*.txt; do
    grep "runTcoffee" "$resource" | grep "FAILED" | grep "$seqsToAlign" | awk '{print $2}' | while read -r wdir; do
        log_file=$(find ../../"${wdir}"* -name ".command.log" -print -quit)
        if [[ -f "$log_file" ]]; then
            echo "Processing $log_file"

            # Attempt 3 failures
            grep "attempt 3" "$log_file" | awk -F'\\[' '{print $2}' | awk -F']' '{print $1}' | tr ' ' '\n' >> "$entries"

            # UNSPECIFIED structure errors (embedded awk)
            awk '
                /-- ERROR: UNSPECIFIED UNSPECIFIED/ { in_block = 1; next }
                /\*{5,}/ { in_block = 0 }
                in_block && /structures\/[^ ]+\.pdb/ {
                    for (i = 1; i <= NF; i++) {
                        if ($i ~ /structures\/[^ ]+\.pdb/) {
                            split($i, a, "/")
                            print a[length(a)]
                        }
                    }
                }
            ' "$log_file" >> "$entries"
        fi
    done
done

# Remove duplicates
sort -u "$entries" -o "$entries"

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
