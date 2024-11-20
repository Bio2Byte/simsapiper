import sys
import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt 

msa_file = sys.argv[1]
output_csv = "pairwise_identity_matrix.csv" 
output_csv_avg = "average_identity_per_sequence.csv"
histogram_plot = "identity_histogram.png"
histogram_plot_av = "average_identity_histogram.png"

# Read file and split into labels and sequences
found_df = pd.read_csv(msa_file, header=None)
found_df = pd.DataFrame({'label': found_df[0].iloc[::2].values, 'seq': found_df[0].iloc[1::2].values})
seq_df = pd.DataFrame(found_df.seq.apply(list).tolist())

# Fully conserved positions
a = seq_df.to_numpy()
common = np.all(a == a[0], axis=0) 
common_counter = np.sum(common)

# Fully occupied (gapless) positions
gaps = (seq_df == '-').sum()
occupied = (gaps == 0).sum() 

# General sequence identity calculation with matrix output
labels = found_df['label'].values
seqlist = found_df.seq.tolist()
n = len(seqlist)
identity_matrix = np.zeros((n, n))  

for i, seq1 in enumerate(seqlist):
    for j, seq2 in enumerate(seqlist):
        if i < j: 
            non_gap_positions = [(x != '-') and (y != '-') for x, y in zip(seq1, seq2)]
            matches = sum(x == y for x, y in zip(seq1, seq2) if x != '-' and y != '-')
            non_gap_count = sum(non_gap_positions)
            identity = matches / non_gap_count if non_gap_count > 0 else 0
            identity_matrix[i, j] = identity
            identity_matrix[j, i] = identity 

identity_matrix = np.round(identity_matrix * 100, 2) 

# Create a DataFrame for the matrix and save as CSV
identity_df = pd.DataFrame(identity_matrix, index=labels, columns=labels) 
identity_df.to_csv(output_csv)

# Calculate average sequence identity per sequence and save to CSV
average_identities = identity_matrix.mean(axis=1)  # Average per row (sequence)
average_identities_df = pd.DataFrame({"Label": labels, "Average_Identity (%)": average_identities})
average_identities_df.to_csv(output_csv_avg, index=False)


# Calculate and print summary information
av_id = np.round(identity_matrix[np.triu_indices(n, 1)].mean(), 2) 

line1 = f"* The final alignment contains {common_counter} conserved positions and {occupied} gapless positions of {seq_df.shape[1]} total positions."
line2 = f"* The average pairwise sequence identity in the final alignment is {av_id}%. Gaps were treated as mismatches."

print(line1)
print(line2)


# Generate histogram plot for pairwise sequence identities
pairwise_identities = identity_matrix[np.triu_indices(n, 1)]  # Extract upper triangle values only
plt.hist(pairwise_identities, bins=range(0, 105, 5), edgecolor='black')  # Bins in 5% increments
plt.xlabel("Pairwise Sequence Identity (%)")
plt.ylabel("Frequency")
plt.title("Distribution of Pairwise Sequence Identities")
plt.savefig(histogram_plot)
plt.close()

# Generate histogram plot for average pairwise sequence identities
plt.hist(average_identities_df["Average_Identity (%)"], bins=range(0, 101, 1), edgecolor='black')  # Bins in 1% increments
plt.xlabel("Average Pairwise Sequence Identity (%)")
plt.ylabel("Frequency")
plt.title("Distribution of Average Pairwise Sequence Identities")
plt.savefig(histogram_plot_av)
plt.close()