import sys
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def freq_DSSP_codes(dssp_file):
    read_dsspMSA = AlignIO.read(dssp_file, 'fasta')

    dsspMSA_array = np.array([list(rec.seq) for rec in read_dsspMSA])
    transposed = np.transpose(dsspMSA_array)

    codes_interest = {"Helix":["H","I","E"], "Sheet":["E","B"], "Loop":["X","T","S"]}

    tot_freq_SS = {}
    for SS,codes in codes_interest.items():
        freq = []
        for i, column in enumerate(transposed):
            tot = 0
            for code in codes:
                    tot +=list(column).count(code)
            freq.append(tot / len(column))
        tot_freq_SS[SS]=freq
    return tot_freq_SS



def plot_dssp_frequencies(tot_freq_SS, output_file="dssp_freq_plot.pdf"):
    length = len(next(iter(tot_freq_SS.values())))  # Number of alignment columns
    x = np.arange(1, length + 1)

    fig_width = min(max(length / 25, 12), 50)
    fig, axes = plt.subplots(3, 1, figsize=(fig_width, 6))

    ss_types = list(tot_freq_SS.keys())
    label_interval = 10

    colors = ["green","blue","orange"]
    for i, (ss_type,col) in enumerate(zip(ss_types,colors)):
        ax = axes[i]
        freq = tot_freq_SS[ss_type]
        ax.plot(x, freq, label=ss_type, linewidth=1,color=col)

        ax.set_xlim(0.5, length + 0.5)
        ax.set_ylim(-0.1,1.1)
        ax.set_ylabel("Frequency", fontsize=8)

        # Major ticks
        major_tick_positions = [i for i in range(1, length + 1) if i % label_interval == 0 or i == 1]
        ax.set_xticks(major_tick_positions)
        ax.set_xticklabels([str(i) for i in major_tick_positions], fontsize=6)

        # Minor ticks
        ax.set_xticks(np.arange(1, length + 1), minor=True)
        ax.tick_params(axis='x', which='major', length=4, color='gray', labelsize=6)
        ax.tick_params(axis='x', which='minor', length=2, color='lightgray')

        # ax.legend(loc='upper right', fontsize=8)
        ax.set_title(f"Secondary Structure: {ss_type}", fontsize=10)

    axes[-1].set_xlabel("Column position", fontsize=9)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


dssp_file = sys.argv[1]
output_folder=sys.argv[2] #"DSSPcodes_aligned.pdf"
total_frequencies_SS = freq_DSSP_codes(dssp_file)
plot_dssp_frequencies(total_frequencies_SS, 'DSSPcodes_aligned.pdf')

import csv
with open(output_folder+"/DSSPcodes_frequencies.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Position", "Helix residues", "Sheet residues","Loop residues"])
    for i, (helix, sheet, loop) in enumerate(zip(total_frequencies_SS["Helix"], total_frequencies_SS["Sheet"],total_frequencies_SS["Loop"]), start=1):
        writer.writerow([i, helix, sheet, loop])