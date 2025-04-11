import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib

def shannon_entropy(column):
    """
    Compute normalized Shannon entropy and gap frequency for a column.
    Gaps are excluded from entropy calculation, but included in gap frequency.
    """
    total = len(column)
    gap_count = column.count("-")
    gap_freq = gap_count / total

    residues = [res for res in column if res != "-"]
    if not residues:
        return 0.0, 1.0  # All gaps, no informative residues

    total_residues = len(residues)
    residue_types = set(residues)

    entropy = 0.0
    for res in residue_types:
        freq = residues.count(res) / total_residues
        entropy -= freq * math.log(freq, 2)

    M = len(residue_types)
    if M > 1:
        entropy /= math.log(M, 2)  # Normalize to [0, 1]

    return entropy, gap_freq

def compute_conservation(alignment):
    """
    Calculate conservation score per position: 
    C(x) = (1 - E(x)) * (1 - G(x))
    """
    alignment_length = alignment.get_alignment_length()
    return [
        (1 - entropy) * (1 - gap)
        for entropy, gap in (shannon_entropy(list(alignment[:, i])) for i in range(alignment_length))
    ]

def compute_occupancy(alignment):
    """
    Compute occupancy (non-gap frequency) for each position.
    """
    alignment_length = alignment.get_alignment_length()
    total_sequences = len(alignment)
    counter = np.zeros(alignment_length)

    for record in alignment:
        for i, residue in enumerate(record.seq):
            if residue != "-":
                counter[i] += 1

    return counter / total_sequences

def plot_conservation(conservation_scores, occupancy, output_file):
    """
    Plot sequence occupancy and conservation (from entropy).
    Auto-adjusts label density for long MSAs (ticks every residue, labels downsampled).
    """
    length = len(conservation_scores)
    x = np.arange(1, length + 1)

    fig_width = min(max(length / 25, 12), 50)
    fig, ax = plt.subplots(figsize=(fig_width, 2.5))
    
    ax.set_title('Alignment Occupancy and Conservation')
    ax.plot(x, occupancy, label="Occupancy", color="black", linewidth=0.75)
    ax.set_ylabel("Sequence\occupancy")
    ax.set_xlabel("Column position")

    extent = [0.5, length + 0.5, 0, 1.05]
    cmap = matplotlib.cm.Blues(np.linspace(0, 1, 20))
    cmap = matplotlib.colors.ListedColormap(cmap[:5])

    img = ax.imshow(
        np.array(conservation_scores)[np.newaxis, :],
        cmap=cmap,
        extent=extent,
        aspect="auto",
        vmin=0,
        vmax=1,
    )
    
    label_interval = 10

    # Define major tick positions (for labels)
    major_tick_positions = [i for i in range(1, length + 1) if i % label_interval == 0 or i == 1]
    ax.set_xticks(major_tick_positions)
    ax.set_xticklabels([str(i) for i in major_tick_positions], fontsize=6)

    # Set minor ticks at every position (for alignment)
    ax.set_xticks(np.arange(1, length + 1), minor=True)

    # Style major and minor ticks separately
    ax.tick_params(axis='x', which='major', length=4, color='gray', labelsize=6)
    ax.tick_params(axis='x', which='minor', length=2, color='lightgray')


    # Add vertical grid at label positions
    #for xpos in range(1, length + 1):
    #    if xpos % label_interval == 0 or xpos == 1:
    #        ax.axvline(x=xpos, color='gray', linestyle='--', linewidth=0.3)

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=0.2, pad=0.3)
    cbar = plt.colorbar(img, cax=cax, orientation="vertical", ticks=[0, 1])
    cbar.ax.set_yticklabels(["No\nconservation", "Conserved"], fontsize=8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

msa_filename = sys.argv[1]
output_plot=sys.argv[2] #"ShannonEntropyConservation_occupency.pdf"
output_folder=sys.argv[3]

alignment = AlignIO.read(msa_filename, "fasta")
conservation_scores = compute_conservation(alignment)
occupancy_scores = compute_occupancy(alignment)

plot_conservation(conservation_scores, occupancy_scores,output_plot)
