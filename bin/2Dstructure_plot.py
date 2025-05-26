import sys
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import secstructartist as ssa


def limits_secondary_structure(dssp_file, conserved_2structure_dssp, anchor_point):
    read_dsspMSA = AlignIO.read(dssp_file, 'fasta')
    alignment_length = read_dsspMSA.get_alignment_length()

    dsspMSA_array = np.array([list(rec.seq) for rec in read_dsspMSA])
    transposed = np.transpose(dsspMSA_array)

    structure_elements = []
    last_structured = 0

    for conserved in conserved_2structure_dssp:
        positions = []
        for i, column in enumerate(transposed):
            freq = list(column).count(conserved) / len(column)
            if freq >= anchor_point:
                positions.append(i)

        if positions:
            start = positions[0]
            stop = None
            counter = 1
            for i in range(1, len(positions) - 1):
                if positions[i] - positions[i - 1] != 1:
                    stop = positions[i - 1]
                    structure_elements.append((conserved + str(counter), (start, stop)))
                    counter += 1
                    start = positions[i]
            structure_elements.append((conserved + str(counter), (start, positions[-1])))
            if last_structured < positions[-1]:
                last_structured = positions[-1]

    structure_elements = sorted(structure_elements, key=lambda x: x[1][0])
    print(structure_elements)
    dssp_code_plotss = {"E": "S", "H": "H", "G": "H", "I": "H" }
    secstruct_str = ""
    first = True
    if structure_elements == []:
        raise ValueError('No shared secondary structure elements could be identified')
    for ss in structure_elements:
        ss_name = ss[0][0]
        start, end = ss[1]
        plot_symbol = dssp_code_plotss[ss_name]
        length = end - start + 1

        if first:
            secstruct_str += "L" * start
            first = False
        else:
            loop_len = start - stop_previous - 1
            secstruct_str += "L" * loop_len

        secstruct_str += plot_symbol * length
        stop_previous = end
    secstruct_str += "L" * (alignment_length - stop_previous - 1)

    return secstruct_str


def plot_secondary_structure(secstruct_str, output_file="ConsensusSecondaryStructure_alignment.pdf"):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    length = len(secstruct_str)
    x = np.arange(1, length + 1)

    fig_width = min(max(length / 25, 12), 80)
    fig, ax = plt.subplots(figsize=(fig_width, 1.5))
    ax.set_title("Secondary Structure Alignment")
    ax.margins(0)

    # Draw secondary structure
    artist = ssa.SecStructArtist()
    artist["H"].fillcolor = (.9, 0., 0.)
    artist["S"].fillcolor = "#dc0"
    artist.draw(secstruct_str, x, ypos=0, ax=ax)

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


    # Remove Y axis
    ax.set_yticks([])
    ax.set_xlabel("Column Position")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

if __name__ == "__main__":
    dssp_file = sys.argv[1]
    squeeze =sys.argv[2]
    if squeeze == "true": #user hasn't specified towards what should be squeezed
        squeeze = "H,E"    
    if squeeze == "false": #squeeze is optional
        squeeze = "H,E"
    conserved_2structure_dssp = squeeze.split(',') # ["H","E"]
    output_folder = sys.argv[3]
    anchor_point=0.5 #to show consensus secondary structure
    secstruct_str = limits_secondary_structure(dssp_file,conserved_2structure_dssp, anchor_point)
    plot_secondary_structure(secstruct_str)
    
    # Save to CSV
    import csv
    with open(output_folder+"/SecondaryStructure_rawdata.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position", "consensus_SecondaryStructure"])
        for idx, ss in enumerate(secstruct_str, start=1):
            writer.writerow([idx, ss])
