import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import sys

dsspfile=sys.argv[1]#'dssp_squeezed_merged_magicMSA.fasta'
outname='DSSPcodes_sequenceView'

def read_fasta_to_csv(seq_file):
    df = pd.read_csv(seq_file, header=None)
    seq_df = pd.DataFrame({'label':df[0].iloc[::2].values, 'msa':df[0].iloc[1::2].values})

    seq_df.label= seq_df.label.str.lstrip(' ')
    seq_df.label= seq_df.label.str.lstrip('>')
    return seq_df

def plot_dssp_ali(df,col,outname):

    pdf=df[['label',col]]
    pdf.set_index('label',inplace=True, drop=True)
    pdf=pdf[col].str.split('',expand=True)
    pdf = pdf.iloc[:, 1:-1]  #remove inserted /final col


    codes_interest = {"Helix":["H","I","G"], "Sheet":["E","B"], "Loop":["X","T","S"], 'Gap':['-']}
    colordic={'Helix':'green', 'Sheet':'blue','Loop':'darkorange','Gap':'lightgrey'}

    symbol_to_cat = {s: cat for cat, symbols in codes_interest.items() for s in symbols}
    cat_to_int = {cat: i for i, cat in enumerate(codes_interest.keys())}
    int_to_color = [colordic[cat] for cat in codes_interest.keys()]

    matrix = pdf.replace(symbol_to_cat).replace(cat_to_int).to_numpy()

    length = pdf.shape[1]

    fig_width = min(max(length / 25, 12), 80)
    fig_len = min(max(pdf.shape[0] / 25, 5), 80)

    fig, ax = plt.subplots(figsize=(fig_width, fig_len))

    cmap = ListedColormap(int_to_color)

    im = ax.imshow(matrix, aspect="auto", cmap=cmap)

    # Y axis labels = sequence IDs
    ax.set_yticks(range(len(pdf)))
    ax.set_yticklabels(pdf.index)

    ax.set_xlabel("Column Position")
    label_interval = 10
    # Major ticks
    major_tick_positions = [i for i in range(1, length + 1) if i % label_interval == 0 or i == 1]
    ax.set_xticks(major_tick_positions)
    ax.set_xticklabels([str(i) for i in major_tick_positions], fontsize=6)

    # Minor ticks
    ax.set_xticks(range(1, length + 1), minor=True)
    ax.tick_params(axis='x', which='major', length=4, color='gray', labelsize=6)
    ax.tick_params(axis='x', which='minor', length=2, color='lightgray')

    # Add legend manually
    handles = [plt.Rectangle((0,0),1,1, color=colordic[cat]) for cat in codes_interest.keys()]
    ax.legend(handles, codes_interest.keys(), bbox_to_anchor=(1.05, 1), loc='upper left')


    ax.set_title("Secondary Structure Alignment")

    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.close()

df=read_fasta_to_csv(dsspfile)
plot_dssp_ali(df, 'msa', outname+'.pdf')
