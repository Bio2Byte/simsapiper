
#Map DSSP annotations to MSA    
import sys
import pandas as pd 
from Bio.PDB.DSSP import make_dssp_dict
import numpy as np
import glob
from Bio import pairwise2, SeqIO


msa_file = sys.argv[1] #input file - finalMSA.fasta?
annotated_dssp_msa = sys.argv[2] #output file that will be written- finalMSA_dssp.fasta?

dssp_files = glob.glob("*/*.dssp")

def dssp_parse(dssp_file):
    """
    Parses a DSSP file and extracts sequence and secondary structure annotations.

    Parameters:
        dssp_file (str): Path to the DSSP file.

    Returns:
        list: A list containing two elements:
            - modelseq (str): The amino acid sequence from the 3D structure.
            - dssp (str): The corresponding DSSP secondary structure annotations,
                         with '-' replaced by 'X' to indicate unassigned elements.
    """
    dssp_tup = make_dssp_dict(dssp_file)
    dssp_dic = dssp_tup[0]

    df = pd.DataFrame.from_dict(dssp_dic, orient='index').reset_index(drop=True)
    df_filter = df.iloc[:, [0, 1]]

    modelseq = "".join(df_filter[0].values.tolist())
    dssp_orig = df_filter[1].values.tolist()
    dssp = "".join([elem.replace('-', 'X') for elem in dssp_orig])

    return [modelseq, dssp]

    
def match_msa(dssp_parsed, msa_file):
    """
    Matches DSSP-parsed structural annotations with sequences from MSA.

    Parameters:
        dssp_parsed (dict): A dictionary where keys are sequence labels and values contain model sequences and DSSP annotations.
        msa_file (str): Path to the MSA file in FASTA-like format.

    Returns:
        pandas.DataFrame: A merged DataFrame with columns:
            - 'label': Sequence identifiers
            - 'seq': Sequence from MSA
            - 'model': Sequence from 3D structure (if available)
            - 'dssp': DSSP annotation (if available)
    """
    # Read MSA file
    mdf = pd.read_csv(msa_file, header=None)

    msa_df = pd.DataFrame({
        'label': mdf[0].iloc[::2].values,
        'seq': mdf[0].iloc[1::2].values  # sequences are on alternating lines
    })

    # Convert DSSP-parsed dict to DataFrame
    dssp_df = pd.DataFrame.from_dict(dssp_parsed, orient='index').reset_index()
    dssp_df.columns = ['label', 'model', 'dssp']

    # Merge MSA and DSSP data on sequence label
    df = msa_df.merge(dssp_df, how='outer', on='label')
    return df

def handle_unmappable(df):
    """
    Attempts to fix mismatches between ungapped model sequences and aligned MSA sequences
    by allowing limited point mutations and aligning sequences with Biopython’s local alignment.

    Parameters:
        df (pandas.DataFrame): DataFrame containing at least the following columns:
            - 'label': sequence identifier
            - 'model': model (ungapped) sequence
            - 'unali': aligned sequence without gaps
            - 'dssp': DSSP annotation
            - 'same': boolean indicating perfect match (can be updated if fixed)

    Returns:
        pandas.DataFrame: Modified DataFrame with updated 'same' and potentially corrected 'dssp'.
    """

    for index, row in df[df['same'] == False].iterrows():
        model_seq = row['model']
        unali_seq = row['unali']
        dssp = row['dssp']
        label = row['label']
        fixed = False

        # Case 1: Allow small mutations (up to 5% difference)
        if len(model_seq) == len(unali_seq):
            mismatch_count = sum(1 for a, b in zip(model_seq, unali_seq) if a != b)
            if mismatch_count <= 0.05 * len(unali_seq):
                df.at[index, 'same'] = True
                continue
            else:
                print(f"[Warning] >5% mutations in {label} sequence VS sequence 3D structure. Check './msas/unmappable.txt'.")

        # Case 2: Use alignment to resolve differences in length
        alignments = pairwise2.align.localms(model_seq, unali_seq, 2, -1, -10, -0.1)
        aligned_model_seq, aligned_unali_seq = alignments[0][:2]

        # Fix Met alignment bug
        if aligned_unali_seq.startswith("M-"):
            aligned_unali_seq = "-M" + aligned_unali_seq[2:]
        elif aligned_model_seq.startswith("M-"):
            aligned_model_seq = "-M" + aligned_model_seq[2:]

        gap_model = aligned_model_seq.count('-')
        gap_unali = aligned_unali_seq.count('-')
        
        def find_front_index(seq):
            for i in range(len(seq)-1):
                if seq[i] == '-' and seq[i+1] != '-':
                    return i + 1
            return 0

        def find_back_index(seq):
            for i in range(len(seq)-1):
                if seq[i] != '-' and seq[i+1] == '-':
                    return i + 1
            return len(seq)

        # Handle model shorter than unali
        if gap_unali == 0 and gap_model >= 1:
            if aligned_model_seq.startswith('-') and not aligned_model_seq.endswith('-'):
                i = find_front_index(aligned_model_seq)
                dssp = 'X'*i + dssp
                fixed = True
            elif not aligned_model_seq.startswith('-') and aligned_model_seq.endswith('-'):
                i = find_back_index(aligned_model_seq)
                dssp = dssp + (len(aligned_unali_seq)-i-1)*'X'
                fixed = True
            elif aligned_model_seq.startswith('-') and aligned_model_seq.endswith('-'):
                front = find_front_index(aligned_model_seq)
                back = find_back_index(aligned_model_seq)
                dssp = 'X'*front + dssp + (len(aligned_unali_seq)-back-1)*'X'
                fixed = True

        # Handle unali shorter than model
        elif gap_unali >= 1 and gap_model == 0:
            if aligned_unali_seq.startswith('-') and not aligned_unali_seq.endswith('-'):
                i = find_front_index(aligned_unali_seq)
                dssp = dssp[i:]
                fixed = True
            elif not aligned_unali_seq.startswith('-') and aligned_unali_seq.endswith('-'):
                i = find_back_index(aligned_unali_seq)
                dssp = dssp[:i]
                fixed = True
            elif aligned_unali_seq.startswith('-') and aligned_unali_seq.endswith('-'):
                front = find_front_index(aligned_unali_seq)
                back = find_back_index(aligned_unali_seq)
                dssp = dssp[front:back]
                fixed = True

        # Handle gaps at front of one and gaps at end of the other one
        else:
            if aligned_model_seq.startswith('-'): #then seq must have '-' at the end otherwise would have entered previous condition
                front = find_front_index(aligned_model_seq)
                back = find_back_index(aligned_unali_seq)
                dssp = 'X'*front + dssp[:back]
                fixed = True

            if aligned_unali_seq.startswith('-'): #then 3D seq must have '-' at the end otherwise would have entered previous condition
                front = find_front_index(aligned_unali_seq)
                back = find_back_index(aligned_model_seq)
                dssp = dssp[front:]+ (len(aligned_unali_seq)-back-1)*'X'
                fixed = True

        if fixed:
            df.at[index, 'dssp'] = dssp
            df.at[index, 'same'] = True
        else: #i give up
            print(f"[Warning] Could not match sequence of {label} to its 3D model. Please check './msas/unmappable.txt'.")

    return df


def map_msa_dssp(df):
    """
    Aligns DSSP secondary structure annotations to an MSA-aligned sequence.

    Parameters:
        df (pandas.DataFrame): A DataFrame containing at least the following columns:
            - 'label': sequence identifier
            - 'seq': MSA-aligned sequence (with gaps)
            - 'model': ungapped model sequence from 3D structure
            - 'dssp': DSSP annotation for the model sequence

    Outputs:
        - Saves 'dssp_csv.csv' with mapped annotations
        - Writes two FASTA files:
            - 'annotated_dssp_msa.fasta': successfully mapped DSSP annotations
            - 'unmappable_dssp.fasta': sequences that couldn't be mapped
    """

    df['unali'] = df.seq.str.replace('-', '', regex=False)
    df['same'] = np.where(df['model'] == df['unali'], True, False)

    df.to_csv("dssp_csv_before_handling.csv", index=False)

    # Remove protein models without matching sequence
    df = df.dropna(subset=['seq'])

    # Handle sequences that may require additional fixing
    df = handle_unmappable(df)

    df.to_csv("dssp_csv.csv", index=False)

    mapped_dssp = []
    failed_indices = []
    for row in range(len(df)):
        if df.iloc[row]['same']:
            dssp = df.iloc[row]['dssp']
            seq_msa = df.iloc[row]['seq']
            pos_dssp = 0
            aligned_dssp = ""
            try:
                for char in seq_msa:
                    if char != '-' and pos_dssp < len(dssp):
                        aligned_dssp += dssp[pos_dssp]
                        pos_dssp += 1
                    else:
                        aligned_dssp += '-'
                mapped_dssp.append(aligned_dssp)
            except Exception:
                df.at[row, "same"] = False
                failed_indices.append(row)

    # Build DataFrame for successfully mapped sequences
    df_same = df[df.same == True].copy()
    df_same['mapped_dssp'] = mapped_dssp
    df_same = df_same[['label', 'mapped_dssp']]

    # Handle failed mappings
    df_false = df[df.same == False]
    if not df_false.empty:
        error_df = df_false[['label', 'unali', 'model', 'dssp']]
        write_fasta_from_df(error_df, 'unmappable_dssp')

        # # Reuse the input sequence as fallback
        #Wim doesn't like this feature
        # nomap = df_false[['label', 'seq']].rename(columns={"seq": "mapped_dssp"})
        # out_df = pd.concat([df_same, nomap])
        # out_df = out_df.dropna(subset=['mapped_dssp'])

    out_df = df_same
    write_fasta_from_df(out_df, annotated_dssp_msa)


def write_fasta_from_df(df, outname):
    """
    Writes a FASTA file from a DataFrame.

    Parameters:
        df (pandas.DataFrame): The data to write.
        outname (str): The name of the output file (without .fasta extension).
    """
    out_path = outname + '.fasta'
    lines = []

    if outname == "unmappable_dssp":
        for i in range(len(df)):
            label = df.iloc[i, 0]
            lines.extend([
                f"{label}_seq_inputted\n{df.iloc[i, 1]}",
                f"{label}_seq_model\n{df.iloc[i, 2]}",
                f"{label}_dssp\n{df.iloc[i, 3]}"
            ])
    else:
        for _, row in df.iterrows():
            label, seq = row[0], row[1]
            lines.append(f"{label}\n{seq}")

    with open(out_path, 'w') as f:
        for line in lines:
            f.write(line.strip() + '\n')


#Extract info from dssp files and sequence 3D structure
dssp_parsed = {}
for dssp_file in dssp_files: 
    id = dssp_file.split(".")[0].split("/")[-1]
    dssp_parsed['>'+id] = dssp_parse(dssp_file)

#Extract aligned sequence and store with sequence 3D structure and DSSP sequence
df = match_msa(dssp_parsed,msa_file)

#Check if 3D structure sequence and AA sequence match, convert DSSP sequence to an aligned sequence.
map_msa_dssp(df)