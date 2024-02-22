
#Map DSSP annotations to MSA    
import sys
import pandas as pd 
from Bio.PDB.DSSP import make_dssp_dict
import numpy as np
import glob
from Bio import pairwise2, SeqIO


msa_file = sys.argv[1] #input file - finalMSA.fasta?
annotated_dssp_msa = sys.argv[2] #output file that will be written- finalMSA_dssp.fasta?

#dssp_files = glob.glob("dssp/*.dssp")
dssp_files = glob.glob(sys.argv[3]+"/*.dssp") 
dssp_folder = sys.argv[3]

def dssp_parse(dssp_file):
        dssp_tup = make_dssp_dict(dssp_file)
        dssp_dic = dssp_tup[0]

        df = pd.DataFrame.from_dict(dssp_dic, orient='index').reset_index(drop=True)
        df_filter = df.iloc[:,[0,1]]

        modelseq="".join(df_filter[0].values.tolist())
        
        dssp_orig=df_filter[1].values.tolist()
        dssp = "".join([elem.replace('-','X') for elem in dssp_orig])

        outtup= [modelseq, dssp]
        return outtup
    
def match_msa(dssp_parsed,msa_file):
    mdf = pd.read_csv(msa_file, header=None)
    msa_df = pd.DataFrame({'label':mdf[0].iloc[::2].values, 'seq':mdf[0].iloc[1::2].values})

    dssp_df = pd.DataFrame.from_dict(dssp_parsed, orient='index').reset_index()
    dssp_df = dssp_df.set_axis(['label', 'model', 'dssp'], axis=1)

    df = msa_df.merge(dssp_df, how='outer', on='label')

    print(df.iloc[[2]])
    return (df)

def handle_unmappable(df):
    # Iterate over rows where df.same == False
    for index, row in df[df['same'] == False].iterrows():
        print(row)
        fixed=False
        model_seq = row['model']
        unali_seq = row['unali']
        dssp = row['dssp']

        print(unali_seq)
        print(model_seq)
        print(dssp)

        if len(model_seq)==len(unali_seq): #unmatching is because of point mutations, we ignore point mutations
            df.at[index, 'same'] = True
        else: #either the sequence of entry of the sequence of the 3D structure/model have different length
            # align
            alignments = pairwise2.align.globalms(model_seq, unali_seq,2, -1, -10, -0.1)
            print(alignments)
            aligned_model_seq = alignments[0][0]
            aligned_unali_seq = alignments[0][1]
            # Count gaps
            gap_count_unali = aligned_unali_seq.count('-')
            gap_count_model = aligned_model_seq.count('-')

            if gap_count_unali == 0 and gap_count_model>1: #seq model is shorter than seq
                if aligned_model_seq[0] == "-" and aligned_model_seq[-1]!= "-": #longer at the front
                    for i in range(len(aligned_model_seq)-1):
                        if aligned_model_seq[i]=="-" and aligned_model_seq[i+1]!='-':
                            end_front_alignment = i+1
                            break
                    dssp=unali_seq[:end_front_alignment]+dssp
                    fixed = True
                elif aligned_model_seq[0] != "-" and aligned_model_seq[-1]== "-": #longer at the end 
                    for i in range(len(aligned_model_seq)-1):
                        if aligned_model_seq[i]!="-" and aligned_model_seq[i+1]=='-':
                            start_back_alignment = i+1
                            break
                    dssp=dssp+unali_seq[start_back_alignment:]
                    fixed = True
                elif aligned_model_seq[0] == "-" and aligned_model_seq[-1]== "-": #longer both at front and end of protein
                    for i in range(len(aligned_model_seq)-1):
                        if aligned_model_seq[i]=="-" and aligned_model_seq[i+1]!='-':
                            end_front_alignment = i+1
                        elif aligned_model_seq[i]!="-" and aligned_model_seq[i+1]=='-':
                            start_back_alignment = i+1
                            break
                    dssp = unali_seq[:end_front_alignment] + dssp + unali_seq[start_back_alignment:]
                    fixed = True
                if fixed:
                    df.at[index, 'dssp'] = dssp
                    df.at[index, 'same'] = True
            elif gap_count_unali > 1 and gap_count_model==0: #seq is shorter than seq model
                if aligned_unali_seq[0] == "-" and aligned_unali_seq[-1] != "-": #longer at the front
                    for i in range(len(aligned_unali_seq)-1):
                        if aligned_unali_seq[i]=="-" and aligned_unali_seq[i+1]!='-':
                            end_front_alignment = i+1
                            break
                    dssp=dssp[end_front_alignment:]
                    fixed = True
                elif aligned_unali_seq[0] != "-" and aligned_unali_seq[-1] == "-": #longer at the end 
                    for i in range(len(aligned_unali_seq)-1):
                        if aligned_unali_seq[i]!="-" and aligned_unali_seq[i+1]=='-':
                            start_back_alignment = i+1
                            break
                    dssp=dssp[:start_back_alignment]
                    fixed = True
                elif aligned_unali_seq[0] == "-" and aligned_unali_seq[-1]== "-": #longer both at front and end of protein
                    for i in range(len(aligned_unali_seq)-1):
                        if aligned_unali_seq[i]=="-" and aligned_unali_seq[i+1]!='-':
                            end_front_alignment = i+1
                        elif aligned_unali_seq[i]!="-" and aligned_unali_seq[i+1]=='-':
                            start_back_alignment = i+1
                            break
                    dssp = dssp[end_front_alignment:start_back_alignment]
                    fixed = True
                if fixed:
                    df.at[index, 'dssp'] = dssp
                    df.at[index, 'same'] = True


def map_msa_dssp(df):
    #cases that make modelseq != msaseq:
    #   -used esm fold, but sequence had X? (add more gaps)
    #   -beginning, end, loops not resolved
    #   -sg nuc mutation (easy)
    # we need to rethink that, this is a lot of work. 
    # probably align seqs and only map common regions? 
    # for now return fasta file with non-matching seqs and their dssp
    df['unali'] =  df.seq.str.replace('-','')

    df['same'] = np.where((df.model == df.seq.str.replace('-','')), True, False)

    df.to_csv("before_clean_dssp_csv.csv")

    handle_unmappable(df)

    df.to_csv("dssp_csv.csv")

    df_same = df[df.same==True]
    mapped_dssp = []
    for row in range(len(df_same)):
        dssp = df_same.iloc[row]['dssp']
        seq_msa =df_same.iloc[row]['seq']
        pos_dssp =0
        seq= ""
        for i in range(len(seq_msa)):
            if seq_msa[i] != '-' and pos_dssp <= len(dssp):
                seq = seq+ dssp[pos_dssp]
                pos_dssp+=1
            else:
                seq = seq +'-'
        mapped_dssp.append(seq)

    df_same['mapped_dssp'] = mapped_dssp
    df_same = df_same.filter(items=['label', 'mapped_dssp'])

    df_false = df[df.same==False]
    if len(df_false) >0: 
        error_df = df_false.filter(items=['label', 'dssp'])
        write_fasta_from_df(error_df,'unmappable_dssp')

        #for seqs without dssp match, just add the sequence back
        nomap=df_false.filter(items=['label','seq']).rename(columns={'label':'label', "seq": "mapped_dssp"})
        out_df= pd.concat([df_same,nomap])
        #drop seqs that have been collapsed with cd-hit
        out_df = out_df.dropna(subset = ['mapped_dssp'])
        #print (out_df)
    else:
        out_df = df_same

    write_fasta_from_df(out_df,annotated_dssp_msa)
    #print (out_df)

def write_fasta_from_df(df,outname):
    out_name = outname + '.fasta'

    df.iloc[:,0]= df.iloc[:,0] + '\n'

    rows = df.to_string(header=False,index=False,index_names=False).split('\n')

    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ','')
            m.write( row + '\n')


#dssp_list = dssp_files.split(' ')
dssp_parsed = {}
# for dssp_file in dssp_files:
for record in SeqIO.parse(msa_file,"fasta"):
    id=record.id
    # id = dssp_file.split(".")[0].split("/")[-1]
    dssp_parsed['>'+id] = dssp_parse(dssp_folder+"/"+id+".dssp")

df = match_msa(dssp_parsed,msa_file)
map_msa_dssp(df)