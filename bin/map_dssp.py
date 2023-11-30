
#Map DSSP annotations to MSA    
import sys
import pandas as pd 
from Bio.PDB.DSSP import make_dssp_dict
import numpy as np
import glob

#dssp_files = sys.argv[1]
dssp_files = glob.glob("dssp/*.dssp")

msa_file = sys.argv[1] #input file - finalMSA.fasta?
annotated_dssp_msa = sys.argv[2] #output file that will be written- finalMSA_dssp.fasta?



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
        print (out_df)
    else:
        out_df = df_same

    write_fasta_from_df(out_df,annotated_dssp_msa)
    print (out_df)

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
for dssp_file in dssp_files:
    id = dssp_file.split(".")[0].split("/")[1]
    if len(id.rsplit("_",1)) >1:
        if 'active' in id.rsplit("_",1)[1] : #remove _active or _inactive (only for AFms) 
            id = id.rsplit("_",1)[0]
    
    dssp_parsed['>'+id] = dssp_parse(dssp_file)

df = match_msa(dssp_parsed,msa_file)
map_msa_dssp(df)