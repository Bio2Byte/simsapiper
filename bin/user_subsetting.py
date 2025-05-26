import sys
import pandas as pd

foundModels = sys.argv[1] #individual fasta files
sequenceFastas = sys.argv[2] #one fasta file with all structure seqs

outfilename= foundModels[:-6] +'_matchedModel'


def write_fasta_from_df(df,labelcol,seqcol,outname):
    out_name = outname + '.fasta'
    out_df = df.filter(items=[labelcol, seqcol])
    out_df[labelcol] =  out_df[labelcol]+"\n"
    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')
    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ',"")
            m.write( '>'+row + '\n')


found_df = pd.read_csv(foundModels, header=None)
found_df = pd.DataFrame({'label':found_df[0].iloc[::2].values, 'seq_f':found_df[0].iloc[1::2].values})

found_df.label= found_df.label.str.lstrip('>')


seq_df = pd.read_csv(sequenceFastas, header=None)
seq_df = pd.DataFrame({'label':seq_df[0].iloc[::2].values, 'seq_s':seq_df[0].iloc[1::2].values})

seq_df.label= seq_df.label.str.lstrip(' ')
seq_df.label= seq_df.label.str.lstrip('>')
seq_df.label = seq_df.label.str.replace('[^a-zA-Z0-9]', '_')

match_df = seq_df.merge(found_df, on='label')
print(match_df)
match_df.drop_duplicates(inplace=True)

print(match_df)
write_fasta_from_df(match_df, 'label', 'seq_f', outfilename)
