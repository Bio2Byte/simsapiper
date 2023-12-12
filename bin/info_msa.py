import sys
import pandas as pd 
import Levenshtein
import itertools


msa_file = sys.argv[1]
out_file = 'simsapiper_summary.md' #append line to file


found_df = pd.read_csv(msa_file, header=None)
found_df = pd.DataFrame({'label':found_df[0].iloc[::2].values, 'seq':found_df[0].iloc[1::2].values})
seq_df = pd.DataFrame(found_df.seq.apply(list).tolist())

#fully conserved positions
a = seq_df.to_numpy()
common = (a[0] == a).all(0)
common_counter = (common==True).sum() 

#fully occupied positions
nseqs, positions = seq_df.shape
gaps = (seq_df == '-').sum()
occupied = 0
for i in range (0,positions):
    if gaps[i] == 0:
        occupied+=1


#general sequence similarity
lev_ratios=[]
seqlist=found_df.seq.tolist()
combis = itertools.combinations(seqlist, 2)
for elem in combis:
    seq1,seq2 = elem
    lessgap_seq1=''
    lessgap_seq2=''
    common_gap=0
    for n in range(0,len(seq1)):
        if seq1[n] == '-':
            if seq2[n] !='-':
                lessgap_seq2+=seq2[n]
                lessgap_seq1+=seq1[n]
            else:
                common_gap+=1
        else:
            lessgap_seq2+=seq2[n]
            lessgap_seq1+=seq1[n]
    
    if len(seq2)-len(lessgap_seq2) !=common_gap:
        raise ValueError ("common gap removal failed")
  
    rat = Levenshtein.ratio(lessgap_seq1,lessgap_seq2)
    lev_ratios.append(rat)

av_id = round(100*sum(lev_ratios)/len(lev_ratios) ,2)
print (av_id)

with open (out_file, 'a') as f:
    line = 'The final alignment contains' +str(common_counter)+ 'conserved positions and ' +str(occupied)+ 'gapless positions of '+str(positions)+ 'total positions.'
    line2 = 'The average pairwise sequence identity (Levenshtein) in the final alignment is' +str(av_id)
    f.write(line+'\n'+line2)