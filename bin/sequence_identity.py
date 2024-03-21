import sys
import pandas as pd 
import itertools


msa_file = sys.argv[1]


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


#general sequence id
id_ratios=[]
seqlist=found_df.seq.tolist()
combis = itertools.combinations(seqlist, 2)
for elem in combis:
    seq1,seq2 = elem
    lenNoCommonGap=0
    matchCounter=0
    for n in range(0,len(seq1)):
        if seq1[n]==seq2[n]:
            if seq1[n] != '-':
                lenNoCommonGap +=1
                matchCounter+=1
        else:
            lenNoCommonGap+=1
    rat=matchCounter/lenNoCommonGap
    id_ratios.append(rat)

av_id = round(100*sum(id_ratios)/len(id_ratios) ,2)
#print (av_id)

line = '* The final alignment contains ' +str(common_counter)+ ' conserved positions and ' +str(occupied)+ ' gapless positions of '+str(positions)+ ' total positions.'
line2 = '* The average pairwise sequence identity in the final alignment is ' +str(av_id) + '%. Gaps were treated as missmatch.'

print(line)
print(line2)
    