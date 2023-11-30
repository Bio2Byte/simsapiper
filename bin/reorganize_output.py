# if the user wants an organized output file, they will need to fill in the following flag with the order in which 
# the cluster files will be added to the finalMSA file
# --order cluster3.fasta, cluster4.fasta, cluster5.fasta
# if 1 or more cluster sequence files are missing, they will be added after these files

#from Bio import SeqIO
import os
import sys
import pandas as pd

input_seqs = sys.argv[1].split(' ')
order_files_by= sys.argv[2]
msa_file = sys.argv[3]
outfile =sys.argv[4]+'.fasta'


if order_files_by != 'true':
    order_list = order_files_by.split(',')
    order_files_by_list=[]
    #need to be fasta2line formatted > use converted files
    for elem in order_list:
        nel = elem + '_converted.fasta'
        order_files_by_list.append(nel)
    #and imagine there is also a cluster 1_1 and 1_2 that didn't got mentionned in the command line so they need to be added after these 3 files
    for filename in input_seqs:
        if filename not in order_files_by_list:
            order_files_by_list.append(filename)
else:
    #order by input files (alphabetically)
    order_files_by_list = sorted(input_seqs)

print (order_files_by_list)

#create dataframe with all input sequences in correct order
order_by_df = pd.DataFrame()
for seqfile in order_files_by_list:
    print ('This step fails if there are "," in your original sequence labels, please remove them')
    sdf = pd.read_csv(seqfile, header=None)
    seq_df = pd.DataFrame({'label':sdf[0].iloc[::2].values, 'seq':sdf[0].iloc[1::2].values})
    
    seq_df.label = seq_df.label.str.replace('>','').str.replace('[^a-zA-Z0-9]', '_')

    if order_by_df.empty:
        order_by_df = seq_df
    else:
        order_by_df = pd.concat([order_by_df,seq_df])

#match with msa file, keep order of input files
mdf =  pd.read_csv(msa_file, header=None)
msa_df = pd.DataFrame({'label':mdf[0].iloc[::2].values, 'seq':mdf[0].iloc[1::2].values})
msa_df.label = msa_df.label.str.replace('>','')

reorga_df = order_by_df.merge(msa_df, how='left', on='label')
reorga_df = reorga_df.filter(items=['label','seq_y'])

out_df = reorga_df[reorga_df.seq_y.notna()]
print (out_df)


#Write the new organized MSA
out_df.iloc[:,0]= '>'+ out_df.iloc[:,0] + '\n'
rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')

with open (outfile, 'w') as m:
    for row in rows:
        row = row.replace('\\n','\n').replace(' ','')
        m.write( row + '\n')