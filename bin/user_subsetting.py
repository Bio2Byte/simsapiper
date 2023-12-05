import sys
import pandas as pd
from Bio import SeqIO

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
write_fasta_from_df(match_df, 'label', 'seq_f', outfilename)




"""


seq_in_dict = SeqIO.to_dict(SeqIO.parse(sequenceFastas, "fasta"))

for seq_record in SeqIO.parse(foundModels, "fasta"):
    if seq_record.id in seq_in_dict.keys():      
        print (seq_record)
        #use this line to write standard fasta record using SeqIO.write
        SeqIO.write(seq_record, outfilename, "fasta")
    else:
        print('la')




if newdict == {}:
    print ("nothing was matched!")
else:
    if 'orphan' in sequenceFastas:
        print ('orphan sequence from ',sequenceFastas )
        for id,seq in newdict.items():
            outfilename = 'orphan_' +id +'_matchedModel.fasta'
            with open (outfilename, 'w') as out:
                out.write('>'+id+'\n'+seq)
    else:
        print (outfilename)
        with open (outfilename, 'w') as out:
            for id in newdict :     
                SeqIO.write(id, out, "fasta-2line")

           # for id,seq in newdict.items():
            #    out.write('>'+id+'\n'+seq+'\n')






my_id_file     = open(foundModels,'r')
my_fasta_file  = open(sequenceFastas,'r')
result_file    = open(outfilename, "w")

my_dictionary = {} # fasta IDs are keys, value can be anything.



ogLabels = []
for record in SeqIO.parse(sequenceFastas, "fasta"):
    ogLabels.append(record.id)

records =  [''.join([ c if c.isalnum() else "_" for c in elem ]) for elem in ogLabels]
print (records)

for elem in records:
    my_dictionary[elem] = 'value'


for seq_record in SeqIO.parse(my_fasta_file, "fasta"):
    if seq_record.id in my_dictionary:      
        #use this line to write standard fasta record using SeqIO.write
        SeqIO.write(seq_record, result_file, "fasta")

my_fasta_file.close()
my_id_file.close()
result_file.close()



mydict = SeqIO.to_dict(SeqIO.parse(foundModels, "fasta"))
newdict={}

for rkey in records:
    try:
        nreq = mydict[rkey]
        newdict.setdefault(rkey, nreq)
    except:
        print(rkey +" is not in this file or no structure has been found")

print(newdict)
if newdict == {}:
    print ("nothing was matched!")
else:
    if 'orphan' in sequenceFastas:
        print ('orphan sequence from ',sequenceFastas )
        for id,seq in newdict.items():
            outfilename = 'orphan_' +id +'_matchedModel.fasta'
            with open (outfilename, 'w') as out:
                out.write('>'+id+'\n'+seq)
    else:
        print (outfilename)
        with open (outfilename, 'w') as out:
            for id in newdict :     
                SeqIO.write(id, out, "fasta-2line")

           # for id,seq in newdict.items():
            #    out.write('>'+id+'\n'+seq+'\n')



   # if seq_record.id in newdict:      
        #use this line to write standard fasta record using SeqIO.write
       



"""