#This code gathers all the steps you have to run in bash to run MAFFT merge locally (and avoids running ruby)
#See their "documentation" https://mafft.cbrc.jp/alignment/software/merge.html

from Bio import SeqIO
import sys

input_files_list = sys.argv[1].split(' ')

subclasses = []
#subclasses = ["alpha.fasta","beta.fasta","gamma.fasta","delta.fasta"] #can also be 1 file in list if no subclassing needed but structureless_seqs present so MAFFT needed
#non-aligned files need to be at the end of the list!
tail = []
for elem in input_files_list:
    if 'structureless_seqs' in elem:
        #file with all the structureless_seqs - unaligned
        structureless_seqs = elem
        tail.append(elem)
    elif 'orphan'in elem:
        tail.append(elem)
    elif 'removed'in elem:
        tail.append(elem)
    else:
        subclasses.append(elem)

all_files = subclasses + tail
print (all_files)
#add white line at the end of EVERY file
for myfile in all_files:
    with open(myfile, "a+") as f:
        f.write("\n")

#create input file
input_file = "input"
with open(input_file, 'w') as outfile:
    for fname in all_files:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

#create table: 1 line per file, start at 1 and per sequence add +1, continue incrementing so the final number 
# of your final sequence of your final file is equal to the total number of sequences you have in your entire dataset (without the structureless_seqs)
table = "tableMAFFT.txt"
lengths={}
#retrieve lengths files (the info from the structureless_seqs file doesn't need to be in tableMAFFT accoding to their documentation)
for filesub in subclasses:
    ids=0
    for record in SeqIO.parse(filesub, "fasta"):
        ids +=1
    lengths[filesub]=ids

start =0
with open(table, 'w') as f:
    for subfile in subclasses:
        print('submsa:',subfile)
        length = lengths[subfile]
        for i in range(start,length+start):
            f.write(str(i+1)+" ")
        f.write("#"+subfile+"\n")
        start = i+1