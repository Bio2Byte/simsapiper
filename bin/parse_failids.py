
import sys
import glob

p=sys.argv[1]
#p='results/simsa_all_groel_80_enr_2024_01_30_15_12_49/collected_fail_ids.txt'

fulldatafiles=sys.argv[2]
#fulldata='data/seqs/'
f#ulldatafiles=glob.glob(fulldata+'*.fasta')

with open(p,'r') as pf:
    l=pf.read().splitlines()
    
l = [elem.lstrip() for elem in l]
print (l)



for f in fulldatafiles:
    print(f)
    with open(f,'r') as ff:
        fasta=ff.read().splitlines()

    corrseq=False
    outfasta=[]
    for elem in fasta:
        if corrseq==False:
            sid=elem.split('>')[-1].split('_')[0]
            result = list(filter(lambda x: sid == x, l))
            if result==[]:
                outfasta.append(elem) 
                corrseq=False
            else:
                corrseq=True
                print (elem)

        else:
            corrseq=False

    for nelem in outfasta:
        with open(f[:-6]+'_filtered.fasta', 'a+') as of:
            of.write(nelem+'\n')

    print('output file:', f[:-6]+'_filtered.fasta')
    print ("input seqs", int(len(fasta)/2-0.5), "filter:" ,len(l), 'ouput seqs:', int(len(outfasta)/2-0.5))

