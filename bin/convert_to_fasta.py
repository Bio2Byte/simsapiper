from Bio import SeqIO
import sys

sequenceFiles = sys.argv[1]
format = sys.argv[2]
outname = sys.argv[3] + '.fasta'


message=   'converting %s with format %s into %s' %(sequenceFiles, format , outname)
print (message)


count = SeqIO.convert(sequenceFiles, format, outname, "fasta-2line")
print("Converted %i records" % count)