
import math
from Bio import AlignIO
import sys
#Shannon_entropy calculates the conservation of the amino acids at each position. It takes into account
#the gaps in the columns. Check https://doi.org/10.1002/prot.10146 for info formula

msa_filename_noDesigned = sys.argv[1]

def shannon_entropy(list_input):
    tot = len(list_input) #total number of AA at particular position in MSA
    gaps = list_input.count("-") / tot #count frequency of gaps in that column of MSA
    unique_base = set(list_input) #remove duplicates, "-" is seen as an amino acid type
    unique_base_len = len(unique_base) #total number of AA at particular position in MSA with no duplicates
    entropy_list = [] # entropy of AA at particular position
    for base in unique_base:
        n_i = list_input.count(base)
        P_i = n_i / tot
        entropy_i = P_i * (math.log(P_i, 2))
        entropy_list.append(entropy_i)
    #sum entropy of every residue at a position and normalize it so it is between 0 and 1
    entropy_sum = math.fsum(entropy_list)
    if unique_base_len == 1: # log(1, 2) = 0; n/0 throws ZeroDivisionError
        shannon_entropy = entropy_sum
    else:
        unique_base_len_log = math.log(unique_base_len, 2)
        shannon_entropy = (-1 / unique_base_len_log) * entropy_sum
    #Return entropy AA and entropy gaps at 1 position
    #If entropy is high than there are many possible arrangements (high variability)
    return shannon_entropy, gaps
def conservation(alignment_file):
    conservation_AA_list = []
    for col_no in range(len(list(alignment_file[0]))):
        list_input = list(alignment_file[:, col_no])
        sh_entropy_AA, sh_entropy_gaps = shannon_entropy(list_input)
        #Translate entropy into conservation. 0 means no conservation, 1 highly conserved
        #As e took into account the gaps: if only gaps: entropy=0 freq_gaps=1 conservation = 0
        conservation_AA = (1 - sh_entropy_AA) * (1 - sh_entropy_gaps)
        conservation_AA_list.append(conservation_AA)
    return conservation_AA_list
#Retrieve sequence conservation (Shannon entropy)
msa_type = "fasta"
alignment_file_natural = AlignIO.read(msa_filename_noDesigned, msa_type)
conservation_AA_list = conservation(alignment_file_natural)


average = (sum(conservation_AA_list) / len(conservation_AA_list))*100


sys.stdout.write(str(round(average,4)))
sys.exit(0)



"""
Alternative:totalCount is number of seqs, aaCount the count for one amino acid type.

value = 0.0
totalCount = sum(aaList) * 1.0
for aaCount in aaList:
  if aaCount:
    pj = aaCount / totalCount
    value += pj * math.log(pj) # Is natural logarithm ln by default!

"""