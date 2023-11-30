from Bio import AlignIO, SeqIO
import numpy as np
import sys

dssp_filter= ['H','B','E','G','I','T','S','-','X']

input_file = sys.argv[1] #final_MSA.fasta
dssp_file = sys.argv[2] #final_MSA_dssp.fasta = MSA with dssp code instead of residues
squeeze = sys.argv[3]
conserved_2structure_dssp = squeeze.split(',') # ["H","E"]

anchor_point = int(sys.argv[4]) / 100
#anchor_point = 0.8 #if 80% of your proteins have a conserved secondary structure dssp code in a particular column of your MSA, then the column will be considered as part of the anchor region

output_file = sys.argv[5]+'.fasta' #MSA with residues (not dssp code) but that has been squeezed towards anchor points


read_dsspMSA=AlignIO.read(dssp_file,'fasta')
alignment_length = read_dsspMSA.get_alignment_length()

dsspMSA_array = np.array([list(rec.seq) for rec in read_dsspMSA])

#filter non-dssp rows
sieve = [x for x in np.unique(dsspMSA_array) if x not in dssp_filter]
if sieve != []:
    filtered_dsspMSA_array_list=[]
    for item in dsspMSA_array:
        m = [x for x in item if x not in dssp_filter]
        if m ==[]:
            filtered_dsspMSA_array_list.append(item)
    print (dsspMSA_array)
    dsspMSA_array = np.stack( filtered_dsspMSA_array_list, axis=0 )


# Transpose the MSA so we can read per column instead of per protein (row)
transposed_dsspMSA_array = np.transpose(dsspMSA_array)

#find positions conserved 2nd structure elements
structure_elements = []
last_structured = 0
for conserved in conserved_2structure_dssp:
    positions = []
    for position,column in enumerate(transposed_dsspMSA_array):
        count = list(column).count(conserved)
        freq = count/len(list(column))
        if freq >= anchor_point:
            positions.append(position) #positions part of the anchor region
    if positions == []:
        print('The secondary structure element',conserved,'was not found in all sequences')
        continue
    else:
    #find start and end anchor points
        start = positions[0]
        stop = None
        counter = 1
        for i in range(1, len(positions)-1):
            if positions[i] - positions[i - 1] != 1:
                stop = positions[i-1]
            if stop != None:
                structure_elements.append((conserved+str(counter),(start,stop+1)))
                counter +=1
                stop = None
                start = positions[i]
        structure_elements.append((conserved+str(counter),(start,positions[-1])))
        if last_structured < positions[-1]:
            last_structured = positions[-1]

structure_elements = sorted(structure_elements, key=lambda x: x[1][0]) #sort it in ascending order

#add the unstructured regions
#here, with loop in mean all 2nd structure elements that are not in conserved_2structure_dssp 
loop_regions = []
loop_regions.append(('loop0', (0, structure_elements[0][1][0])))
# Iterate through the existing secondary structure elements and add the loops in between
for i in range(len(structure_elements) - 1):
    ss_start = structure_elements[i][1][1]
    ss_end = structure_elements[i + 1][1][0]
    loop_name = f'loop{i + 1}'
    loop_regions.append((loop_name, (ss_start, ss_end)))
loop_regions.append((f'loop_final', (last_structured, alignment_length)))

#all 2nd structure elements
all_regions = loop_regions
all_regions.extend(structure_elements)
all_regions = sorted(all_regions, key=lambda x: x[1][0])

# Get the minimum loop's length
loop_len = []
for region in loop_regions:
    gaps_list = []
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        gaps_list.append(sequence[region[1][0]:region[1][1]].count("-")) # count number of gaps in one loop
    loop_len.append(len(sequence[region[1][0]:region[1][1]]) - min(gaps_list)) #give the number of gaps in the sequence with the most AA in one loop
    loop_len.append(0)

# Squeeze the alignment
with open(output_file, "w+") as outfile:
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        outfile.write(">" + name + "\n")
        new_align = str()
        for region, length in zip(all_regions, loop_len):
            if region[0].startswith("loop"):
                nogaps_seq = sequence[region[1][0]:region[1][1]].replace("-", "")
                if region[0] == "loop0":
                    new_align += "-"*(length - len(nogaps_seq)) + nogaps_seq
                elif region[0] == "loop_final":
                    new_align += nogaps_seq + "-"*(length - len(nogaps_seq))
                else:
                    new_align += nogaps_seq[0:len(nogaps_seq)//2] + "-"*(length - len(nogaps_seq)) + nogaps_seq[len(nogaps_seq)//2:]
            else:
                new_align += sequence[region[1][0]:region[1][1]]
        outfile.write(new_align + "\n")