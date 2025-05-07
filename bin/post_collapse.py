import pandas as pd 
import sys

inputFile=sys.argv[1]
selectedFile=sys.argv[2]
favs=sys.argv[3]
collapsedFileName=sys.argv[4]


def write_fasta_from_df(df,labelcol,seqcol,outname):
    out_name = outname + '.fasta'
    out_df = df.filter(items=[labelcol, seqcol])
    out_df[labelcol] =  out_df[labelcol]+"\n"
    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')
    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ',"")
            m.write( '>'+row + '\n')

def read_fasta_to_df(seq_file):
    df = pd.read_csv(seq_file, header=None)#, sep=' ') #sometimes helps 
    seq_df = pd.DataFrame({'label':df[0].iloc[::2].values, 'seq':df[0].iloc[1::2].values})

    seq_df.label= seq_df.label.str.lstrip(' ')
    seq_df.label= seq_df.label.str.lstrip('>')

    return seq_df


findFavs=False
#check if (optional) favs are present
if favs!='false':
    favsl=favs.split(',')
    findFavs=True

inputDf=read_fasta_to_df(inputFile)
selectedDf=read_fasta_to_df(selectedFile)

##if no, add them
missingfavs=[]
if findFavs:
    for elem in favsl:
        if elem in selectedDf.label.tolist():
            print ("found ", elem)
        else:
            favseq=inputDf[inputDf['label'] ==elem ]
            missingfavs.append([elem,favseq.iloc[0]['seq']])
if missingfavs !=[]:
    favselectedDf=pd.concat([selectedDf,pd.DataFrame(missingfavs, columns=selectedDf.columns)])
else :
    favselectedDf=selectedDf

write_fasta_from_df(favselectedDf,'label','seq',selectedFile)


#find the sequences that were not selected / are too similar and collect in output file 
inputLabels=inputDf.label.tolist()
selectedLabels=favselectedDf.label.tolist()



collapsedDf=inputDf[~inputDf.label.isin(favselectedDf.label)].dropna()

print(inputDf)
print (inputDf[~inputDf.label.isin(favselectedDf.label)])
print(collapsedDf)



if len(inputDf)-len(favselectedDf) != len(collapsedDf):
    print ("something went wrong here")
else:
    if len(collapsedDf)>0:
        write_fasta_from_df(collapsedDf,'label','seq',collapsedFileName)
    else:
        print("No sequence was excluded here!")