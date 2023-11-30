import pandas as pd
import sys
import math


seqs_file = sys.argv[1]
clusters_file = sys.argv[2]
perc = float(sys.argv[3])

minSubset = sys.argv[4]
clusterTheClusters = False
if "min" in minSubset:
    clusterTheClusters = True
else:
    try:
        minSubsetSimilarity = int(sys.argv[4])/100
    except:
        print("--minSubsetId can only be 'min' or a percentage between 20-100, not " + minSubset)

maxSubsetSize = sys.argv[5]

orphanPerc = 0.05 #if too many orphans, no benefit of structure alignment

def cd_hit_to_df(clusters_file, seq_file):
    #read cdhit
    with open(clusters_file, 'r') as reps:
        match_dic = {}
        for line in reps:
            if line.startswith('>Cluster'):
                rep='Subset'+line.split(' ')[1].replace('\n','')
            else:
                id='>'+line.split('>')[1].split('..')[0]
                match_dic[id] = rep
    cluster_df=pd.DataFrame.from_dict(match_dic, orient='index')

    cluster_df = cluster_df[0].rename('subset')
    cluster_df= cluster_df.reset_index()
    
    #read fasta
    df = pd.read_csv(seq_file, header=None)
    seq_df = pd.DataFrame({'index':df[0].iloc[::2].values, 'seq':df[0].iloc[1::2].values})

    df = pd.merge(seq_df, cluster_df, on='index')
    return(df)

def cluster_fasta(rep_df, perc): 
    #general data info
    total_no_seqs = len(rep_df)
    all_subsets = rep_df.subset.unique().tolist()

    #select cutoff
    rep_df['seq_len'] = rep_df.seq.str.len()
    no_smallbois = len(rep_df[rep_df.seq_len < 400]) 
    print ('seqs longer then 400AA: ', total_no_seqs-no_smallbois)

    percc = no_smallbois/total_no_seqs
    print(total_no_seqs)
    if maxSubsetSize == 'true':
        if percc < 0.9:
            print (100-round(percc*100),'% of sequences are longer then 400 AA, therefore the subsets must be small!')
            cutoff = 50
        else:
            cutoff = 100
    else:
        cutoff = int(maxSubsetSize)
        
    #check if cluster size appropriate
    info = rep_df.groupby('subset').agg(counter=pd.NamedAgg(column='seq', aggfunc='count'), seqlen = pd.NamedAgg(column='seq_len', aggfunc='mean') )
    info = info.reset_index()

    info.sort_values("seqlen", ascending=[False],inplace=True)

    giantclusters = info[info.counter > cutoff] 
    no_giantclusters = len(giantclusters)

    singletons = info[info.counter < 2] 
    no_singletons = len(singletons)

    if clusterTheClusters:
        #make clusters as big as they can
        new_big_cluster_ids = []
        cluclusters={}
        members=0
        setnum=0
        all_min_subsets=['Subset1']

        for row in range(0,len(info)):
            members += info.counter.iloc[row]
            if members <= cutoff:
                cluclusters[info.subset.iloc[row]]= 'Subset'+str(setnum)
            else:
                setnum+=1
                all_min_subsets.append('Subset'+str(setnum))
                cluclusters[info.subset.iloc[row]]= 'Subset'+str(setnum)

                members=info.counter.iloc[row]
                if members > cutoff:
                    print (info.subset.iloc[row])
                    new_big_cluster_ids.append('Subset'+str(setnum))
    
        cluclusters[info.subset.iloc[row]]= 'Subset'+str(setnum)

        rep_df.rename({"subset":"cluster"},inplace=True, axis='columns')
        rep_df['subset']=rep_df['cluster'].map(cluclusters)

        print(new_big_cluster_ids)
        split_per_subset(all_min_subsets,rep_df,new_big_cluster_ids,cutoff)

    else:
        ##fail if clusters too big to recluster at lower seq id; unless already at 20%; then random split
        toobig_cluster_ids = []
        if no_giantclusters != 0:
            if perc > minSubsetSimilarity:        
                raise ValueError('Cluster with too many members detected:', no_giantclusters)
            else:
                toobig_cluster_ids = giantclusters.subset.values.tolist()
                print ('clusters with more than ', cutoff,' members: ', no_giantclusters)
        
        if no_singletons/len(all_subsets) > orphanPerc:
            if perc > minSubsetSimilarity: 
                raise ValueError('Too many singleton clusters detected:', no_singletons )
     
        split_per_subset(all_subsets,rep_df,toobig_cluster_ids,cutoff)
            
            

def split_per_subset(all_subsets,rep_df,toobig_cluster_ids,cutoff):
       
    print ('all subsets:',all_subsets)
    for id in all_subsets:

        subset_df = rep_df[rep_df.subset == id]
        if subset_df.empty:
            print ('Skipping because first cluster too big')
            continue
        #create subsubsets if cluster still too large
        size = len(subset_df)
        print(subset_df)
        factor = math.ceil(size/cutoff)
        print(factor)
        if factor <1:
            factor=1
        chunk = math.ceil(size/factor)
        i =1
        #print (chunk)

        if id in toobig_cluster_ids:
            print(subset_df)
            print(size)
            for start in range(0, size, chunk):
                print("subset:", i, 'from ',start, start + chunk )
                sub_subset_df = subset_df.iloc[start:start + chunk]
                sub_id = id +'_' +str(i)
                write_fasta_from_df(sub_subset_df,sub_id)
                i +=1
        elif len(subset_df) < 2:
            write_fasta_from_df(subset_df,'orphan_'+id)
        else:
            write_fasta_from_df(subset_df,id)   

def write_fasta_from_df(df,outname):
    out_name = outname + '.fasta'

    out_df = df[['index', 'seq']]
    out_df.iloc[:,0]= out_df.iloc[:,0] + '\n'

    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')

    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ','')
            m.write( row + '\n')



df = cd_hit_to_df(clusters_file, seqs_file)
print (df)
cluster_fasta(df, perc)

