
process cdHitCollapse{
    tag "${sequencesValid.name}"
    publishDir "$params.outFolder/seqs/CD-HIT", mode: "copy"

    input:
    path sequencesValid
    val clustering

    output:
    path "*_collapsed.clstr", emit: clusters
    path "*_collapsed", emit: seqs
    env num , emit: num

    script:
    """
    if (($clustering >=70)); then
        n=5
    elif (($clustering >=60)); then
        n=4
    elif (($clustering >=50)); then
        n=3
    else
        n=2
    fi

    result_perc=\$(echo $clustering 100 | awk '{ print \$1/\$2 }')
    result_perc=\$(echo "\$result_perc" |tr ',' '.')
    cd-hit -i ${sequencesValid} -o ${sequencesValid.baseName}_${clustering}_collapsed  -c \$result_perc -n \$n -d 200 

    num=\$(grep "total seq:" .command.log | awk '{print \$NF}')
    """

}

process cdHitSubsetting{
    tag "${sequencesValid.name}"
    publishDir "$params.outFolder/seqs/CD-HIT_for_t-coffee", mode: "copy"
    
    errorStrategy 'retry'
    maxRetries 15

    input:
    path sequencesValid
    val initclustering
    val minSubsetSimilarity
    val maxSubsetSize

    output:
    path "*.clstr", emit: clusters
    path "*.fasta", emit: seq

    script:
    """
    result_perc=\$(echo $initclustering 5 | awk '{ print \$1+\$2-\$1*\$2 }')
    result_perc=\$(echo "\$result_perc" |tr ',' '.')

    clustering=\$(echo \$result_perc 100 | awk '{ print \$1/\$2 }')
    clustering=\$(echo "\$clustering" |tr ',' '.')
    
    outname=${sequencesValid.baseName}_\${clustering}_clustered

    if ((\$clustering >=40)); then
         if ((\$clustering >=70)); then
            n=5
        elif ((\$clustering >=60)); then
            n=4
        elif ((\$clustering >=50)); then
            n=3
        else
            n=2
        fi
        cd-hit -i ${sequencesValid} -o \$outname  -c \$result_perc -n \$n -d 200
    else
        $projectDir/bin/psi-cd-hit.pl -i ${sequencesValid} -o \$outname  -c \$result_perc 
    fi

    python3 $projectDir/bin/cluster_to_subset.py ${sequencesValid} \${outname}.clstr \$result_perc $minSubsetSimilarity $maxSubsetSize
    """
    //INFO: 
    //python script will fail if clusters are too big until minSubsetSimilarity (default 20%), then it will split randomly OR just make the largest possible subsets
}


process userSubsetting {
    publishDir "$params.outFolder/seqs/for_t-coffee/", mode: "copy"
    input:
    path sequenceFastas
    path foundModels

    output:
    path "*.fasta" , emit: subsetSeqsForTC

    script:

    """
    python3 $projectDir/bin/user_subsetting.py $sequenceFastas $foundModels
    """
}
