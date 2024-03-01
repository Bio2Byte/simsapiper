process writeFastaFromChannel{
    publishDir "$params.outFolder/seqs", mode: "copy"
    tag "${outname}"

    input:
    val seqlist
    val outname

    output:
    path "*.fasta", emit: found

    script:

    """
    cleanseq=\$( echo "$seqlist" | cut -c2- | rev | cut -c2- |rev | tr -d ' ' )
    echo "\$cleanseq" | tr ',' '\n'  > $outname
    """
}



process attendance{
    publishDir "$params.outFolder", mode: "copy"

    input:
    val collapsedSequencesCount
    val seqsInvalidCount
    val structurelessCount
    path finalMsa
    val foundSequencesCount

    output:
    path '*.txt'

    script:
    """
    echo $finalMsa
    fin=\$(grep -c ">" $finalMsa)

    echo 'No. of sequences in inputfile(s): ' $foundSequencesCount >> "sequence_report.txt"

    echo 'No. of invalid sequences: ' $seqsInvalidCount >> "sequence_report.txt"
    echo 'No. of sequences not matched to a model: ' $structurelessCount >> "sequence_report.txt"
    echo 'No. of sequences matched to a model: ' $collapsedSequencesCount >> "sequence_report.txt"

    echo 'No. of sequences in final alignment: ' \$fin >> "sequence_report.txt"

    av_conservation=`python3 $projectDir/bin/shannons_entropy.py $finalMsa`
    echo 'Average sequence conservation (Shannons Entropy): ' \$av_conservation 


    if (( \$fin !=$collapsedSequencesCount + $structurelessCount ))
    then
        echo "ERROR: Not all valid sequences are found in the output file, please check $finalMsa in in the output directory!"
        exit 1
    fi

    """
}

process missingQC {
    errorStrategy 'finish'

    input:
    val all 
    val found
    val strucQC

    output:
    val true, emit: gate

    script:

   """
    result=\$(echo $found $all | awk '{ print \$1/\$2 }')
    result=\$(echo "\$result" |tr ',' '.')

    result_perc=\$(awk -v res=\$result -v perc=100 'BEGIN{print (res*perc); }')
    result_perc=\$(echo "\$result_perc" |tr ',' '.')
    echo \$result_perc

    qc=\$(echo 100 $strucQC | awk '{ print \$1-\$2 }')
    echo \$qc

    message="\$result_perc % of sequences not could be matched to a model, please check the structureless_sequences.fasta in the output directory! Try using Colabfold (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) to provide the missing models!"

    if (( \$result_perc < \$qc )) ; then
        echo "ERROR!"
        echo \$message
        exit 1
    else 
        echo  \$message > structurecheck.txt
    fi

    """
}


process readSeqs {
    publishDir "$params.outFolder/seqs", mode: "copy"
    input:
    val format
    each sequenceFiles

    output:
    path "*.fasta", emit: convertedSeqs


    script:
    """ 
    echo  'converting $sequenceFiles with format $format into ${sequenceFiles.baseName}' 
    python3 $projectDir/bin/convert_to_fasta.py $sequenceFiles $format converted_${sequenceFiles.baseName}
    """
}


process createSummary{
    publishDir "$params.outFolder", mode:"copy"

    input:
    path outdir 

    path inputSeqFiles //allSequences
    val foundSequencesCount //fullInputSeqsNum
    val seqsdir // params.seqs

    val dropSimilar //param
    val collapsedSequencesCount //allSequencesCount
    val favoriteSeqs //param
    path similarClstr //cdHitCollapse.out.clusters

    val seqQC
    val seqsInvalidCount
    path seqsInvalidFile



    output:
    path "*.md"

    script:

    """
    outfile=simsapiper_summary.md
    
    echo "# Input sequences files " >> \$outfile
    echo ' No. of sequences in inputfile(s): ' $foundSequencesCount >> \$outfile
    echo '\n Input files can be found in' $seqsdir >> \$outfile
    echo '\n SIMSApiper found these files:' $inputSeqFiles >> \$outfile


    echo ' # 1 Data preparation' >> \$outfile
    echo " ## 1.2 Data reduction with CD-Hit" >> \$outfile
    echo ' + Similarity cutoff for CD-Hit: ' $dropSimilar >> \$outfile
    echo ' + No. of sequences after collapsing: ' $collapsedSequencesCount >> \$outfile
	echo ' + Included selected sequences to respresenting sequence file: ' $favoriteSeqs >> \$outfile
    echo ' + Find representing and represented sequences here: ' >> \$outfile
    echo \$(readlink -f $similarClstr ) >> \$outfile


    echo " ## 1.3 Invalid sequences " >>\$outfile
    echo '+ $seqsInvalidCount sequences removed because they contain more then $seqQC % non-standard or unresolved amino acids' >> \$outfile
    echo '+ Find invalid sequences here: ' >> \$outfile


    echo "$outdir"  >> \$outfile

    """    

    //echo  \$(readlink -f $seqsInvalidFile) >> \$outfile
    //inputSeqFilePath=\$(readlink -f $inputSeqFiles )
    //echo '[$inputSeqFiles]('\$inputSeqFilePath')' >> \$outfile

}