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

process summary {
    publishDir "$params.outFolder", mode: "copy"

    input:
    path inputSeqFiles
    val foundSequencesCount

    val dropSimilar
    val collapsedSequencesCount
    path similarClstr 

    val seqsInvalidCount
    val structurelessCount
    path finalMsa
    

    output:
    path '*.md'

    script:
    """
    touch simsapiper_summary.md

    echo "# Input sequences files" >> simsapiper_summary.md
    echo 'No. of sequences in inputfile(s): ' $foundSequencesCount >> simsapiper_summary.md
    echo $inputSeqFiles >> simsapiper_summary.md

    echo "# Data reduction with CD-Hit" >> simsapiper_summary.md
    echo 'Similarity cutoff for CD-Hit: ' $dropSimilar
    echo 'No. of sequences after collapsing: ' $collapsedSequencesCount >> simsapiper_summary.md
    echo [$similarClstr.baseName]($similarClstr) >> simsapiper_summary.md



    """
    //check if all sequences been aligned 
    //input = collapsed(input-clusters,sum all input files) (paths:90clstr,who represents who) + invalid (path)+ structureless (path)+ matched to model(path seqs to align)
    //final_msa(path) = matched_to_model + structureless 
    //squeeze impact = gap counter: pre, post squeeze
    //seq id
    //[link](file:///d:/absolute.md)
}

process missingQC {
    errorStrategy 'finish'

    input:
    val all 
    val missing
    val strucQC

    output:
    val true, emit: gate

    script:

    """

    result=\$(echo $missing $all | awk '{ print \$1/\$2 }')
    result=\$(echo "\$result" |tr ',' '.')

    result_perc=\$(awk -v res=\$result -v perc=100 'BEGIN{print (res*perc); }')
    result_perc=\$(echo "\$result_perc" |tr ',' '.')

    message="\$result_perc % of sequences not could be matched to a model, please check the structureless_sequences.fasta in the output directory! Try using Colabfold (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) to provide the missing models!"

    if (( \$(echo "\$result_perc $strucQC" | awk '{ print (\$1 >\$2)}') )); then
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
