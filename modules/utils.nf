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

    if(( $all == 0 )) ; then
        echo "ERROR! No sequences found!"
        exit 1
    fi

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

    val structuresFolder //params.structures
    val foundModelsCount // joined.modelFound.count()
    val foundModelsAF //afjoined.modelFound.count()
    val foundModelsESM //esmfjoined.modelFound.count()
    val structurelessCount 
    path structurelessFile //structureless_seqs

    val createSubsets
    val useSubsets
    path subsetClstr
    path seqsToAlignFile //foundSeqs
    val tcparams
    path mergedMSA  //finalMsa
    val mafftparams


    output:
    path "*.md"

    script:

    """
    outname=simsapiper_summary
    outfile=\${outname}.md
    
    echo "# Summary file for SIMSApiper " >> \$outfile
    echo "Find detailed explations for each step on [GitHub](https://github.com/Bio2Byte/simsapiper) "  >> \$outfile
    echo "# Input sequences files " >> \$outfile
    echo ' No. of sequences in inputfile(s): ' $foundSequencesCount >> \$outfile
    echo '\n Input files can be found in $seqsdir' >> \$outfile

    echo '\n SIMSApiper found these files: '$inputSeqFiles >> \$outfile
    

    echo ' # 1 Data preparation' >> \$outfile
    echo " ## 1.2 Data reduction with CD-Hit" >> \$outfile
    echo ' + Similarity cutoff for CD-Hit: ' $dropSimilar >> \$outfile
    echo ' + No. of sequences after collapsing: ' $collapsedSequencesCount >> \$outfile
    echo ' + Included selected sequences to respresenting sequence file: ' $favoriteSeqs >> \$outfile
    echo ' + Find representing and represented sequences here: ' >> \$outfile
    echo \$(readlink -f $similarClstr ) >> \$outfile


    echo " ## 1.3 Invalid sequences " >>\$outfile
    echo '+ $seqsInvalidCount sequences removed because they contained more then $seqQC % non-standard or unresolved amino acids' >> \$outfile

    echo '+ Find invalid sequences here:' >> \$outfile
    echo \$(readlink -f $seqsInvalidFile) >> \$outfile

    echo '# 2 Structural information'  >>  \$outfile
    echo "## Identify missing models:" >>  \$outfile
    echo '+ Models found in '$structuresFolder' folder at start: '$foundModelsCount >>  \$outfile

    echo "## Retrieve structural information"  >>  \$outfile
    echo '+ Models found in AlphaFold Database: '$foundModelsAF >>  \$outfile
    echo '+ Models generated by ESMF: '$foundModelsESM >>  \$outfile
    echo '+ No. of sequences not matched to a model: ' $structurelessCount  >> \$outfile

    echo '++ Find sequences not matched to a model here:' >> \$outfile
    echo \$(readlink -f $structurelessFile) >> \$outfile
	echo '+ Models generated by local ESMF: '$structurelessCount >>  \$outfile
    
    echo '# 3 Generation of subsets' >> simsapiper_summary.md
    echo '+ Find which subset includes which sequence here: '  >> \$outfile
    echo \$(readlink -f $subsetClstr) >> simsapiper_summary.md
    ## echo '+ orphan seqs'

    echo '## 3.1 Subsets were generated from input files:' $useSubsets >> simsapiper_summary.md
    if [ "$createSubsets" != "false" ] ; then
        echo '## 3.2 Automatically generated subsets : true ' >> \$outfile
        echo '+ Subsets were generated' $createSubsets '% sequence similarity cutoff' >> simsapiper_summary.md
    else
        echo '## 3.2 Automatically generated subsets : false ' >> \$outfile
    fi
    
    echo '# 4 Align each subsets with T-Coffee' >> simsapiper_summary.md
    echo '+ Cleaned, matched sequences to submit to T-coffee: ' >> \$outfile
    echo \$(readlink -f $seqsToAlignFile) >> simsapiper_summary.md
    echo 'Additional T-Coffee parameters: ' $tcparams >> \$outfile

    echo '# 5 Align subset alignments with MAFFT' >> simsapiper_summary.md
    echo 'Alignment of subset alignments, structureless and orphan sequences '  >> \$outfile
    echo \$(readlink -f $mergedMSA) >> simsapiper_summary.md
    echo 'Additional MAFFT parameters: ' $mafftparams >> \$outfile
        



    echo "$outdir"  >> \$outfile
    

    
    """    

    //echo  \$(readlink -f $seqsInvalidFile) >> \$outfile
    //inputSeqFilePath=\$(readlink -f $inputSeqFiles )
    //echo '[$inputSeqFiles]('\$inputSeqFilePath')' >> \$outfile

    //[link](file:///Users/matb/Desktop/cat.gif) 
    //echo '<a href="file://'\$inputSeqFilesPath'">link</a>' >>\$outfile



    //inputSeqFilesPath=\$(readlink -f $inputSeqFiles) 
    //echo '\n SIMSApiper found these files: ![$inputSeqFiles](file://'\$inputSeqFilesPath')' >> \$outfile
    

}