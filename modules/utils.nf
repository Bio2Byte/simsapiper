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
    val gate
    val outdir 

    path inputSeqFiles //allSequences
    val foundSequencesCount //fullInputSeqsNum
    val seqsdir // params.seqs

    val dropSimilar //param
    val collapsedSequencesCount //allSequencesCount
    val favoriteSeqs //param

    val seqQC
    val seqsInvalidCount

    val structuresFolder //params.structures
    val foundModelsCount // joined.modelFound.count()
    val retrieve
    val foundModelsAF //afjoined.modelFound.count()
    val model
    val foundModelsESM //esmfjoined.modelFound.count()
    val structurelessCount 

    val createSubsets
    val useSubsets
    val tcparams
    val mafftparams

    val dssp
    val dsspfiles //params.dsspPath
    val squeeze
    val squeezePerc

    val reorder
    val convert  // convertMSA

    output:
    path "*.md"

    script:

    """
    outfile=simsapiper_summary.md

    
    echo "# Summary file for SIMSApiper "
    echo "Find detailed explations for each step on [GitHub](https://github.com/Bio2Byte/simsapiper) "

    echo '\n All output files can be found here: ' $outdir
    
    echo "# Input sequences files "
    echo '* No. of sequences in inputfile(s): ' $foundSequencesCount
    echo '* Input files can be found in $seqsdir'

    echo '* SIMSApiper found these files: '$inputSeqFiles
    

    echo '# 1 Data preparation'
    echo "## 1.2 Data reduction with CD-Hit"
    echo '* Similarity cutoff for CD-Hit: ' $dropSimilar
    
    if [ "$dropSimilar" != "false" ] ; then
        echo '* Included selected sequences to respresenting sequence file: ' $favoriteSeqs
        echo '* No. of sequences after collapsing: ' $collapsedSequencesCount
        echo '* Find collapsed sequences here: '
        collapsed=$outdir/seqs/CD-HIT/removed*.fasta
        for i in \$collapsed; do
            echo '    * '\$i
        done

        echo '* Find representing and represented sequences here: '
        similarClstr=$outdir/seqs/CD-HIT/*.clstr
        for i in \$similarClstr; do
            echo '    * '\$i
        done
    fi

    echo "## 1.3 Invalid sequences "
    echo '* $seqsInvalidCount sequences removed because they contained more then $seqQC % non-standard or unresolved amino acids'
    if (( $seqsInvalidCount != 0 ));then
        echo '* Find invalid sequences here:' $outdir/seqs/too_many_unknown_characters.fasta
    fi

    echo '# 2 Structural information'
    echo "## 2.1 Identify missing models"
    echo '* Models found in '$structuresFolder
    echo '* Models found at start of the pipeline: ' $foundModelsCount

    echo '## 2.2 Search AFDB: ' $retrieve
    if [ "$retrieve" != "false" ] ; then
        foundAF=\$(( $foundModelsAF - $foundModelsCount ))
        echo '* Models found in AlphaFold Database: ' \$foundAF
    else
        foundAF=0
    fi
    
    echo '## 2.4 Use ESM Atlas: ' $model
    if [ "$model" != "false" ] ; then
        echo '* Models generated by ESMF: ' \$(( $foundModelsESM - \$foundAF - $foundModelsCount))
        fi
    echo '## 2.6 No. of sequences not matched to a model: ' $structurelessCount

    if  (( $structurelessCount != 0 )) ; then
        echo 'Find sequences not matched to a model here:' $outdir/seqs/structureless_seqs.fasta
        fi
    
    echo '# 3 Generation of subsets' 
    echo '## 3.1 Subsets were generated from input files:' $useSubsets 

    if [ "$createSubsets" != "false" ] ; then
        echo '## 3.2 Automatically generated subsets : true '
        echo '* Subsets were generated' $createSubsets '% sequence similarity cutoff' 
        echo '* Find which subset includes which sequence here: ' $outdir/seqs/CD-HIT_for_t-coffee/*.clstr

        echo '# 4 Align each subsets with T-Coffee'
        echo '* Cleaned, matched sequences to submit to T-coffee: '
        for i in $outdir/seqs/CD-HIT_for_t-coffee/*.fasta ;do
            echo '    * ' \$i
            done
    else
        echo '## 3.2 Automatically generated subsets : false '

        echo '# 4 Align each subsets with T-Coffee'
        echo '* Cleaned, matched sequences to submit to T-coffee: '
        for i in $outdir/seqs/for_t-coffee/*.fasta ;do
            echo '    * ' \$i
            done
    fi
    
    echo '* Additional T-Coffee parameters: ' $tcparams

    echo '# 5 Align subset alignments with MAFFT'
    echo '* Alignment of subset alignments, structureless and orphan sequences ' $outdir/msas/merged*.fasta
    echo '* Additional MAFFT parameters: ' $mafftparams
        
    echo '# 6 Run DSSP:' $dssp
    if [ "$dssp" != "false" ] ; then
        echo '* DSSP file can be found in' $dsspfiles

        echo '# 7 Improve MSA '
        echo '## 7.1 Map DSSP to MSA' 
        echo '* DSSP codes mapped to merged alignment:' $outdir/msas/dssp_merged*.fasta
    fi

    if [ "$squeeze" != "false" ] ; then
        echo '## 7.2 Squeeze MSA towards conserved secondary structure elements' 
        echo '* Categories selected for squeezing: ' $squeeze
        echo '* Threshold for region to be considered conserved: ' $squeezePerc
        echo '* Squeezed alignemnt: '$outdir/msas/squeezed_merged*.fasta

        echo '## 7.3 Map DSSP to squeezed MSA'
        echo '* DSSP codes mapped to merged alignment:' $outdir/msas/dssp_squeezed_merged*.fasta

        echo "## Gap reduction"
        count_file1=\$(grep -o '-' $outdir/msas/merged*.fasta | wc -l)
        count_file2=\$(grep -o '-' $outdir/msas/squeezed_merged*.fasta | wc -l)
        difference=\$((count_file1 - count_file2))
        echo '* The inital merged alignment has ' \$count_file1 ' gaps.'
        echo '* The the squeezed alignment has ' \$count_file2 ' gaps.'
        echo '* The squeezing step removed ' \$difference ' gaps from the alignment.'
    fi

    echo '# 8 Reorder MSA ' 
    if [ "$reorder" != "false" ] ; then
        echo '* Alignment reorded by: ' $reorder
        echo '* Reordered alignment: ' $outdir/reordered*.fasta
    fi

    
    if [ "$convert" != "false" ] ; then
        echo '# 9 Convert MSA'
        echo '* Converted MSA to ' $convert 'format' 
        echo '* Converted MSA: ' $outdir/msas/converted*.fasta
    else
        echo '# 9 Convert MSA: false'
    fi


    cp .command.out \$outfile
    
    """    
    
    // add new dependency
    //python3  $projectDir/bin/sequence_sim.py $finalMSA


    //av_conservation=`python3 $projectDir/bin/shannons_entropy.py $finalMsa`
   // echo 'Average sequence conservation (Shannons Entropy): ' \$av_conservation 

    //md improvemends
    //echo  \$(readlink -f $seqsInvalidFile)
    //inputSeqFilePath=\$(readlink -f $inputSeqFiles )
    //echo '[$inputSeqFiles]('\$inputSeqFilePath')'

    //[link](file:///Users/matb/Desktop/cat.gif) 
    //echo '<a href="file://'\$inputSeqFilesPath'">link</a>'

    //inputSeqFilesPath=\$(readlink -f $inputSeqFiles) 
    //echo '\n SIMSApiper found these files: ![$inputSeqFiles](file://'\$inputSeqFilesPath')'
    

}