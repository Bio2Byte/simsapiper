process runTcoffee {
    publishDir "$params.outFolder/msas/t-coffee", mode: "copy" 
    tag "$seqsToAlign"

    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate' }

    input:
    path seqsToAlign
    path strucsToAlign
    val tcoffeeParams
    val gate

    output:
    path "*.aln", emit : msa

    script:
    """
    # Set up environment variables
    export PATH="/software:/software/.t_coffee/bin/linux:/usr/local/bin:\$PATH"
    export PLUGINS_4_TCOFFEE="/software/.t_coffee/plugins/linux"
    export MAX_N_PID_4_TCOFFEE="\${MAX_N_PID_4_TCOFFEE:-5000000}"

    # Default values
    WORKING_DIRECTORY_4_TCOFFEE=\$PWD
    SEQUENCE=$seqsToAlign

    CACHE_DIRECTORY="\$WORKING_DIRECTORY_4_TCOFFEE/cache"
    mkdir -p "\$CACHE_DIRECTORY"
    echo "Cache directory: \$CACHE_DIRECTORY"

    TMP_DIRECTORY="\$WORKING_DIRECTORY_4_TCOFFEE/tmp"
    mkdir -p "\$TMP_DIRECTORY"
    echo "Temp directory: \$TMP_DIRECTORY"

    PDB_DIR=$strucsToAlign
    echo "PDB directory: \$PDB_DIR"
    export PDB_DIR=\$PDB_DIR

    echo "Running T-Coffee..."
    t_coffee ${params.tcoffeeParams ? tcoffeeParams : ''} -thread ${task.cpus}  -in="\$SEQUENCE" \
        -method TMalign_pair \
        -mode=3dcoffee \
        -outfile=${seqsToAlign.baseName}.aln  -debug  

    echo "T-Coffee execution completed successfully."

    if grep -q "proba_pair" .command.log ; then
        echo "Error: File contains proba_pair, a HMM-based alignment tool that T-Coffee automatically runs when something is wrong with the structure models. If you want to proceed anyways, remove these lines from modules/msas.nf"
        exit 1
        fi
    """
}

process mergeMafft {
    publishDir "$params.outFolder/msas", mode: "copy"  

    input:
    path seqs
    val mafftParams
    val outName

    output:
    path "merged_${outName}_alignment.fasta" , emit: finalMsa
    
    script:
    """
    numSeqFiles=\$(echo "$seqs" | wc -w)
    echo \$numSeqFiles
    echo $mafftParams

    if (( \$numSeqFiles > 1 ))
    then
        python3 $projectDir/bin/MAFFT_merge.py '$seqs'
        mafft ${mafftParams ? mafftParams : ''} --merge tableMAFFT.txt input > prefinalMSA.fasta
    else
        echo "Only one alignment has been generated: $seqs"  
        cp $seqs prefinalMSA.fasta
    fi

    python3 $projectDir/bin/convert_to_fasta.py prefinalMSA.fasta fasta merged_${outName}_alignment
    """
    //INFO: mafftParams = "" This is a recommended option, MAFFT merge changes everything otherwise
    //--inputorder enables to keep sequence order of inputted files (recommended)
    //--op gap penalty (default 1.53)

    //different modes (only iterative allignment methods (iterative + progressive), possible to use progressive alignment only see MAFFT):
    //see their documentation https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html#GLE
    //FFT-NS-i
    //mafft --merge $table --maxiterate 1000 $input > $output

    //E-INS-i =>  recommended for <200 sequences with multiple conserved domains and long gaps 
    //mafft --merge $table --genafpair --maxiterate 1000 $input > $output 

    //L-INS-i (Very slow; recommended for <200 sequences with one conserved domain and long gaps)
    //mafft --merge $table --localpair --maxiterate 1000 $input > $output

    //G-INS-i (Very slow; recommended for <200 sequences with global homology)
    //mafft --merge $table --globalpair --maxiterate 1000 $input > $output 

    //FFT-NS-2 recommended for very big datasets >2,000 sequences
    //mafft --merge $table --retree 2 $input > $output

}

process mapDssp{
    publishDir "$params.outFolder/msas", mode: "copy"

    input:
    path dssps
    path msa

    output:
    path "dssp_${msa.baseName}.fasta" , emit: mmsa
    path "*.fasta"

    script:
    """
    python3 $projectDir/bin/map_dssp.py $msa dssp_${msa.baseName} dssp
    """
}


process squeeze{
    publishDir "$params.outFolder/msas", mode: "copy"

    input:
    path msa
    path dssp
    val squeeze
    val squeezePerc

    output:
    path "*.fasta" , emit: msa

    script:
    """
    python3 $projectDir/bin/squeeze_msa.py $msa $dssp "$squeeze" $squeezePerc squeezed_${msa.baseName} 
    """
    //INFO: 
    //choose conserved secondary structure elements according to dssp across your dataset as it is does elements that TCOFFEE will use to align your proteins
}

process reorder{
    publishDir "$params.outFolder", mode: "copy"

    input:
    path finalMsa
    path inputSeqs
    val reorder

    output:
    path "reordered_${finalMsa.baseName}.fasta", emit: msaOrga
    
    script:
    """
    for file in *; do
    if [ -f "\$file" ]; then
        mv "\$file" "\$file.tmp"
        sed 's/,/-/g' "\$file.tmp" > "\$file"
        rm "\$file.tmp"
    fi
    done

    python3 $projectDir/bin/reorganize_output.py "$inputSeqs" "$reorder" $finalMsa reordered_${finalMsa.baseName} 
    """

}