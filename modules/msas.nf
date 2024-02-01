
process runTcoffee {
    publishDir "$params.outFolder/msas/t-coffee", mode: "copy" 
    tag "$seqsToAlign"

    errorStrategy {task.attempt < 3 ? 'retry' : 'finish' }

    input:
    path seqsToAlign
    path strucsToAlign
    val tcoffeeParams
    val gate

    output:
    path "*.aln", emit : msa

    script:
    """
    head $seqsToAlign

    echo $tcoffeeParams
    echo '\n ---------------------- \n'

    set -e # Exit immediately if any command fails

    # Set up environment variables
    export PATH="/software:/software/.t_coffee/bin/linux:/usr/local/bin:\$PATH"
    export PLUGINS_4_TCOFFEE="/software/.t_coffee/plugins/linux"
    export MAX_N_PID_4_TCOFFEE="\${MAX_N_PID_4_TCOFFEE:-5000000}"

    # Default values
    WORKING_DIRECTORY_4_TCOFFEE=\$PWD
    SEQUENCE=$seqsToAlign
    PDB_DIR=$strucsToAlign
    CACHE_DIRECTORY=""
    TMP_DIRECTORY=""
    OUTPUT_DIRECTORY=\$PWD

    # Check if all mandatory parameters are provided
    if [[ -z "\$WORKING_DIRECTORY_4_TCOFFEE" || -z "\$SEQUENCE" || -z "\$PDB_DIR" || -z "\$OUTPUT_DIRECTORY" ]]; then
    echo "Error: Missing arguments."
    display_usage
    exit 1
    fi

    # Set default cache and tmp directories if not provided
    if [ -z "\$CACHE_DIRECTORY" ]; then
    CACHE_DIRECTORY="\$WORKING_DIRECTORY_4_TCOFFEE/cache"
    fi

    if [ -z "\$TMP_DIRECTORY" ]; then
    TMP_DIRECTORY="\$WORKING_DIRECTORY_4_TCOFFEE/tmp"
    fi

    # Create necessary directories
    mkdir -p "\$CACHE_DIRECTORY"
    mkdir -p "\$TMP_DIRECTORY"
    OUTPUT_PATH="\$OUTPUT_DIRECTORY"
    mkdir -p "\$OUTPUT_PATH"

    export PDB_DIR=\$PDB_DIR

    echo "Working directory: \$WORKING_DIRECTORY_4_TCOFFEE"
    echo "Cache directory: \$CACHE_DIRECTORY"
    echo "Temp directory: \$TMP_DIRECTORY"
    echo "Output directory: \$OUTPUT_PATH"
    echo "Input sequence: \$SEQUENCE"
    echo "PDB directory: \$PDB_DIR"
    echo "Environment variables:"
    echo "PATH: \$PATH"
    echo "PLUGINS_4_TCOFFEE: \$PLUGINS_4_TCOFFEE"
    echo "MAX_N_PID_4_TCOFFEE: \$MAX_N_PID_4_TCOFFEE"

    # Run T-Coffee and save output to a unique log file
    TIMESTAMP=\$(date +"%Y%m%d%H%M%S")
    LOG_FILE="\$WORKING_DIRECTORY_4_TCOFFEE/tcoffee_\${TIMESTAMP}.log"

    echo "Running T-Coffee..."
    t_coffee ${params.tcoffeeParams ? tcoffeeParams : ''} -thread 0 -in="\$SEQUENCE" \
        -method TMalign_pair \
        -evaluate_mode=t_coffee_slow \
        -mode=3dcoffee \
        -pdb_min_cov=1 \
        -outfile="\$OUTPUT_PATH/aligned_${seqsToAlign.baseName}.aln" > "\$LOG_FILE" 2>&1

    echo "T-Coffee execution completed successfully."

    if grep -q "proba_pair" "\$LOG_FILE"; then
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