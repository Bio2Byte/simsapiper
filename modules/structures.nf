process fetchEsmAtlasStructure {
    tag "${header}"
    publishDir "$params.structures", mode: "copy"
    
    errorStrategy { task.attempt > 3 ? 'retry' : 'ignore' }

    input:
    tuple val(header), val(sequence)
    val noOfSeqs

    output:
    path "${header}.pdb", emit: esmStructures

    script:
    """
    echo Folding sequence ${header} using ESM Atlas

    sleep \$((RANDOM % $noOfSeqs))
    curl -X POST -k -connect-timeout 8 --data "${sequence}" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${header}.pdb

    find . -type f -name "*.pdb" -size -2k -exec bash -c 'mv "\$1" "\${1%.pdb}".fail' - '{}' +
    """
    //INFO:
    //sleep is needed to prevent failure to to api rate limit
}


process getAFmodels {
    publishDir "$params.structures", mode: "copy"
    errorStrategy 'ignore'
    tag "${id}"

    input:
    each id 

    output:
    path "${id}.pdb"

    script:
    """
    echo Fetching model ${id} from AFDB

    result=\$(echo ${id} | cut -d '_' -f 1)
    echo \$result

    curl https://alphafold.ebi.ac.uk/files/AF-\${result}-F1-model_v4.pdb -o ${id}.pdb

    find . -type f -name "*.pdb" -size -2k -exec bash -c 'mv "\$1" "\${1%.pdb}".fail' - '{}' +
    """
}

process runDssp{
    
    publishDir "$params.dsspPath", mode: "copy"

    input:
    path model
    val gate

    output:
    path "*.dssp" , emit:dsspout

    script:
    """
    echo Gate is open $gate
    mkdssp -i $model -o ${model.baseName}.dssp   
    """
}


process esmFolds{
    publishDir "$params.structures", mode: "copy"
    
    input:
    path structureless
    val apptainerPath

    output:
    path "*.pdb", emit: esmFoldsStructures 
    val true, emit: gate
    path "esm_fold_statistics.csv"  , emit: structurestats
    path "*.pae"
    

    script:
    """
    export TORCH_HOME=$apptainerPath
    echo $structureless
    
    if [ -z "$structureless" ] ; then
        echo "There is no proteins to fold here"
    else
        head $structureless
        python3 $projectDir/bin/esmfold_inference.py --chunk-size 32  --max-tokens-per-batch 0 -i $structureless -o .
        python3 $projectDir/bin/extract_plddt.py . esm_fold_statistics.csv .
    fi
    """
    //export TORCH_HOME=\$VSC_SCRATCH_VO/ESMFold
    //--max-tokens-per-batch 1024 is standard --max-tokens-per-batch 0 to turn off  short-seq-grouping --max-tokens-per-batch 512
}