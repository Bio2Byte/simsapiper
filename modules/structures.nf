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

    curl https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v4.pdb -o ${id}.pdb

    find . -type f -name "*.pdb" -size -2k -exec bash -c 'mv "\$1" "\${1%.pdb}".fail' - '{}' +
    """
}

process runDssp{
    publishDir "$params.outFolder/dssp", mode: "copy"

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
    //INFO: secondary structure elements according to dssp  
    //H = alpha-helix
    //B = beta-bridge residue
    //E = extended strand (in beta ladder)
    //G = 3/10-helix
    //I = 5-helix
    //T = H-bonded turn
    //S = beta-bend or beta-turn
}