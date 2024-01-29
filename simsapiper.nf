params.targetSequences  = "$params.seqs/*.$params.seqFormat"
targetSequencesFile     = file(params.targetSequences)
allSequences            = Channel.fromPath(params.targetSequences)
userStructures          = Channel.fromPath("$params.structures/*.pdb")


log.info """\

Usage:

\$ nextflow run simsapiper.nf -profile server,withsingularity --magic 

================================================================================
                                LIST OF PARAMETERS
================================================================================
                                GENERAL

Launch dir      : $launchDir
Project dir     : $projectDir
Execution time  : $params.executionTimestamp
================================================================================
                                INPUT FILES

Input file folder (--data): $params.data
Input sequence format (--seqFormat fasta): $params.seqFormat
    !Find all possible formats at https://biopython.org/wiki/SeqIO
Ignore sequences with % unknown characters (--seqQC 5): $params.seqQC    
Collapse sequences with % sequence identity (--dropSimilar 90): $params.dropSimilar
================================================================================
                                OUTPUT FILES

Output folder: $params.outFolder
Final MSA file name (--outName mymsa): $params.outName
================================================================================
                                SUBSETTING

Creates subsets of % sequence identity (--createSubsets 30): $params.createSubsets 
Set minimal % sequence identity in subsets (--minSubsetID 20): $params.minSubsetID
    !Reduce the number of subsets by collating with --minSubsetID "min"
Set maximal size of subset(--maxSubsetSize 100): $params.maxSubsetSize

User provides fasta files with subsets (--useSubsets): $params.useSubsets
    !Provide sequences not fitting any subset in a file containing 'orphan' in filename
================================================================================
                                ADD STRUCTURAL INFORMATION

Retrieve protein structure models from AFDB (--retrieve): $params.retrieve
Predict protein structure models with ESM Atlas (--model): $params.model
Predict protein structure models with Local ESMfold (--localModel): $params.localModel
    !Add factor for each multiple of 300 you expect to model to allocate more time
Maximal % of sequences not matched to a model (--strucQC 5): $params.strucQC
================================================================================
                                ALINGMENT PARAMETERS

Additional parameters for Tcoffee (--tcoffeeParams): $params.tcoffeeParams
Additional parameters for MAFFT (--mafftParams "--localpair --maxiterate 1000"): $params.mafftParams
Map DSSP to alignment (--dssp): $params.dssp
Squeeze alignment towards anchors (--squeeze "H,E"): $params.squeeze
    !Select secondary structure elements for anchor points according to DSSP annotation
Set minimal occurence % of anchor element in MSA: (--squeezePerc 80): $params.squeezePerc
Order MSA by input files (--reorder "cluster3.fasta,cluster4.fasta"): $params.reorder
Convert final MSA output from fasta (--convertMSA clustal): $params.convertMSA
================================================================================
"""

include {readSeqs as convertSeqs;
            readSeqs as convertTCMsa;
            readSeqs as convertFinalMsa;
            missingQC;
            attendance;
            writeFastaFromChannel as writeFastaFromMissing ;
            writeFastaFromChannel as writeFastaFromFound ;
            writeFastaFromChannel as writeFastaFromSeqsInvalid ;
            writeFastaFromChannel as writeFastaFromSeqsValid ;
} from "$projectDir/modules/utils"

include {userSubsetting;
            cdHitSubsetting;
            cdHitCollapse;
} from "$projectDir/modules/subsetting"

include {
    fetchEsmAtlasStructure;
    getAFmodels;
    runDssp;
    esmFolds;
} from "$projectDir/modules/structures"

include{
    runTcoffee ;
    mergeMafft ;
    reorder;
    mapDssp as mapDsspRough;
    mapDssp as mapDsspSqueeze;
    squeeze;
}from "$projectDir/modules/msas"


workflow {
    //convert input files to fasta format

    allSequences.view{'Sequence file found:' + it}
    convertSeqs(params.seqFormat, allSequences)
    sequenceFastas =  convertSeqs.out.convertedSeqs

    //data reduction with cdhit
    if (params.dropSimilar){
        cdHitCollapse(sequenceFastas, params.dropSimilar)
        reducedSeqs = cdHitCollapse.out.seqs
    } else {
        reducedSeqs = sequenceFastas
    }
    
    //Quality control input sequences 
    seqsRelabeled = reducedSeqs
        .splitFasta( record: [header: true,  sequence: true ])
        .map { record -> [header: record.header.replaceAll("[^a-zA-Z0-9]", "_"),
                sequence: record.sequence.replaceAll("\n","").replaceAll("[^ARNDCQEGHILKMFPOSUTWYVarndcqeghilkmfposutwvy]", "X")] }

    seqsQC = seqsRelabeled
        .branch{
            valid: it.sequence.count('X')*100 / it.sequence.size()  <= params.seqQC 
            invalid: it.sequence.count('X') *100/ it.sequence.size() > params.seqQC
        }.set { seqsFiltered}

    seqsFiltered.invalid.view { "INVALID >${it.header}" }
    seqsInvalidCount = seqsFiltered.invalid.count()
    writeFastaFromSeqsInvalid (seqsFiltered.invalid.map{record -> '>' + record.header + ',' + record.sequence}.collect(), "too_many_unknown_characters.fasta")
    
    //compare sequence and structure labels
    seqIDs =seqsFiltered.valid.map{tuple(it.header , it.sequence)}
    allSequencesCount = seqIDs.count()
    allSequencesCount.view{"Sequences to be aligned:" + it}

    if (params.dropSimilar){
        fullInputSeqsNum = cdHitCollapse.out.num.toInteger().sum().view{"Number of sequences before collapsing:" + it} 
    }else{
        fullInputSeqsNum = allSequencesCount
    }
    
    strucIDs = userStructures.map{strucL -> tuple(strucL.baseName, strucL.name)}

    strucIDs
        .join(seqIDs, remainder: true)
        .branch {
            modelNotFound:      it[1] == null
            sequenceNotFound:   it[2] == null
            modelFound:         true //it[1] != null || it[2] !=null
        }.set{joined}

    //joined.modelFound.view{"Provided model matched:" + it[0]}
    joined.modelFound.count().view{"Provided model matched:" + it}
    joined.sequenceNotFound.count().view {"Model without matching sequence:" +it}
    nosequenceModels = joined.sequenceNotFound.map{it[1]}


    //sumbit list of modelNotFound IDs to AFDB
    if (params.retrieve) {
        uniprotIDs = joined.modelNotFound.map{ it[0]}
        getAFmodels (uniprotIDs)
        protsFromAF = getAFmodels.out

        //update list
        getAFout = getAFmodels.out
            .map{it -> tuple(it.baseName, it.name)}

        getAFout
            .mix(strucIDs)
            .join(seqIDs, remainder: true)
            .branch {
                modelNotFound:      it[1] == null
                sequenceNotFound:   it[2] == null
                modelFound:         true //it[1] != null || it[2] !=null
            }.set{afjoined}
            //afjoined.modelFound.view{"Post AF Models matched:" + it[0]} 
            //afjoined.modelNotFound.view{"Missing after AF2DB search:" + it[0]}
        
    } else {
        getAFout = Channel.empty()
        protsFromAF= Channel.empty()
    }

    
    //sumbit list of modelNotFound sequences from AFDB or folder search to ESM Atlas
    if (params.model){
    
        if (params.retrieve) { missingModels =afjoined.modelNotFound }
        else {missingModels =joined.modelNotFound }

        missingModels
            .branch{
                toolong:    it[2].length() >= 400
                good:       it[2].length() < 400
            }.set{tomodel}

        tomodel.toolong.view{"Sequence to long to be predicted with ESM Atlas:" + it[0]}
        fetchEsmAtlasStructure(tomodel.good.map{tuple(it[0],it[2])},tomodel.good.count())

        protsFromESM = fetchEsmAtlasStructure.out

        //update list
        fetchEsmAtlasStructure.out
            .map{it -> tuple(it.baseName, it.name)}
            .mix(strucIDs,getAFout)
            .join(seqIDs, remainder: true)
            .branch {
                modelNotFound:      it[1] == null
                sequenceNotFound:   it[2] == null
                modelFound:         true //it[1] != null || it[2] !=null
            }.set{esmfjoined}
            //esmfjoined.modelFound.view{"Post ESMF Models matched:" + it[0]}      
            //esmfjoined.modelNotFound.view{"Missing after ESM Atlas search:" + it[0]}
        }else{
            protsFromESM = Channel.empty()
        }


    // assess which models could not be found in the folder, AFDB or ESM Atlas
    if (params.model) {
        finalMissingModels = esmfjoined.modelNotFound
        finalModelFound = esmfjoined.modelFound
    } else if (params.retrieve) { 
        finalMissingModels = afjoined.modelNotFound
        finalModelFound = afjoined.modelFound
    } else {
        finalMissingModels = joined.modelNotFound
        finalModelFound = joined.modelFound
    }
    
    structurelessCount=finalMissingModels.count()
    structurelessCount.view{"Missing models: " + it}
    foundSequencesCount = finalModelFound.count()
    
    writeFastaFromMissing(finalMissingModels.map{record -> ">"+ record[0] + ',' + record[2]}.collect(), 'structureless_seqs.fasta')
    structureless_seqs = writeFastaFromMissing.out.found

    //run local esmfold
    if (params.localModel){
        esmFolds(structureless_seqs)
        esmFoldsGate = esmFolds.out.gate 
        
        writeFastaFromSeqsValid (seqIDs.map{record -> ">"+ record[0]+ ',' + record[1]}.collect(),'seqs_to_align.fasta')

        foundSeqs = writeFastaFromSeqsValid.out.found
        missingQC (allSequencesCount, esmFoldsGate, params.strucQC)
    
    }else{
        missingQC (allSequencesCount, structurelessCount, params.strucQC)
        writeFastaFromFound(finalModelFound.map{record ->  ">"+ record[0] + ',' + record[2]}.collect(), 'seqs_to_align.fasta')
        foundSeqs = writeFastaFromFound.out.found
    }


    //subsetting
    if (params.useSubsets){
        userSubsetting(sequenceFastas,foundSeqs)
        subsets = userSubsetting.out.subsetSeqsForTC.collect().flatten()

    }else if (params.createSubsets){
        cdHitSubsetting(foundSeqs, params.createSubsets, params.minSubsetID, params.maxSubsetSize)
        subsets =  cdHitSubsetting.out.seq.collect().flatten()
    }else{
        subsets = foundSeqs 
    }

    subsets.branch{
        orphans:  it =~ /orphan/
        subsets: true
    }.set{subSeqs}

    seqsToAlign = subSeqs.subsets

    //submit to tcoffee
    runTcoffee(seqsToAlign, params.structures, params.tcoffeeParams, missingQC.out.gate)
    strucMsa =runTcoffee.out.msa.flatten()

    convertTCMsa('clustal', strucMsa)
    convertedMsa =  convertTCMsa.out.convertedSeqs.mix(structureless_seqs,subSeqs.orphans).collect()

    //create final alignment
    mergeMafft(convertedMsa, params.mafftParams, params.outName)
    finalMsa = mergeMafft.out.finalMsa


    //check if all sequences been aligned 
    attendance(foundSequencesCount,seqsInvalidCount,structurelessCount,finalMsa,fullInputSeqsNum)


    //map to dssp
    if (params.dssp){
        foundModels = userStructures.mix(protsFromAF,protsFromESM )
        runDssp(foundModels, missingQC.out.gate)
        dssps = runDssp.out.dsspout.collect()

        //relies on dssp being finished WAY before tcoffee alignment. 
        mapDsspRough("$params.outFolder/dssp", finalMsa)
        mappedFinalMsa = mapDsspRough.out.mmsa
    }

    //squeeze MSA
    if (params.squeeze){
        squeeze(finalMsa,mappedFinalMsa,params.squeeze,params.squeezePerc)
        squeezedMsa = squeeze.out.msa

        //map dssp to final msa
        mapDsspSqueeze("$params.outFolder/dssp", squeezedMsa)
        mappedFinalMsaSqueeze = mapDsspSqueeze.out.mmsa
    }else{squeezedMsa=finalMsa}


    //reorder final MSA
    if (params.reorder){
        reorder(squeezedMsa,sequenceFastas.collect(), params.reorder)
        reorderedFinalMsa = reorder.out.msaOrga
    }else{reorderedFinalMsa= squeezedMsa}


    //select output format for MSA
    if (params.convertMSA){
        convertFinalMsa(params.convertMSA, reorderedFinalMsa)

    }

}

workflow.onComplete {
    println "Pipeline completed at               : $workflow.complete"
    println "Time to complete workflow execution : $workflow.duration"
    println "Execution status                    : ${workflow.success ? 'Success' : 'Failed' }"
    println "Output folder                       : $params.outFolder"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    println "Details: \n ${workflow.errorReport}"
    
}

