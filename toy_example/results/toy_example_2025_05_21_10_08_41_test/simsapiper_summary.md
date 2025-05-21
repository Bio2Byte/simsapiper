# Summary file for SIMSApiper 
Find detailed explations for each step on [GitHub](https://github.com/Bio2Byte/simsapiper) 

All output files can be found here:  /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test

The command executed was:  nextflow run simsapiper.nf -profile standard,withdocker --data /Users/sophie/workspace/simsapiper/toy_example/data --magic --outFolder /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test
# Input sequences files 
* No. of sequences in inputfile(s):  13
* Input files can be found in  /Users/sophie/workspace/simsapiper/toy_example/data/seqs
* SIMSApiper found these files:  toy_example_seqs.fasta
# 1 Data preparation
## 1.2 Data reduction with CD-Hit
* Similarity cutoff for CD-Hit:  90
* Included selected sequences to respresenting sequence file:  false
* No. of sequences after collapsing:  11
* Find collapsed sequences here: 
    * /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/seqs/CD-HIT_collapse/removed_converted_toy_example_seqs_90_perc_similar.fasta
* Find representing and represented sequences here: 
    * /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/seqs/CD-HIT_collapse/converted_toy_example_seqs_90_collapsed.clstr
## 1.3 Invalid sequences 
* Sequences removed because they contained more then 5 % non-standard or unresolved amino acids:  0
# 2 Structural information
## 2.1 Identify missing models
* Models found in /Users/sophie/workspace/simsapiper/toy_example/data/structures
* Models found at start of the pipeline:  0
## 2.2 Search AFDB:  true
* Models found in AlphaFold Database:  11
## 2.4 Use ESM Atlas:  true
* Models generated with ESMAtlas:  0
## 2.4 Alternative: Use local ESMfold:  false
## 2.6 No. of sequences not matched to a model:  0
# 3 Generation of subsets
## 3.1 Subsets were generated from input files: false
## 3.2 Automatically generated subsets : true 
* Subsets were generated 30 % sequence similarity cutoff
* Find which subset includes which sequence here:  /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/seqs/CD-HIT_for_t-coffee/seqs_to_align_30_clustered.clstr
# 4 Align each subsets with T-Coffee
* Cleaned, matched sequences to submit to T-coffee: 
    *  /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/seqs/CD-HIT_for_t-coffee/Subset0.fasta
* Additional T-Coffee parameters:  false
# 5 Align subset alignments with MAFFT
* Alignment of subset alignments, structureless and orphan sequences  /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/msas/merged_magicMSA.fasta
* Additional MAFFT parameters:  false
# 6 Run DSSP: true
* DSSP file can be found in /Users/sophie/workspace/simsapiper/toy_example/data/dssp
# 7 Improve MSA 
## 7.1 Map DSSP to MSA
* DSSP codes mapped to merged alignment: /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/msas/dssp_merged_magicMSA.fasta
## 7.2 Squeeze MSA towards conserved secondary structure elements
* Categories selected for squeezing:  H,E
* Threshold for region to be considered conserved:  80
* Squeezed alignment:  /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/msas/squeezed_merged_magicMSA.fasta
## 7.3 Map DSSP to squeezed MSA
* DSSP codes mapped to merged alignment: /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/msas/dssp_squeezed_merged_magicMSA.fasta
## Gap reduction
* The inital merged alignment has  507  gaps.
* The the squeezed alignment has  463  gaps.
* The squeezing step removed  44  gaps from the alignment.
## Estimated sequence identity
* The final alignment contains 15 conserved positions and 213 gapless positions of 293 total positions.
* The average pairwise sequence identity in the final alignment is 32.89%. Gaps were treated as mismatches.
# 8 Reorder MSA 
* Alignment reorded by:  true
* Reordered alignment:  /Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/msas/reordered_squeezed_merged_magicMSA.fasta
# Final MSA file
/Users/sophie/workspace/simsapiper/toy_example/results/toy_example_2025_05_21_10_08_41_test/magicMSA_simsa.fasta

# Execution time
* Total runtime (without queue time): 3 minutes
* Total runtime runTcoffee only: 2 minutes
