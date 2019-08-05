#!/bin/bash
source activate motifs
meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_OnlyPalmTreeTx_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_OnlyNoPalmTreeTx_BreakpointSequences.fa -objfun de -revcomp -nmotifs 10 -dna -oc ~/Desktop/PalmTrees/Results/MEME/UnionPTTxvsNonPTTx

meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_PTvsNoPT_BreakpointSequences.fa -revcomp -nmotifs 10 -dna -oc ~/Desktop/PalmTrees/Results/MEME/UnionAllBreakpoints
meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_PTvsNoPT_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Data/1MRandomSeq_61bp.fa -objfun de -revcomp -nmotifs 10 -dna -oc ~/Desktop/PalmTrees/Results/MEME/UnionAllBreakpointsVsRandom

meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_CircleSeq_genome-genome_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_CircleSeq_non-genome-genome_BreakpointSequences.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/MEME/Union_CircleSeq_genome-genome
meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_CircleSeq_circle-genome_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_CircleSeq_non-circle-genome_BreakpointSequences.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/MEME/Union_CircleSeq_circle-genome
meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_CircleSeq_circle-circle_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_CircleSeq_non-circle-circle_BreakpointSequences.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/MEME/Union_CircleSeq_circle-circle

meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_WGSCircles_genome-genome_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_WGSCircles_non-genome-genome_BreakpointSequences.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/MEME/Union_WGSCircles_genome-genome
meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_WGSCircles_circle-genome_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_WGSCircles_non-circle-genome_BreakpointSequences.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/MEME/Union_WGSCircles_circle-genome
meme ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_WGSCircles_circle-circle_BreakpointSequences.fa -neg ~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/Union_WGSCircles_non-circle-circle_BreakpointSequences.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/MEME/Union_WGSCircles_circle-circle

meme ~/Desktop/PalmTrees/Results/Figures/Homologies/SvabaTx_Homologies.fa -minw 5 -revcomp -nmotifs 10 -dna -oc ~/Desktop/PalmTrees/Results/Figures/Homologies/SvabaTx_Homologies-MEME
meme ~/Desktop/PalmTrees/Results/Figures/Homologies/SvabaTx_Homologies.fa -minw 5 -revcomp -nmotifs 10 -dna -markov_order 3 -oc ~/Desktop/PalmTrees/Results/Figures/Homologies/SvabaTx_Homologies-MEME-Markov3

source deactivate
