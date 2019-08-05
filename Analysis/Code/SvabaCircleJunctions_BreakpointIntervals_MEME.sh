#!/bin/bash
source activate motifs
meme ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_BreakpointIntervals_41bp.fa -neg ~/Desktop/PalmTrees/Data/1MRandomSeq.fa -objfun de -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_BreakpointIntervals_vs_RandomSeqs_MEME/
meme ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Homology_Min5.fa -minw 5 -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Homology_Min5_MEME
meme ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Homology_Min5.fa -markov_order 3 -minw 5 -revcomp -nmotifs 5 -dna -oc ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Homology_Min5_MEME_MarkovOrder3
source deactivate motifs
