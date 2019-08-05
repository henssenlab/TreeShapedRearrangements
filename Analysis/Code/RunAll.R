setwd("~/Desktop/PalmTrees/Analysis/Code/")

source("~/Desktop/PalmTrees/Analysis/Code/ParseCopyNumberData.R")
source("~/Desktop/PalmTrees/Analysis/Code/ParseCircles.R.R")

print("Palm Tree Calling")
source("~/Desktop/PalmTrees/Analysis/Code/CallPalmTrees.R")
source("~/Desktop/PalmTrees/Analysis/Code/CircosPlots.R")
source("~/Desktop/PalmTrees/Analysis/Code/CircosPlotsPro.R")

print("Brekpoint Sequences")
source("~/Desktop/PalmTrees/Analysis/Code/GenerateBreakpointSequences.R")
system("~/Desktop/PalmTrees/Analysis/Code/MemeAnalysis.sh")

print("General Palm Tree Statistics")
source("~/Desktop/PalmTrees/Analysis/Code/PalmTreeDensity.R")
source("~/Desktop/PalmTrees/Analysis/Code/PalmTreeGenes.R")
source("~/Desktop/PalmTrees/Analysis/Code/GeneralPalmtreeStatistics.R")
source("~/Desktop/PalmTrees/Analysis/Code/ClinicalData.R")

# print("Expression Analysis")
source("~/Desktop/PalmTrees/Analysis/Code/ParseRNASeqData.R")
source("~/Desktop/PalmTrees/Analysis/Code/ExpressionAnalysis.R")
 
# print("Overlap with Circles")
print("CS Ecc")
source("~/Desktop/PalmTrees/Analysis/Code/OverlapCircleSeqEccCircleCalls.R")
print("CS Ec")
source("~/Desktop/PalmTrees/Analysis/Code/OverlapCircleSeqEcCircleCalls.R")
print("WGS Ecc")
source("~/Desktop/PalmTrees/Analysis/Code/OverlapWGSCirclesEccCircleCalls.R")
print("WGS Ec")
source("~/Desktop/PalmTrees/Analysis/Code/OverlapWGSCirclesEcCircleCalls.R")
print("Merging")
source("~/Desktop/PalmTrees/Analysis/Code/MergeWGSCSOverlapPlots.R")
source("~/Desktop/PalmTrees/Analysis/Code/PalmTreeStackPlot.R")

# print("Pediatric Pan Cancer Data")
source("~/Desktop/PalmTrees/Analysis/Code/CallPalmTreesPedPanCan.R")
source("~/Desktop/PalmTrees/Analysis/Code/PedPanCanStats.R")

# print("Cell Lines")
source("~/Desktop/PalmTrees/Analysis/Code/CallPalmTreesCellLines.R")

