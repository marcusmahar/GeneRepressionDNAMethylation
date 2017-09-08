#Set the working directory to the folder containing the datasheets
setwd("C:/Users/Marcus Mahar/Desktop/Lab/Lab Work/Sequencing/2nd miR Microarray/Data and Analysis/R Analysis/miR Targets/Down Regulated")

#read in necessary packages
library(xlsx)

#File Name
filename = "mmu-miR-9-5p"

#read in datasets - change to the appropriate miRNA
arraydata <- read.xlsx2("3 and 8 Hour Fold Analysis Edited.xls", sheetIndex = 1)
mirTargets <- read.csv(paste0(filename,".csv"), header = TRUE, na.strings = c("NA", ""))

#Grab the miR Targets
TargetScan <- tolower(mirTargets$TargetScan_Targets[!is.na(mirTargets$TargetScan_Targets)])
DianaTools <- tolower(mirTargets$DianaTools_Targets[!is.na(mirTargets$DianaTools_Targets)])
MetaCore <- tolower(mirTargets$MetaCore_Targets[!is.na(mirTargets$MetaCore_Targets)])

#Find the overlapping and uniquely predicted genes
predOverlap <- intersect(TargetScan, DianaTools)
allOverlap <- union(predOverlap, MetaCore)
TargetUnique <- unique(subset(TargetScan, !(TargetScan %in% predOverlap)))
DianaUnique <- unique(subset(DianaTools, !(DianaTools %in% predOverlap)))

#Union of all 3 datasets
TargetUnion <- union(MetaCore, union(TargetScan, DianaTools))

#Subset Targets that are in the array
ArrayTargets <- arraydata[tolower(arraydata$SYMBOL) %in% TargetUnion,]

#Identify targets with incorrect directionality
NotTargets <- unique(subset(ArrayTargets, as.numeric(as.character(X3.Hours)) < 1.1 & as.numeric(as.character(X8.Hours)) < 1.1))

#Eliminate NotTargets from Union dataset
PossTargets <- setdiff(TargetUnion, tolower(NotTargets$SYMBOL))
OverlapLikely <- setdiff(allOverlap, tolower(NotTargets$SYMBOL))

#Change outputs into dataframes and change colnames
OLDF <- data.frame(OverlapLikely)
LikelyTarDF <- data.frame(PossTargets)
MCDF <- data.frame(MetaCore)
OverlapDF <- data.frame(predOverlap)
TSDF <- data.frame(TargetUnique)
DTDF <- data.frame(DianaUnique)
colnames(OLDF) <- c("High_Confidence_Targets")
colnames(LikelyTarDF) <- c("Likely_Targets")
colnames(MCDF) <- c("MetaCore_Targets")
colnames(OverlapDF) <- c("Overlapping_Predictions")
colnames(TSDF) <- c("Unique_TargetScan_Predictions")
colnames(DTDF) <- c("Unique_DianaTools_Predictions")

#Print out the xls page in a new directory

write.xlsx2(OLDF, file = paste0(filename," Targets.xls"), col.names = TRUE, row.names = FALSE, append = FALSE, sheetName = "High Confidence Targets")
write.xlsx2(LikelyTarDF, file = paste0(filename," Targets.xls"), col.names = TRUE, row.names = FALSE, append = TRUE, sheetName = "Possible Targets")
write.xlsx2(MCDF, file = paste0(filename," Targets.xls"), col.names = TRUE, row.names = FALSE, append = TRUE, sheetName = "Target Lists" )
write.xlsx2(OverlapDF, file = paste0(filename," Targets.xls"), col.names = TRUE, row.names = FALSE, append = TRUE, sheetName = "Overlap")
write.xlsx2(TSDF, file = paste0(filename," Targets.xls"), col.names = TRUE, row.names = FALSE, append = TRUE, sheetName = "TargetScan")
write.xlsx2(DTDF, file = paste0(filename," Targets.xls"), col.names = TRUE, row.names = FALSE, append = TRUE, sheetName = "DianaTools")