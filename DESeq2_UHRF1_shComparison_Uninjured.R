#Set the Working Directory

directory <- "C:/Users/Marcus Mahar/Desktop/Lab/Lab Work/Sequencing/UHRF1_RNAseq/HTseq-count Files"
setwd(directory)

#Call all necessary packages

library("DESeq2")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("grDevices")
library("org.Mm.eg.db")
library("genefilter")

#Set the matrix variables
sampleFiles <- c("shControl_Uninjured_R1.counts", "shControl_Uninjured_R2.counts", "shControl_Uninjured_R3.counts", "shUHRF1_Uninjured_R1.counts", 
                 "shUHRF1_Uninjured_R2.counts", "shUHRF1_Uninjured_R3.counts")

sampleNames <- c("shControl_Uninjured_R1", "shControl_Uninjured_R2", "shControl_Uninjured_R3", "shUHRF1_Uninjured_R1", "shUHRF1_Uninjured_R2", 
                 "shUHRF1_Uninjured_R3")

sampleCondition <- c("shControl", "shControl", "shControl", "shUHRF1", "shUHRF1", "shUHRF1")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

#Diffential Expression Equation
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = c("shControl", "shUHRF1"))
ddsHT <- ddsHTSeq[rowMeans(counts(ddsHTSeq)) > 1, ]
dds <- DESeq(ddsHT)
res <- results(dds)

#Rlog transformation of values
rld <- rlog(dds)

#MA Plot
pdf("RNAseq_MAPlot_shComparison_Uninjured.pdf")
plotMA(res, main="MA Plot", ylim = c(-5,5))
dev.off()
#identify(res$baseMean, res$log2FoldChange)

#Dispersion Plot
pdf("RNAseq_Dispersion_shComparison_Uninjured.pdf")
plotDispEsts(dds, ylim = c(1e-6, 1e1), main = "Dispersion Plot")
dev.off()

#P Value Histogram
pdf("RNAseq_PvalHisto_shComparison_Uninjured.pdf")
hist(res$pvalue, breaks=20, col="grey")
dev.off()

#Samples Heatmap - can also be done with particular genes
pdf("RNAseq_SampleHeatmap_shComparison_Uninjured.pdf")
SampleDists <- dist(t(assay(rld)))
SampleDistMatrix <- as.matrix(SampleDists)
colnames(SampleDistMatrix) <- NULL
colours = colorRampPalette(rev(brewer.pal(6, "Blues")))(255)
heatmap.2(SampleDistMatrix, trace = "none", col = colours, cexRow = .3, cexCol = 1.5, key.title = "", density.info = "none")
dev.off()

#PCA Plot
pdf("RNAseq_PCAPlot_shComparison_Uninjured.pdf")
plotPCA(rld, intgroup = c("condition"))
dev.off()

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, "UHRF1_RNAseq_shComparison_Uninjured.csv")