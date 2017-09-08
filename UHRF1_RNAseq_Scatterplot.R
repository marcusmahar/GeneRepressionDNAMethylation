#Set the working directory
setwd("C:/Users/Marcus Mahar/Desktop/Lab/Lab Work/Sequencing/UHRF1_RNAseq")

#Read in the data

plotdata2 <- read.csv("shComparison_Uninjured_Avg_Counts_2.csv", header = TRUE)
sigdown2 <- subset(plotdata2, plotdata2$padj < .1 & plotdata2$log2FoldChange < 0)
sigup2 <- subset(plotdata2, plotdata2$padj < .1 & plotdata2$log2FoldChange > 0)


#Pick the gene to mark in the graph
Rest <- c("ENSMUSG00000029249")

#Set the data color
plotdata2$color <- "black"
plotdata2$color[plotdata2$padj < .1 & plotdata2$log2FoldChange > 0] <- "green"
plotdata2$color[plotdata2$padj < .1 & plotdata2$log2FoldChange < 0] <- "red"
#plotdata$color[plotdata$ENSEMBL %in% Rest] <- "orange"

#Call the device driver
pdf(file = "UHRF1_Scatterplot.pdf", onefile = TRUE)

#Plot the data
plot(plotdata2$shControl_Uninjured_Average,plotdata2$shUHRF1_Uninjured_Average, log = "xy", pch = 19, cex = .1, 
     xlab = "shControl Counts", ylab = "shUHRF1 Counts", col=plotdata$color)
dev.off()


