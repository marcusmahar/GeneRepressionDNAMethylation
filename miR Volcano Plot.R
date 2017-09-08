#Set the working directory to the folder containing the datasheets
setwd("C:/Users/Marcus Mahar/Desktop/Lab/Lab Work/Sequencing/2nd miR Microarray/Data and Analysis/R Analysis")

#read in necessary packages
library(xlsx)
library(calibrate)

#read in datasets
data <- read.csv("miRNA Array Full Results.csv")
mir1H <- read.xlsx("Differentially Regulated.xlsx", sheetIndex = 1)
mir3H <- read.xlsx("Differentially Regulated.xlsx", sheetIndex = 2)
mir8H <- read.xlsx("Differentially Regulated.xlsx", sheetIndex = 3)

#transform p-values to log scale
pval1H <- -log10(data$X1H.p.val)
pval3H <- -log10(data$X3H.p.val)
pval8H <- -log10(data$X8H.p.val)

#Select down regulated miRs
mir1Hdown <- mir1H[mir1H$Fold.Change.vs..Control < 0, ]
mir3Hdown <- mir3H[mir3H$Fold.Change.vs..Control < 0, ]
mir8Hdown <- mir8H[mir8H$Fold.Change.vs..Control < 0, ]

#Select just the down regulated miR names
mir1Hdownnames <- mir1Hdown$Transcript.ID
mir3Hdownnames <- mir3Hdown$Transcript.ID
mir8Hdownnames <- mir8Hdown$Transcript.ID

#Select upregulated miRs
mir1Hup <- mir1H[mir1H$Fold.Change.vs..Control > 0, ]
mir3Hup <- mir3H[mir3H$Fold.Change.vs..Control > 0, ]
mir8Hup <- mir8H[mir8H$Fold.Change.vs..Control > 0, ]

#Select just the upregulated miR names
mir1Hupnames <- mir1Hup$Transcript.ID
mir3Hupnames <- mir3Hup$Transcript.ID
mir8Hupnames <- mir8Hup$Transcript.ID

#Compare full dataset to down names
sigdown1H <- data$Transcript.ID.Array.Design. %in% mir1Hdownnames
sigdown3H <- data$Transcript.ID.Array.Design. %in% mir3Hdownnames
sigdown8H <- data$Transcript.ID.Array.Design. %in% mir8Hdownnames

#Compare full dataset to up names
sigup1H <- data$Transcript.ID.Array.Design. %in% mir1Hupnames
sigup3H <- data$Transcript.ID.Array.Design. %in% mir3Hupnames
sigup8H <- data$Transcript.ID.Array.Design. %in% mir8Hupnames

#miRNAs to be labeled with text in the plots
mir9 <- data$Transcript.ID == "mmu-miR-9-3p"

#Call the device driver
pdf(file = "miR_Volcano_test.pdf", onefile = TRUE)

#Global plot settings - Change to format location/height
par(mfcol = c(1,3)) 
par(oma = c(2,3,2,0))
par(mar = c(4,2,4,2))

#1 HPI Plot
with(data, plot(X1H.Log2.FC, pval1H, pch = 20, xlab = "1 HPI Log2 Fold Change", ylab = "-log10(p-value)", col = "gray", ylim = c(0, 5.5), xlim = c(-3.5,3.5)))
points(data$X1H.Log2.FC[sigup1H], pval1H[sigup1H], pch = 20, col = "green")
points(data$X1H.Log2.FC[sigdown1H], pval1H[sigdown1H], pch = 20, col = "red")
points(data$X1H.Log2.FC[mir9], pval1H[mir9], pch = 20, col = "blue")
#textxy(data$X1H.Log2.FC[mir9], pval1H[mir9], labs = "miR-9-3p", cex = .8, offset = .6)

#3 HPI Plot
with(data, plot(X3H.Log2.FC, pval3H, pch = 20, xlab = "3 HPI Log2 Fold Change", col = "gray", ylim = c(0, 5.5), xlim = c(-3.5,3.5)))
points(data$X3H.Log2.FC[sigup3H], pval3H[sigup3H], pch = 20, col = "green")
points(data$X3H.Log2.FC[sigdown3H], pval3H[sigdown3H], pch = 20, col = "red")
points(data$X3H.Log2.FC[mir9], pval3H[mir9], pch = 20, col = "blue")
#textxy(data$X3H.Log2.FC[mir9], pval3H[mir9], labs = "miR-9-3p", cex = .8, offset = .6)

#8 HPI Plot
with(data, plot(X8H.Log2.FC, pval8H, pch = 20, xlab = "8 HPI Log2 Fold Change", col = "gray", ylim = c(0, 5.5), xlim = c(-3.5,3.5)))
points(data$X8H.Log2.FC[sigup8H], pval8H[sigup8H], pch = 20, col = "green")
points(data$X8H.Log2.FC[sigdown8H], pval8H[sigdown8H], pch = 20, col = "red")
points(data$X8H.Log2.FC[mir9], pval8H[mir9], pch = 20, col = "blue")
#textxy(data$X8H.Log2.FC[mir9], pval8H[mir9], labs = "miR-9-3p", cex = .8, offset = .6)

#Turn off device
dev.off()
     
