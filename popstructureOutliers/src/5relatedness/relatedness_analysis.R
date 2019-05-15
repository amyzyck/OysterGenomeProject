library(reshape2)

# Load the data from the relatedness data folder
rel <- read.csv("data/relatedness_data/relatedness_stats.relatedness",
                sep="\t",
                header=TRUE,
                stringsAsFactors=FALSE)

# Transform the data from long to wide format.
new.rel <- dcast(rel, INDV2 ~ INDV1, value.var = "RELATEDNESS_AJK")
# Set the row names to the first row due to them being in the first row.
rownames(new.rel) <- new.rel[,1]
new.rel <- new.rel[,-1]

# Color ramp
col = colorRampPalette(c("white", "red"))(30)

# If you set your working directory correctly, the plot should be saved directly
# to the figures folder.
png("figures/5relatedness/relatedness_plot_allPops.png", height=1500, width=1500)
par(mar = c(8,8,2,1))
image(1:90, 1:91, as.matrix(new.rel), axes = FALSE, xlab="", ylab="", col=col)
axis(1, at=1:90, las=2, labels=rownames(new.rel))
axis(2, at=1:90, las=2, labels=colnames(new.rel))
dev.off()

################################################################################
# Relatedness Analysis for the Atlantic selection populations subset
################################################################################

rel <- read.csv("data/relatedness_data/relatedness_stats_atlanticSelectionPops.relatedness", 
                sep="\t",
                header=TRUE,
                stringsAsFactors=FALSE)

library(reshape2)

new.rel <- dcast(rel, INDV2 ~ INDV1, value.var = "RELATEDNESS_AJK")
rownames(new.rel) <- new.rel[,1]
new.rel <- new.rel[,-1]

col = colorRampPalette(c("white", "red"))(30)

png("relatednessPlot_atlanticSelectionPops.png", height=750, width=750)
par(mar = c(8,8,2,1))
image(1:24, 1:24, as.matrix(new.rel), axes = FALSE, main = "Relatedness: Atlantic Selection Populations", xlab="", ylab="", col=col)
axis(1, at=1:24, las=2, labels=rownames(new.rel))
axis(2, at=1:24, las=2, labels=colnames(new.rel))
dev.off()

################################################################################
# Relatedness Analysis for the Atlantic wild populations subset
################################################################################

rel <- read.csv("data/relatedness_data/relatedness_stats_atlanticWildPops.relatedness", 
                sep="\t",
                header=TRUE,
                stringsAsFactors=FALSE)

library(reshape2)

new.rel <- dcast(rel, INDV2 ~ INDV1, value.var = "RELATEDNESS_AJK")
rownames(new.rel) <- new.rel[,1]
new.rel <- new.rel[,-1]

col = colorRampPalette(c("white", "red"))(30)

png("relatednessPlot_atlanticWildPops.png", height=950, width=950)
par(mar = c(8,8,2,1))
image(1:30, 1:30, as.matrix(new.rel), main = "Relatedness: Atlantic Wild Populations", axes = FALSE, xlab="", ylab="", col=col)
axis(1, at=1:30, las=2, labels=rownames(new.rel))
axis(2, at=1:30, las=2, labels=colnames(new.rel))
dev.off()
