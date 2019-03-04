library(reshape2)
# You will want to set this to your relative path to load data and save the figure correctly.
setwd("~/Documents/lotterhosLab/OysterGenomeProject/popstructureOutliers/")
# Load the data from the relatedness data folder
rel <- read.csv("data/relatedness_data/relatedness_stats.relatedness",
                sep="\t",
                header=TRUE,
                stringsAsFactors=FALSE)

# Transform the data from long to wide format.
New.rel <- dcast(rel, INDV2 ~ INDV1, value.var = "RELATEDNESS_AJK")
# Set the row names to the first row due to them being in the first row.
rownames(new.rel) <- new.rel[,1]
new.rel <- new.rel[,-1]

# Color ramp
col = colorRampPalette(c("white", "red"))(30)

# If you set your working directory correctly, the plot should be saved directly to the figures folder.

png("figures/5relatedness/relatedness_plot.png", height=1500, width=1500)
par(mar = c(8,8,2,1))
image(1:90, 1:91, as.matrix(new.rel), axes = FALSE, xlab="", ylab="", col=col)
axis(1, at=1:90, las=2, labels=rownames(new.rel))
axis(2, at=1:90, las=2, labels=colnames(new.rel))
dev.off()
