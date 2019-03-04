rel <- read.csv("../data/relatedness_data/relatedness_stats.relatedness",
                sep="\t",
                header=TRUE,
                stringsAsFactors=FALSE)

library(reshape2)

new.rel <- dcast(rel, INDV2 ~ INDV1, value.var = "RELATEDNESS_AJK")
rownames(new.rel) <- new.rel[,1]
new.rel <- new.rel[,-1]

col = colorRampPalette(c("white", "red"))(30)

png("relatednessPlot.png", height=1500, width=1500)
par(mar = c(8,8,2,1))
image(1:90, 1:91, as.matrix(new.rel), axes = FALSE, xlab="", ylab="", col=col)
axis(1, at=1:90, las=2, labels=rownames(new.rel))
axis(2, at=1:90, las=2, labels=colnames(new.rel))
dev.off()
