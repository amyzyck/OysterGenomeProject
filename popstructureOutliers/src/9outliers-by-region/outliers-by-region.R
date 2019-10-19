chr.names <- c("NC_035780.1",
               "NC_035781.1",
               "NC_035782.1",
               "NC_035783.1",
               "NC_035784.1",
               "NC_035785.1",
               "NC_035786.1",
               "NC_035787.1",
               "NC_035788.1",
               "NC_035789.1")

fst.outliers <- read.table("FstOutliers_Unrelated.bed", 
                           header = TRUE, 
                           stringsAsFactors=FALSE, 
                           sep="")
exon.data <- read.table("Outliers_exon.txt", 
                        header = FALSE, 
                        stringsAsFactors=FALSE, 
                        sep="\t")
cds.data <- read.table("Outliers_CDS.txt", 
                       header = FALSE, 
                       stringsAsFactors=FALSE, 
                       sep="\t")
gene.data <- read.table("Outliers_gene.txt", 
                        header = FALSE, 
                        stringsAsFactors=FALSE, 
                        sep="\t")
inter.data <- read.table("Outliers_Intergenic.txt", 
                         header = FALSE, 
                         stringsAsFactors=FALSE, 
                         sep="\t")
utr.data <- read.table("Outliers_UTR.txt", 
                       header = FALSE, 
                       stringsAsFactors=FALSE, 
                       sep="\t")

all.regions <- list(exon.data, 
                    cds.data,
                    gene.data, 
                    inter.data,
                    utr.data,
                    fst.outliers)

outliers.by.region <- data.frame("CHR"=chr.names,
                                 "Exon"=NA,
                                 "CDS"=NA,
                                 "Gene"=NA,
                                 "Intergenic"=NA,
                                 "UTR"=NA,
                                 "Total"=NA)

for(i in 1:length(chr.names)) {
    for(j in 1:length(all.regions)) {
        temp <- all.regions[[j]]
        outliers.by.region[i, j+1] <- length(which(temp[,1] == chr.names[i]))
    }
}

write.table(outliers.by.region,
            "OutliersByChrAndRegion_test.txt", 
            row.names = FALSE)

data <- read.table("OutliersByChrAndRegion_test.txt", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

par(mar = c(7.1, 4.1, 4.1, 2.1))
barplot(data$Total, 
        names = data$CHR,
        main = "Outlier SNPs by Chromosome and Region",
        las = 2, 
        cex.names = 0.75,
        ylim = c(0, 8000), 
        border = FALSE, 
        col = rgb(0.5, 0.5, 0.5, 0.5))
barplot(data$Gene, 
        names = data$CHR, 
        las = 2, 
        cex.names = 0.75,
        ylim = c(0, 8000), 
        border = FALSE, 
        col = rgb(0, 0, 0, 0.7),
        add = TRUE)
barplot(data$Intergenic, 
        ylab = "",
        las = 2,
        border = FALSE, 
        col = rgb(1, 0, 0, 0.4),
        add = TRUE)

legend('topright', 
        bty = 'n',
        legend = c("All Outliers", "Gene (Exons, UTR, CDS, etc.)", "Intergenic"), 
        fill = c(rgb(0.5, 0.5, 0.5, 0.5), 
                 rgb(0, 0, 0, 0.7), 
                 rgb(0.5, 0, 0)))
dev.copy(pdf, "OutliersByChrAndRegion_All.pdf")
dev.off()


barplot(data$Gene, 
        names = data$CHR,
        main = "Outliers SNPs by Chromosome and Region (Gene)",
        las = 2, 
        cex.names = 0.75,
        border = FALSE, 
        col = rgb(0.5, 0.5, 0.5, 0.2))
barplot(data$Exon, 
        names = data$CHR, 
        las = 2, 
        cex.names = 0.75,
        border = FALSE, 
        col = rgb(0, 0, 0, 0.7),
        add = TRUE)
barplot(data$CDS, 
        names = data$CHR, 
        las = 2, 
        cex.names = 0.75,
        border = FALSE, 
        col = rgb(0.5, 0, 0),
        add = TRUE)
barplot(data$UTR, 
        names = data$CHR, 
        las = 2, 
        cex.names = 0.75,
        border = FALSE, 
        col = rgb(0.5, 0.5, 0.5),
        add = TRUE)

legend('topright', 
        bty = 'n',
        legend = c("Gene", "Exon", "CDS", "UTR"), 
        fill = c(rgb(0.5, 0.5, 0.5, 0.2), 
                 rgb(0, 0, 0, 0.7), 
                 rgb(0.5, 0, 0), 
                 rgb(0.5, 0.5, 0.5)))
dev.copy(pdf, "OutliersByChrAndRegion_Gene.pdf")
dev.off()


###############################################################################
# Contingency Table 
###############################################################################

gene.total.portions <- data$Gene / data$Total
gene.ind.portions <- data$Gene / sum(data$Gene)

inter.total.portions <- data$Intergenic / data$Total
inter.ind.portions <- data$Intergenic / sum(data$Intergenic)

total.table <- data.frame("Chromosome" = data$CHR, 
                          "Gene" = data$Gene,
                          "GenePercentOfTotal" = gene.total.portions,
                          "GenePercentOfRegion" = gene.ind.portions,
                          "Intergenic" = data$Intergenic,
                          "InterPercentOfTotal" = inter.total.portions,
                          "InterPercentOfRegion" = inter.ind.portions,
                          "Total" = data$Total)

for (i in 2:ncol(total.table)) {
    total.table[11, i] <- sum(total.table[1:10, i])
}

ChrLength <- c(65668440, 61752955, 77061148, 59691872, 98698416, 51258098, 57830854, 75944018, 104168038, 32650045)
total.table$ChrLength <- c(ChrLength, sum(ChrLength))


