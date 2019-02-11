###############################################################################
#
# File      : deeperLdAnalysis.R 
# History   : 12/20/2018 - Created by KBW
#             01/05/2019 - Altered for increased window sizes
#
###############################################################################

###############################################################################
#
# This script is the source code for the deeperLdAnalysis.html document that
# walks through the steps that we went through to find the window size that we
# wanted to use for the bigsnpr::snp_autoSVD() function to prune our SNPs.
#
###############################################################################

thinnedData <- readRDS("thinnedMatrixAndMetaData.rds")
allData <- readRDS("genotypeMatrix.rds")

ld.data.50bp <- read.csv("ldAnalysisFiles\\geno_ld_window_50-50.geno.ld", 
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE)
ld.data.200.250bp <- read.csv("ldAnalysisFiles\\geno_ld_window_200-250.geno.ld", 
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE)                         
ld.data.450.500bp <- read.csv("ldAnalysisFiles\\geno_ld_window_450-500.geno.ld", 
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE)                         
ld.data.4500.5000bp <- read.csv("ldAnalysisFiles\\geno_ld_window_4500-5000.geno.ld", 
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE)                         
ld.data.9500.10kbp <- read.csv("ldAnalysisFiles\\geno_ld_window_9500-10000.geno.ld", 
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE)                         
ld.data.49500.50kbp <- read.csv("ldAnalysisFiles\\geno_ld_window_49500-50000.geno.ld", 
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors=FALSE)                         
ld.data.99500.100kbp <- read.csv("ldAnalysisFiles\\geno_ld_window_99500-100000.geno.ld", 
                                 sep="\t",
                                 header=TRUE,
                                 stringsAsFactors=FALSE)                         
ld.data.499500.500kbp <- read.csv("ldAnalysisFiles\\geno_ld_window_499500-500000.geno.ld", 
                                  sep="\t",
                                  header=TRUE,
                                  stringsAsFactors=FALSE)                         

ld.vars <- list("50" = ld.data.50bp,
                "200" = ld.data.200.250bp, 
                "500" = ld.data.450.500bp, 
                "5000" = ld.data.4500.5000bp,
                "10000" = ld.data.9500.10kbp,
                "50000" = ld.data.49500.50kbp,
                "100000" = ld.data.99500.100kbp,
                "500000" = ld.data.499500.500kbp)

CHR <- c("NC_035789.1", "NC_035780.1", "NC_035781.1", "NC_035782.1", 
         "NC_035783.1", "NC_035784.1", "NC_035785.1", "NC_035786.1", 
         "NC_035787.1", "NC_035788.1")
col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
              '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
              '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
              '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
              '#ffffff', '#000000')
chr.num <- length(unique(ld.vars[[1]]$CHR))
for (i in length(ld.vars)) {
    if (length(unique(ld.vars[[i]]$CHR)) != chr.num) {
        stop("The number of chromosomes in the datasets do not match")
    } else {

    }
}

chr.index <- rep(NA, chr.num)
for (i in 1:chr.num){
    chr.index[i] <- substr(unique(ld.vars[[1]]$CHR)[i], start = 9, stop = 9)
}

ld.stats <- data.frame(matrix(NA, length(ld.vars), chr.num))
ld.stats[,1] <- as.numeric(names(ld.vars))
colnames(ld.stats)[1] <- "windowSize"
for (i in 1:length(ld.vars)) {
    ld.vars[[i]]$R.2 <- as.numeric(ld.vars[[i]]$R.2)
    for (j in 1:chr.num) {
        ld.stats[i, j+1] <- summary(ld.vars[[i]]$R.2[which(ld.vars[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])[4]
        colnames(ld.stats)[j+1] <- paste("NC_03578", chr.index[j], ".1", sep="")
    }
}

chr.length <- list("NC_035789.1" = 32650045,
                   "NC_035780.1" = 65668440,
                   "NC_035781.1" = 61752955,
                   "NC_035782.1" = 77061148,
                   "NC_035783.1" = 59691872,
                   "NC_035784.1" = 98698416,
                   "NC_035785.1" = 51258098,
                   "NC_035786.1" = 57830854,
                   "NC_035787.1" = 75944018,
                   "NC_035788.1" = 104168038)

snpsVsLength <- data.frame(as.integer(chr.length))
rownames(snpsVsLength) <- c("NC_035789.1", "NC_035780.1", "NC_035781.1", "NC_035782.1", "NC_035783.1", "NC_035784.1", "NC_035785.1", "NC_035786.1", "NC_035787.1", "NC_035788.1")


# Number of SNPs per chromosome length for the full thinned data
for (i in 1:10){
    print(c(names(chr.length)[which(as.integer(factor(names(chr.length))) == i)], 
            length(thinnedData$chromosome[which(thinnedData$chromosome == i)]) / as.integer(chr.length[which(as.integer(factor(names(chr.length))) == i)])))
    snpsVsLength[which(as.integer(factor(rownames(snpsVsLength))) == i),2] <- length(thinnedData$chromosome[which(thinnedData$chromosome == i)])
    snpsVsLength[i,3] <- ld.stats[8, i+1]
    snpsVsLength[which(as.integer(factor(rownames(snpsVsLength))) == i),4] <- length(thinnedData$chromosome[which(thinnedData$chromosome == i)]) / as.integer(chr.length[which(as.integer(factor(names(chr.length))) == i)])
    snpsVsLength[which(as.integer(factor(rownames(snpsVsLength))) == i),5] <- length(thinnedData$chromosome[which(thinnedData$chromosome == i)]) / length(allData$map$chromosome[which(allData$map$chromosome == i)])
}

exonPerLength <- c(0.000302358, 0.000533026, 0.000661296, 0.000611488, 0.000611021, 0.000558459, 0.000377209, 0.000405666, 0.000427934, 0.000360494)
snpsVsLength[,6] <- exonPerLength

colnames(snpsVsLength)  <- c("Length", "AmountOfSNPs", "LD", "SnpsPerLength", "ThinnedSnpsPerOriginalSnps", "exonsPerLength")

rownames(snpsVsLength) <- c("NC_035789.1", "NC_035780.1", "NC_035781.1", "NC_035782.1", "NC_035783.1", "NC_035784.1", "NC_035785.1", "NC_035786.1", "NC_035787.1", "NC_035788.1")



col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
              '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
              '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
              '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
              '#ffffff', '#000000')

# Can you plot # pruned SNPs vs. length
par(mar = c(5, 4, 4, 8) + 0.1)
plot(snpsVsLength$AmountOfSNPs ~ snpsVsLength$Length,
     pch = 2:(length(CHR)+1),
     cex = 2,
     col = col_vector[2:(length(CHR)+1)], 
     ylab = "Number of SNPs from Chromosome", 
     xlab = "Length of Chromosome")

legend("topright",
        inset = c(-0.35,0),
        legend = CHR, 
        pch = 2:(length(CHR)+1),
        xpd = TRUE,
        col = col_vector[2:(length(CHR)+1)])
dev.copy(png, "NoOfSNPsAndLengthPlot.png")
dev.off()

# pruned SNPs vs amount of LD at 500K

par(mar = c(5, 4, 4, 8) + 0.1)
plot(snpsVsLength$AmountOfSNPs ~ snpsVsLength$LD,
     pch = 2:(length(CHR)+1),
     cex = 2,
     col = col_vector[2:(length(CHR)+1)], 
     ylab = "Number of SNPs from Chromosome", 
     xlab = "LD @ 500K")

legend("topright",
        inset = c(-0.35,0),
        legend = CHR, 
        pch = 2:(length(CHR)+1),
        xpd = TRUE,
        col = col_vector[2:(length(CHR)+1)])
dev.copy(png, "NoOfSNPsVsLD.png")
dev.off()

# amount of LD at 500K vs length
par(mar = c(5, 4, 4, 8) + 0.1)
plot(snpsVsLength$LD ~ snpsVsLength$Length,
     pch = 2:(length(CHR)+1),
     cex = 2,
     col = col_vector[2:(length(CHR)+1)], 
     ylab = "LD @ 500K", 
     xlab = "Chromosome Length")

legend("topright",
        inset = c(-0.35,0),
        legend = CHR, 
        pch = 2:(length(CHR)+1),
        xpd = TRUE,
        col = col_vector[2:(length(CHR)+1)])
dev.copy(png, "LDVsChrLength.png")
dev.off()

# amount of LD at 500K vs exons/length (see thread on oyster genome slack channel)

par(mar = c(5, 4, 4, 8) + 0.1)
plot(snpsVsLength$LD ~ snpsVsLength$exonsPerLength,
     pch = 2:(length(CHR)+1),
     cex = 2,
     col = col_vector[2:(length(CHR)+1)], 
     ylab = "LD @ 500K", 
     xlab = "Exons per Chromosome Length")

legend("topright",
        inset = c(-0.35,0),
        legend = CHR, 
        pch = 2:(length(CHR)+1),
        xpd = TRUE,
        col = col_vector[2:(length(CHR)+1)])
dev.copy(png, "LDVsExonsPerChrLength.png")
dev.off()

###########################################

par(mar = c(5, 4, 4, 8) + 0.1)
plot(snpsVsLength$ThinnedSnpsPerOriginalSnps ~ snpsVsLength$LD,
     pch = 2:(length(CHR)+1),
     cex = 2,
     col = col_vector[2:(length(CHR)+1)], 
     ylab = "(# SNPs left after pruning)/(# original SNPs) @ chromosome", 
     xlab = "LD @ 500K")

legend("topright",
        inset = c(-0.35,0),
        legend = CHR, 
        pch = 2:(length(CHR)+1),
        xpd = TRUE,
        col = col_vector[2:(length(CHR)+1)])
dev.copy(png, "SnpsPerOriginalSnpsVsLD.png")
dev.off()