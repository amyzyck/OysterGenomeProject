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

load.ld.data <- function (path = NULL) {
    
    if (!dir.exists(path)) {
        stop("ERROR: The path that you have provided is not a directory")
    } else {
        ld.files <- paste(path, list.files(path), sep = "")
    }

    ld.data.50bp <- read.csv(ld.files[which(grepl("_50-50", ld.files, fixed=TRUE) == TRUE)], 
                             sep="\t",
                             header=TRUE,
                             stringsAsFactors=FALSE)

    ld.data.200.250bp <- read.csv(ld.files[which(grepl("_200-250", ld.files, fixed=TRUE) == TRUE)], 
                                  sep="\t",
                                  header=TRUE,
                                  stringsAsFactors=FALSE)

    ld.data.450.500bp <- read.csv(ld.files[which(grepl("_450-500", ld.files, fixed=TRUE) == TRUE)], 
                                  sep="\t",
                                  header=TRUE,
                                  stringsAsFactors=FALSE)

    ld.data.4500.5000bp <- read.csv(ld.files[which(grepl("_4500-5000", ld.files, fixed=TRUE) == TRUE)], 
                                    sep="\t",
                                    header=TRUE,
                                    stringsAsFactors=FALSE)

    ld.data.9500.10kbp <- read.csv(ld.files[which(grepl("_9500-10000", ld.files, fixed=TRUE) == TRUE)], 
                                    sep="\t",
                                    header=TRUE,
                                    stringsAsFactors=FALSE)

    ld.data.49500.50kbp <- read.csv(ld.files[which(grepl("_49500-50000", ld.files, fixed=TRUE) == TRUE)], 
                                    sep="\t",
                                    header=TRUE,
                                    stringsAsFactors=FALSE)

    ld.data.99500.100kbp <- read.csv(ld.files[which(grepl("_99500-100000", ld.files, fixed=TRUE) == TRUE)], 
                                     sep="\t",
                                     header=TRUE,
                                     stringsAsFactors=FALSE)

    ld.data.499500.500kbp <- read.csv(ld.files[which(grepl("500000", ld.files, fixed=TRUE) == TRUE)], 
                                      sep="\t",
                                      header=TRUE,
                                      stringsAsFactors=FALSE)

    ld.vars <- list("50" = ld.data.50bp,
                    "200" = ld.data.200.250bp, 
                    "450" = ld.data.450.500bp, 
                    "4500" = ld.data.4500.5000bp,
                    "9500" = ld.data.9500.10kbp,
                    "49500" = ld.data.49500.50kbp,
                    "99500" = ld.data.99500.100kbp,
                    "499500" = ld.data.499500.500kbp)
    return(ld.vars)
}
    

deep.ld.analysis <- function (ld.path = NULL,
                              thinned.snps = NULL,
                              full.data = NULL,
                              subset = NULL) {

    results.path <- file.path("figures", "1LD_analysis", "deeperLdAnalysis", subset)
    if (!dir.exists(results.path)) {
        dir.create(results.path, recursive = TRUE)
    }

    # Load the LD calculations from all window sizes
    ld.vars <- load.ld.data(path = ld.path)

    thinned.data <- readRDS(thinned.snps)

    allData <- readRDS(full.data)
    
    CHR <- c("NC_035780.1", "NC_035781.1", "NC_035782.1", "NC_035783.1", 
             "NC_035784.1", "NC_035785.1", "NC_035786.1", "NC_035787.1", 
             "NC_035788.1", "NC_035789.1")
     
    col_vector <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
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
    chr.index <- sort(chr.index)

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

    chr.length <- list("NC_035780.1" = 65668440,
                       "NC_035781.1" = 61752955,
                       "NC_035782.1" = 77061148,
                       "NC_035783.1" = 59691872,
                       "NC_035784.1" = 98698416,
                       "NC_035785.1" = 51258098,
                       "NC_035786.1" = 57830854,
                       "NC_035787.1" = 75944018,
                       "NC_035788.1" = 104168038, 
                       "NC_035789.1" = 32650045)

    snpsVsLength <- data.frame(as.integer(chr.length))
    rownames(snpsVsLength) <- CHR


# Number of SNPs per chromosome length for the full thinned data
for (i in 1:10){
    print(c(names(chr.length)[which(as.integer(factor(names(chr.length))) == i)], 
            length(thinned.data$Chr[which(thinned.data$Chr == i)]) / as.integer(chr.length[which(as.integer(factor(names(chr.length))) == i)])))
    snpsVsLength[which(as.integer(factor(rownames(snpsVsLength))) == i),2] <- length(thinned.data$Chr[which(thinned.data$Chr == i)])
    snpsVsLength[i,3] <- ld.stats[8, i+1]
    snpsVsLength[which(as.integer(factor(rownames(snpsVsLength))) == i),4] <- length(thinned.data$Chr[which(thinned.data$Chr == i)]) / as.integer(chr.length[which(as.integer(factor(names(chr.length))) == i)])
    snpsVsLength[which(as.integer(factor(rownames(snpsVsLength))) == i),5] <- length(thinned.data$Chr[which(thinned.data$Chr == i)]) / length(allData$Chr[which(allData$Chr == i)])
}

exonPerLength <- c(0.000533026, 
                   0.000661296, 
                   0.000611488, 
                   0.000611021, 
                   0.000558459, 
                   0.000377209, 
                   0.000405666, 
                   0.000427934, 
                   0.000360494, 
                   0.000302358)
snpsVsLength[,6] <- exonPerLength

colnames(snpsVsLength)  <- c("Length", 
                             "AmountOfSNPs", 
                             "LD", 
                             "SnpsPerLength", 
                             "ThinnedSnpsPerOriginalSnps", 
                             "exonsPerLength")

rownames(snpsVsLength) <- CHR
    
    # Plot # pruned SNPs vs. length
    par(mar = c(5, 4, 4, 8) + 0.1, bg = "white")
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
    dev.copy(png, paste(results.path, "NoOfSNPsAndLengthPlot.png", sep = "/"))
    dev.off()

# pruned SNPs vs amount of LD at 500K

    par(mar = c(5, 4, 4, 8) + 0.1, bg = "white")
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
    dev.copy(png, paste(results.path, "NoOfSNPsVsLD.png", sep = "/"))
    dev.off()

    # amount of LD at 500K vs length
    par(mar = c(5, 4, 4, 8) + 0.1, bg = "white")
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
    dev.copy(png, paste(results.path, "LDVsChrLength.png", sep = "/"))
    dev.off()

    # amount of LD at 500K vs exons/length (see thread on oyster genome slack channel)

    par(mar = c(5, 4, 4, 8) + 0.1, bg = "white")
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
    dev.copy(png, paste(results.path, "LDVsExonsPerChrLength.png", sep = "/"))
    dev.off()

    ###########################################

    par(mar = c(5, 4, 4, 8) + 0.1, bg = "white")
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
    dev.copy(png, paste(results.path, "SnpsPerOriginalSnpsVsLD.png", sep = "/"))
    dev.off()
}

deep.ld.analysis(ld.path = "data/large_data/ldAnalysisData/unrelated/",
                 thinned.snps = "data/thinned_snps/thinnedMatrixAndMetaData50000Window_unrelated.rds",
                 full.data = "data/large_data/genotypeMatrix_unrelated.rds",
                 subset = "unrelated")

deep.ld.analysis(ld.path = "data/large_data/ldAnalysisData/excluding_LM/",
                 thinned.snps = "data/thinned_snps/thinnedMatrixAndMetaData50000Window_exclude_LM.rds",
                 full.data = "data/large_data/genotypeMatrix_exclude_LM.rds",
                 subset = "exclude_LM")