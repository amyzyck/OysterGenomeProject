###############################################################################
#
# File      : LDAnalysisScript.R 
# History   : 10/21/2018  Created by K Bodie Weedop (KBW)
#
###############################################################################
#
# This script works with the various LD analysis files given by VCFtools when
# running --geno-r2 on the original VCF file
# (SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz). 
#
###############################################################################


ld.data.50bp <- read.csv("ldAnalysisFiles\\geno_ld_window_50-50.txt", 
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
                "450" = ld.data.450.500bp, 
                "4500" = ld.data.4500.5000bp,
                "9500" = ld.data.9500.10kbp,
                "49500" = ld.data.49500.50kbp,
                "99500" = ld.data.99500.100kbp,
                "499500" = ld.data.499500.500kbp)

# get the number of loci in the data set associated with each of the chromosomes
num.loci <- function(ld.dataset) {
    chr.num <- length(unique(ld.dataset$CHR))
    for (i in 1:chr.num){
        print(paste(unique(ld.dataset$CHR)[i], 
        "Number of loci in LD calculations:",
        length(ld.dataset$CHR[which(ld.dataset$CHR == unique(ld.dataset$CHR)[i])]), sep=" "))
    }
}

for (i in 1:length(ld.vars)) {
    num.loci(ld.vars[[i]])
    print("--------------")
}

# get the LD decay plots for all of the chromosomes in the ld.datasets

decay.plots <- function(ld.datasets) {
    CHR <- c("NC_035789.1", "NC_035780.1", "NC_035781.1", "NC_035782.1", 
             "NC_035783.1", "NC_035784.1", "NC_035785.1", "NC_035786.1", 
             "NC_035787.1", "NC_035788.1")
    col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                  '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                  '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
                  '#ffffff', '#000000')
    chr.num <- length(unique(ld.datasets[[1]]$CHR))
    for (i in length(ld.datasets)) {
        if (length(unique(ld.datasets[[i]]$CHR)) != chr.num) {
            stop("The number of chromosomes in the datasets do not match")
        } else {

        }
    }
    
    chr.index <- rep(NA, chr.num)
    for (i in 1:chr.num){
        chr.index[i] <- substr(unique(ld.datasets[[1]]$CHR)[i], start = 9, stop = 9)
    }

    ld.stats <- data.frame(matrix(NA, length(ld.datasets), chr.num))
    ld.stats[,1] <- as.numeric(names(ld.datasets))
    for (i in 1:length(ld.datasets)) {
        ld.datasets[[i]]$R.2 <- as.numeric(ld.datasets[[i]]$R.2)
        for (j in 1:chr.num) {
            ld.stats[i, j+1] <- summary(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])[4]
        }
    }

    par(bty="l")
    for (i in 2:ncol(ld.stats)) {
        if (i == 2){
            plot(jitter(ld.stats[,i]) ~ ld.stats[, 1],
                 pch = i,
                 col = col_vector[i],
                 cex = 2,
                 xlab="Distance (bp)", 
                 ylab="LD (r^2)",
                 ylim=c(0,0.24),
                 type = "b",
                 log = "x",
                 main="LD Decay")
        } else {
            points(jitter(ld.stats[,i]) ~ ld.stats[, 1],
                   pch = i,
                   col = col_vector[i],
                   cex = 2,
                   type = "b")
        }
        legend("topright", legend = CHR, pch = 2:(length(CHR)+1), col = col_vector[2:(length(CHR)+1)])
    }
    dev.copy(png, "ldDecayPlot.png")
    dev.off()
}
