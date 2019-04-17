###############################################################################
#
# File      : fstMatrixByChromosome.R 
# History   : 01/30/2019  Created by K Bodie Weedop (KBW)
#           : 02/12/2019  Option to select certain populations rather than just delete
#           
###############################################################################

###############################################################################
#
# This script will use a genotype matrix, separate the matrix by chromosome and
# calculate FST matrices for each chromosome using the OutFLANK package. You
# must run populationStructureScript.R prior to this script in order to get such
# a genotype matrix. Each of the FST matrices will be saved in a folder where
# this script is ran. Each of the FST matrices will be named appropriately given
# the chromosome it is associated with and the populations that have been
# removed (or if any populations have been removed at all).
#
###############################################################################

library(OutFLANK)
source("src/4outlier/OutFLANKv2.R")
library(pcadapt)
library(bigsnpr)
library(dplyr)

getFixedSites <- function(genotype.matrix = NULL) {
    # Convert to format (file-backed matrix; 'FBM') that Prive uses
    # in his walk-through.

    no.samples <- ncol(genotype.matrix)

    G_FBM <- add_code256(big_copy(t(genotype.matrix), type="raw"),
                         code = bigsnpr:::CODE_012)

    # Search all of the SNP sites where there is no variation among individuals (all have
    # a value of 0, 1, or 2)
    counts <- big_counts(G_FBM)
    fixed.sites <- NULL
    for (i in 1:3) {
        fixed.sites <- c(fixed.sites, which(counts[i, ] == no.samples))
    }
    return(fixed.sites)
}

preprocess.fstByChr <- function(data = NULL,
                                metadata = NULL,
                                chr.index = NULL, 
                                subset.name = NULL) {
    # Load the genotype matrix that has been provided
    data <- readRDS(data)

    # Load metadata to check the population identifiers agaisnt that the user has provided.
    metadata <- read.csv(metadata)

    # Selecting just one chromosome to decouple the data processing
    one.chr <- list(G = data$G[which(data$Chr == chr.index),],
                    Pos = data$Pos[which(data$Chr == chr.index)],
                    Chr = data$Chr[which(data$Chr == chr.index)],
                    Pop.ID = data$Pop.ID,
                    Sample.ID = data$Sample.ID)

    # If the user has provided a subset to use, select it and remove
    # the column associated with it from the genotype matrix
    if (!is.null(subset.name)) {
        # Check if the correct data format is being used
        if (!subset.name %in% names(metadata)) {
            stop("Error: the subset that you have provided is not specified in the provided metadata.")
        }
        # Get Sample.ID from 'metadata' that need to be removed
        rm.samples <- metadata$Sample.ID[which(metadata[, subset.name] == 0)]
        #Check to see if any of the samples in the subset are in the data.
        if (!any(one.chr$Sample.ID %in% rm.samples)) {
            # If there are no samples specified by the subset in the data, stop.

        } else {
            one.chr$G <- one.chr$G[, -which(one.chr$Sample.ID %in% rm.samples)]
            one.chr$Pop.ID <- one.chr$Pop.ID[-which(one.chr$Sample.ID %in% rm.samples)]
            one.chr$Sample.ID <- one.chr$Sample.ID[-which(one.chr$Sample.ID %in% rm.samples)]
        }
    }
    
    # Check for any sites that are fixed across all individuals left in the data.
    fixed.sites <- getFixedSites(one.chr$G)
    # Remove any fixed sites found
    if (length(fixed.sites) != 0) {
        one.chr$G <- one.chr$G[-fixed.sites, ]
        one.chr$Pos <- one.chr$Pos[-fixed.sites]
        one.chr$Chr <- one.chr$Chr[-fixed.sites]
    }
    # Throw an error if the number of samples does not correspond to the number
    # of columns provided.
    if (length(one.chr$Pop.ID) != ncol(one.chr$G)) {
        stop("Number of Pop.ID in dataset and labels does not match")
    }
    
    if (length(one.chr$Sample.ID) != ncol(one.chr$G)) {
        stop("Number of Sample.ID in dataset and labels does not match")
    }

    # Return the processed data
    return(one.chr)
}

fstByChromosome <- function(data.file,
                            metadata,
                            subset.name = NULL,
                            chr.index = seq(1, 10),
                            data.path = NULL, 
                            plots.path = NULL) {
    # Create path for output
    data.path <- file.path(data.path, "fstByChromosome")
    if (!dir.exists(data.path)) {
        dir.create(data.path)
    }
    # Run the preprocess funtion and calculate fst matrix for each chromosome 
    for (i in chr.index) {
        preprocessedData <- preprocess.fstByChr(data.file,
                                                metadata,
                                                subset.name = subset.name,
                                                chr.index = i)

        fst.stat <- MakeDiploidFSTMat(t(preprocessedData$G), 
                                      locusNames = preprocessedData$Pos, 
                                      popNames = preprocessedData$Pop.ID)

        fst.stat$Pos <- fst.stat$LocusName

        fst.stat$Chr <- i 

        saveRDS(fst.stat, paste(data.path, 
                                "/fstByChromosome_chr", i, 
                                "_", gsub("\\.", "_", subset.name), 
                                ".rds",
                                sep=""))
    }
}


###############################################################################
#
# File      : fstMatrix_subset.R 
# History   : 01/30/2019  Created by K Bodie Weedop (KBW)
#           
###############################################################################

###############################################################################
#
# This script will use OutFLANK with set of random quasi-independent SNPs in
# order to establish a neutral parameterization for outlier detection. The user
# should already have the subset of SNPs in order to run this script. 
#
###############################################################################

preprocess.neutralParams <- function(rand.snps,
                                     metadata,
                                     subset.name = NULL) {
    # Load data with subset of random independent SNPs
    data <- readRDS(rand.snps)
    # Load metadata to check the population IDs that are given by the user.
    metadata <- read.csv(metadata)

    # Need to remove any populations that should be excluded from the neutral
    # paramaterization.
    if (!is.null(subset.name)) {
        if (!subset.name %in% names(metadata)){
            stop("Error: the subset that you have provided is not specified in the metadata")
        }
        # Get Sample.ID from 'metadata' that need to be removed
        rm.samples <- metadata$Sample.ID[which(metadata[, subset.name] == 0)]
        #Check to see if any of the samples in the subset are in the data.
        if (!any(data$Sample.ID %in% rm.samples)) {
            # If there are no samples specified by the subset in the data        
        } else {
            data$G <- data$G[, -which(data$Sample.ID %in% rm.samples)]
            data$Pop.ID <- data$Pop.ID[-which(data$Sample.ID %in% rm.samples)]
            data$Sample.ID <- data$Sample.ID[-which(data$Sample.ID %in% rm.samples)]
        }
    }

    fixed.sites <- getFixedSites(data$G)
    if (length(fixed.sites) != 0) {
        data$G <- data$G[-fixed.sites, ]
        data$Pos <- data$Pos[-fixed.sites]
        data$Chr <- data$Chr[-fixed.sites]
    }    

    if (length(data$Pop.ID) != ncol(data$G)) {
        stop("Number of Pop.ID in dataset and labels does not match")
    }
    
    if (length(data$Sample.ID) != ncol(data$G)) {
        stop("Number of Sample.ID in dataset and labels does not match")
    }

    return(data)
}

altered.outflank <- function(data, NumberOfSamples) {
    # keep <- apply(fst.stat, 1, function(x) all(!is.na(x)))
    # keep2 <- complete.cases(fst.stat)
    # 
    #fst.stat <- fst.stat[keep, ]
    ## fst.stat <- fst.stat[which(fst.stat$meanAlleleFreq >= 0.05), ]
    # zeros <- which(rowSums(fst.stat[, 2:ncol(fst.stat)]) == 0)
    # zeros <- as.numeric(attr(zeros, "names"))

    outflank.data <- OutFLANK.v2(data, 
                                 NumberOfSamples = NumberOfSamples, 
                                 qthreshold = 0.05, 
                                 Hmin = 0.1)
    return(outflank.data)
}

neutralParams <- function (rand.snps,
                           metadata,
                           subset.name = NULL,
                           data.path = NULL, 
                           plots.path = NULL) {
    if (length(data.path) == 0 | length(plots.path) == 0) {
        stop("Error: No file path given for data or plots")
    }

    data.path <- file.path(data.path, "neutralParams")
    if (!dir.exists(data.path)) {
        dir.create(data.path)
    }

    plots.path <- file.path(plots.path, "neutralParams")
    if (!dir.exists(plots.path)) {
        dir.create(plots.path)
    }

    tmp <- strsplit(rand.snps, split="/")[[1]]
    rand.file <- strsplit(tmp[length(tmp)], split=".rds")[[1]]    
    
    subset.data <- preprocess.neutralParams(rand.snps = rand.snps,
                                            metadata = metadata, 
                                            subset.name = subset.name)


    pops <- length(unique(subset.data$Pop.ID))

    fst.stat <- MakeDiploidFSTMat(t(subset.data$G), 
                                  locusNames = subset.data$Pos, 
                                  popNames = subset.data$Pop.ID)

    outflank.data <- tryCatch(OutFLANK(fst.stat, 
                                       NumberOfSamples = pops, 
                                       qthreshold = 0.05, 
                                       Hmin = 0.1),
                              error = function(cond){
                                  print("Used altered.outflank")
                                  data <- altered.outflank(fst.stat, pops)
                                  return(data)
                              })

    plot(fst.stat$He, fst.stat$FST)
    dev.copy(png, paste(plots.path, 
                        "/heteroVsFst_",
                        gsub("\\.", "_", subset.name), 
                        ".png", sep=""))
    dev.off()

    plot(fst.stat$FST, fst.stat$FSTNoCorr)
    abline(0,1)
    dev.copy(png, paste(plots.path, 
                        "/fstNoCorrVsFst_",
                        gsub("\\.", "_", subset.name), 
                        ".png", sep=""))
    dev.off()
    
    OutFLANKResultsPlotter(outflank.data, 
                           withOutliers = TRUE, 
                           NoCorr = TRUE, 
                           Hmin = 0.1, 
                           binwidth = 0.001, 
                           Zoom = FALSE, 
                           RightZoomFraction = 0.05, 
                           titletext = NULL)
    
    dev.copy(png, paste(plots.path, 
                        "/outflankResults_",
                        gsub("\\.", "_", subset.name), 
                        ".png", sep=""))
    dev.off()

    OutFLANKResultsPlotter(outflank.data, 
                   withOutliers = TRUE, 
                   NoCorr = TRUE, 
                   Hmin = 0.1, 
                   binwidth = 0.001, 
                   Zoom = TRUE, 
                   RightZoomFraction = 0.15, 
                   titletext = NULL)

    dev.copy(png, paste(plots.path, 
                        "/outflankResultsRightTail_",
                        gsub("\\.", "_", subset.name), 
                        ".png", sep=""))
    dev.off()

    saveRDS(fst.stat, paste(data.path, 
                            "/fstMatrix_",
                            gsub("\\.", "_", subset.name), 
                            ".rds", sep=""))
    
    saveRDS(outflank.data, paste(data.path,
                                 "/outflank_",
                                 gsub("\\.", "_", subset.name), 
                                 ".rds", sep=""))
    
}

###############################################################################
#
# File      : outflankAnalysis.R 
# History   : 12/18/2018  Created by K Bodie Weedop (KBW)
#           
###############################################################################

###############################################################################
#
# This script is to be run after fstMatrixByChromosome.R and fstMatrix_subset.R
# as you will need the results from it to proceed with this analysis. This
# analysis will use outFLANK to implement a method developed in Whitlock and
# Lotterhos (2015). The method uses likelihood within the scope of a trimmed
# distribution of Fst values to infer the distribution of Fst for neutral
# markers. This distribution is then used to assign q-values to each locus to
# detect outliers that may be due to spatially heterogeneous selection.
#
###############################################################################

altered.pOutlierFinder <- function(data, dist) {
    return(pOutlierFinderChiSqNoCorr.v2(data, 
                                       Fstbar = dist$FSTNoCorrbar, 
                                       dfInferred = dist$dfInferred, 
                                       qthreshold = 0.00001, 
                                       Hmin = 0.1))
}

add.breakpoints <- function (chr.index) {
       
    breakpoints <- read.table("breakpoints_for_plots.txt")

    breakpoints <- breakpoints[which(breakpoints$Chr == chr.index), ]

    for (i in 1:nrow(breakpoints)) {
        rect(xleft=breakpoints$Low[i], 
             ybottom=-0.25, 
             xright=breakpoints$High[i],
             ytop=1.15,
             col="lightgrey",
             border="transparent")
        x_text <- floor((breakpoints$High[i]+breakpoints$Low[i])/2)
        text(x=x_text, y=1.1 , paste("LG", breakpoints$LG[i], sep=" "))
    }
}

outflank.outlierFinder <- function(data.file = NULL,
                                   null.dist = NULL,
                                   subset.name = NULL,
                                   data.path = NULL,
                                   plots.path = NULL) {
    # Creating folder for outlier plots. If the folder does not exist, create
    # it
    null.file <- strsplit(null.dist, split=".rds")[[1]]
    plots.path <- file.path(plots.path, 
                            "outflank.outlierFinder")
    if (!dir.exists(plots.path)) {
         dir.create(plots.path, recursive = TRUE)
    }
    # Creating folder for outlier data. If the folder does not exist, create
    # it
    results.path <- file.path(data.path, 
                              "outflank.outlierFinder", 
                              "results")
    if (!dir.exists(results.path)) {
         dir.create(results.path, recursive = TRUE)
    }

    # misassembled - misassembled chromosomes that need to have rectangles added.
    misassembled <- c(5,6,9)

    # Load FST matrix from folder
    fst.stat <- readRDS(paste(data.path, "fstByChromosome", data.file, sep="/"))
    # Get the name of the full FST matrix file that is being used.
    data.file <- strsplit(data.file, split=".rds")[[1]]
    # Get the random quasi-independent SNPs in order to be used to establish
    # neutral parameters
    outflank.data <- readRDS(paste(data.path, "neutralParams", null.dist, sep="/"))
    # Find outliers using the neutral parameters
    P1 <- tryCatch(pOutlierFinderChiSqNoCorr(fst.stat, 
                                             Fstbar = outflank.data$FSTNoCorr, 
                                             dfInferred = outflank.data$dfInferred, 
                                             qthreshold = 0.00001, 
                                             Hmin = 0.1),
                    error=function(cond){
                        P1 <- altered.pOutlierFinder(data = fst.stat, 
                                                     dist = outflank.data)
                    })
    
    outlier <- P1$OutlierFlag == TRUE

    # png(paste(folder.path, "/outlierPlotUsing_", data.file , ".png", sep=""), 
    #     height = 512,
    #     width = 2048)
    # plot(P1$He, P1$FST, pch = 19, col = rgb(0,0,0,0.1))
    # points(P1$He[outlier], P1$FST[outlier], col = "blue")
    # dev.off()

     # hist(P1$pvaluesRightTail)

    chr <- as.numeric(gsub("[^0-9]", "",  data.file))

    plot_x <- P1$LocusName[P1$He > 0.1]

    if (nchar(length(plot_x)) == 6) {
        tick_intervals <- round(length(plot_x) / 10000)
    } else {
        tick_intervals <- round(length(plot_x) / 5000)
    }
    
    ticks <- plot_x[seq(1, length(plot_x), by = round(length(plot_x)/tick_intervals))]

    png(paste(plots.path, 
              "/manhattanPlotUsing_", 
              data.file, 
              ".png", sep=""), 
        height = 720,
        width = 2048)
    par(mar = c(8, 4, 2, 1))
    plot(P1$LocusName[P1$He > 0.1], 
         P1$FST[P1$He > 0.1],
         ylab="FST",
         xlab = "",
         ylim=c(-0.2, 1.1),
         xlim=c(min(P1$LocusName), max(P1$LocusName)),
         main = paste("Manhattan Plot: Chromsome", chr), 
         col=rgb(0, 0, 0, 0.2),
         xaxt="n")
    axis(side = 1,
         las = 2,
         at=ticks)
    mtext(text = "Position (BP)",
          side = 1,
          line = 6)
    if (chr %in% misassembled) {
        add.breakpoints(chr.index = chr)
        points(P1$LocusName[P1$He>0.1], P1$FST[P1$He>0.1], pch=20)
    }
    points(P1$LocusName[outlier], P1$FST[outlier], col="red", pch=20)
    dev.off()

    tmp <- strsplit(data.path, split="/")[[1]][3]

    identifier <- strsplit(tmp, split="outlierAnalysis_")[[1]][2]

    pos <- P1$Pos
    chr <- P1$Chr

    P1 <- P1[,-which(names(P1) %in% c("Chr", "Pos"))]

    colnames(P1)[2:ncol(P1)] <- paste("OutFLANK_0.2_bigsnpr_0.8.2", 
                                      colnames(P1)[2:ncol(P1)],
                                      gsub("\\.", "_", subset.name),
                                      sep="_")

    P1$Pos <- pos
    P1$Chr <- chr

    # max.chr - largest number of digits in the chromosome index
    max.chr <- 2
    # max.pos - largest number of digits in the position index
    max.pos <- 9
    P1$unique <- paste(sprintf(paste("%0", max.chr, "d", sep=""), P1$Chr), 
                       sprintf(paste("%0", max.pos, "d", sep=""), P1$Pos), 
                       sep="_")
    P1 <- within(P1, rm("LocusName"))
    
    dir <- getwd()
    setwd(paste(dir, results.path, sep="/"))
    
    saveRDS(P1, paste("outlierFinder_", data.file, ".rds", sep=""))

    setwd(dir)
    return(P1)
}

###############################################################################
#
# File      : pcadaptAnalysis.R 
# History   : 02/02/2019  Created by K Bodie Weedop (KBW)
#           : 02/05/2019  Code altered to run all analyses with one function. 
#
###############################################################################

###############################################################################
#
# This script runs an analysis chromosome by chromosome using the pcadapt
# method. The results will be saved to an rds file in a folder named
# 'pcadaptResults.' 
#
# The default analysis will run on the full genotype matrix using the metadata
# and a random subset of quasi-independent SNPs associated with it. The user can
# specify a few different parameters that alter the analysis a bit. The first is
# the chromosome index (integer) that should be used during the analysis. The
# second is the   
#
###############################################################################

pcadaptAnalysis <- function (data,
                             metadata,
                             rand.snps,
                             subset.name = NULL,
                             chr.index = NULL,
                             data.path = NULL,
                             plots.path = NULL) {
    # Loading the full genotype matrix with associated bp positions, chromosomal
    # positions and population id
    data <- readRDS(data)
    # Load metadata to check the population identifiers agaisnt that the user has provided.
    metadata <- read.csv(metadata)
    # Load the subset of random independent SNPs 
    subset.data <- readRDS(rand.snps)

    # Need to remove any populations that should be excluded from data
    if (!is.null(subset.name)) {
        if (!subset.name %in% names(metadata)){
            stop("Error: the subset that you have provided is not specified in the metadata")
        }
        # Get Sample.ID from 'metadata' that need to be removed
        rm.samples <- metadata$Sample.ID[which(metadata[, subset.name] == 0)]
        #Check to see if any of the samples in the subset are in the data.
        if (!any(data$Sample.ID %in% rm.samples)) {
            # If there are no samples specified by the subset in the data        
        } else {
            # Perform operations for the full dataset
            data$G <- data$G[, -which(data$Sample.ID %in% rm.samples)]
            data$Pop.ID <- data$Pop.ID[-which(data$Sample.ID %in% rm.samples)]
            data$Sample.ID <- data$Sample.ID[-which(data$Sample.ID %in% rm.samples)]
            # And again for the thinned SNP subset.
            subset.data$G <- subset.data$G[, -which(subset.data$Sample.ID %in% rm.samples)]
            subset.data$Pop.ID <- subset.data$Pop.ID[-which(subset.data$Sample.ID %in% rm.samples)]
            subset.data$Sample.ID <- subset.data$Sample.ID[-which(subset.data$Sample.ID %in% rm.samples)]
        }
    }

    # After removing some individuals, you may have fixed sites if the removed
    # individuals were providing the only unique allele type at that site.
    # We need to find these and remove them from the matrix, Pos, and Chr
    fixed.sites <- getFixedSites(data$G)
    if (length(fixed.sites) != 0) {
        data$G <- data$G[-fixed.sites, ]
        data$Pos <- data$Pos[-fixed.sites]
        data$Chr <- data$Chr[-fixed.sites]
    }
    # We will do the same for the thinned SNPs 
    subset.fixed.sites <- getFixedSites(subset.data$G)
    if (length(subset.fixed.sites) != 0) {
        subset.data$G <- subset.data$G[-subset.fixed.sites, ]
        subset.data$Pos <- subset.data$Pos[-subset.fixed.sites]
        subset.data$Chr <- subset.data$Chr[-subset.fixed.sites]
    }

    # Getting the coded genotype matrix of the subset in order to run snp_autoSVD()
    subset.coded <- add_code256(big_copy(t(subset.data$G),
                                         type = "raw"), 
                                code = bigsnpr:::CODE_012)
    subset.data$subset.coded <- subset.coded
    
    # Removing the variables which will not be used any further
    rm(subset.coded)

    newpc <- big_randomSVD(X = subset.data$subset.coded, 
                           fun.scaling = snp_scaleBinom(), 
                           k = 10, 
                           ncores = 4)
    # newpc <- snp_autoSVD(G = subset.data$subset.coded,
    #                      infos.chr = subset.data$Chr,
    #                      infos.pos = subset.data$Pos,
    #                      size = 50000,
    #                      is.size.in.bp = TRUE)        

    if (length(data$Pop.ID) != ncol(data$G)) {
        stop("Number of Pop.ID in dataset and labels does not match")
    }
    
    if (length(data$Sample.ID) != ncol(data$G)) {
        stop("Number of Sample.ID in dataset and labels does not match")
    }
    
    combined.data <- data.frame()

    for (i in chr.index) {
       # Subsetting the full data by chromosome
        
        if (!i %in% unique(data$Chr)) {
            stop("Error: The chromosome index which has been specified is not an index found in the data")
        }

        # Select only the SNPs associated with the specified chromosome
        one.chr <- list(G = data$G[which(data$Chr == i),],
                        Pos = data$Pos[which(data$Chr == i)],
                        Chr = data$Chr[which(data$Chr == i)],
                        Pop.ID = data$Pop.ID,
                        Sample.ID = data$Sample.ID)
        
        # keep <- apply(allData$G, 1, function(x) length(unique(x[!is.na(x)])) != 1)
        # data$genotype <- data$genotype[keep, ]
        # data$Positions <- data$Positions[keep]
        # data$chromosome <- data$chromosome[keep]
        # 
        # subset.keep <- apply(subset.data$G, 1, function(x) length(unique(x[!is.na(x)])) != 1)
        # subset.data$G <- subset.data$G[subset.keep, ]
        # subset.data$Positions <- subset.data$Positions[subset.keep]
        # subset.data$chromosomes <- subset.data$chromosome[subset.keep]

        # Getting the coded genotype matrix in order to run snp_autoSVD()
        G.coded <- add_code256(big_copy(t(one.chr$G),
                                        type="raw"), 
                               code=bigsnpr:::CODE_012)

        # appending the coded genotype matrix to the original data
        one.chr$G.coded <- G.coded

        max.pos <- 9
        # Max char for the chromosome will have to be hardcoded here.
        # This is due to the fact that we are 
        max.chr <- 2 
        unique <- paste(sprintf(paste("%0", max.chr, "d", sep=""), one.chr$Chr), 
                        sprintf(paste("%0", max.pos, "d", sep=""), one.chr$Pos), 
                        sep="_")

        # pcadapt on all loci but with the pruned PCs
        pcadapt_4.0.3_bigsnpr_0.8.2_Md <- snp_pcadapt(one.chr$G.coded, U.row = newpc$u[, 1:5])
        full_stats_pcadapt_4.0.3_bigsnpr_0.8.2 <- snp_gc(pcadapt_4.0.3_bigsnpr_0.8.2_Md)
        # get the negative log10 values of the results.
        pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p <- -predict(full_stats_pcadapt_4.0.3_bigsnpr_0.8.2, log10=T)
        # put the full results and the negative log10 values into a dataframe.

        # I could put the entire dataset into the file but not needed at the time of writing.
        # "full_stats_pcadapt_4.0.3_bigsnpr_0.8.2" = full_stats_pcadapt_4.0.3_bigsnpr_0.8.2, 

        pcadapt.results <- data.frame("Pos" = one.chr$Pos,
                                      "Chr" = one.chr$Chr,
                                      "unique" = unique,
                                      "pcadapt_4.0.3_bigsnpr_0.8.2_Md" = pcadapt_4.0.3_bigsnpr_0.8.2_Md$score,
                                      "pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p" = pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p)

        names(pcadapt.results)[4] <- paste("pcadapt_4.0.3_bigsnpr_0.8.2_Md", 
                                            gsub("\\.", "_", subset.name), 
                                            sep="_")
        
        names(pcadapt.results)[5] <- paste("pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p", 
                                            gsub("\\.", "_", subset.name), 
                                            sep="_")

        results.path <- file.path(data.path, 
                                  paste("/pcadaptAnalysis_", 
                                        gsub("\\.", "_", subset.name), 
                                        sep=""))
        if (!dir.exists(results.path)) {
            dir.create(results.path)
        }
        
        saveRDS(as.data.frame(pcadapt.results), 
                paste(results.path, 
                      "/pcadaptResults_chr", 
                      i, 
                      ".rds",
                      sep = ""))
        
        combined.data <- rbind(combined.data, pcadapt.results)
    }
    combined.data$unique <- levels(droplevels(combined.data$unique))
    return(combined.data)
}

merge.results <- function (outflank.data, pcadapt.data) {
    joined.data <- full_join(out.data, pc.data, by = c("unique", "Pos", "Chr"))
}
# Below is a wrapper function that will allow you to run all of the analyses
# that you would like while just pressing 'enter'. If you would like to run the
# analyses for all populations, leave the rm.pops and select.pops parameters as
# NULL. If you would like to run analyses where populations are removed or only
# certain populations are selected, you can concatenate those parameters into
# those parameters, respectively. However, do not put NULL into your
# concatenated parameters (i.e rm.pops = list(c("LM"), NULL, c("CL")) 

wrapper <- function (data.file = "data/large_data/genotypeMatrix_unrelated.rds",
                     metadata = "data/modified_samplemetadata.csv",
                     rand.snps = "data/thinned_snps/thinnedMatrixAndMetaData50000Window_unrelated.rds",
                     subset.name = subset.name,
                     chr.index = seq(1, 10)) {
    # Create folder structure for the analyses and their outputs depending upon
    # whether they remove or select certain populations
    
    data.path <- file.path("data",
                            "large_outputs", 
                            "outlierAnalysis", 
                            paste("outlierAnalysis_", gsub("\\.", "_", subset.name), sep=""),
                            "data")
    plots.path <- file.path("data",
                            "large_outputs", 
                             "outlierAnalysis", 
                             paste("outlierAnalysis_", gsub("\\.", "_", subset.name), sep=""),
                             "plots")
    
    if (!dir.exists(data.path)) {
        dir.create(data.path, recursive = TRUE)
    }
    if (!dir.exists(plots.path)) {
        dir.create(plots.path, recursive = TRUE)
    }

    fstByChromosome(data.file = data.file,
                    metadata = metadata, 
                    chr.index = chr.index,  
                    subset.name = subset.name,
                    data.path = data.path, 
                    plots.path = plots.path)
    
    neutralParams(rand.snps = rand.snps,
                  metadata = metadata, 
                  subset.name = subset.name,
                  data.path = data.path, 
                  plots.path = plots.path)
    
    # Check to see if the folder which has been provided has been created.
    data.folder <- file.path(data.path, "fstByChromosome")
    if (!dir.exists(data.folder)) {
         stop("The folder holding the data files does not exist")
    # Check to see if there are any files in the folder provided. If there are
    # none, return with an error     
    } 
        
    if (length(list.files(path = data.folder)) == 0) {
         stop("There are no data files in the folder that you have provided")
    }

    # Run the function above across the files which the user has provided.
    files <- list.files(path = data.folder)

    rand.file <- strsplit(rand.snps, split=".rds")[[1]]    
    null.folder <- file.path(data.path, "neutralParams")
    null.dist <- list.files(null.folder)[which(grepl("outflank", list.files(null.folder)))]
    
    if (length(null.dist) > 1) {
        null.dist <- null.dist[which(grepl(rand.file, null.dist))]
    }

    outflank.data <- data.frame()
    for (i in files) {
        tmp.outflank <- outflank.outlierFinder(data.file = i, 
                                               null.dist = null.dist,
                                               subset.name = subset.name,
                                               data.path = data.path,
                                               plots.path = plots.path)
        outflank.data <- rbind(outflank.data, tmp.outflank)
    }

    pcadapt.data <- pcadaptAnalysis(data = data.file, 
                                    chr.index = chr.index, 
                                    rand.snps = rand.snps,
                                    metadata = metadata, 
                                    subset.name = subset.name,
                                    data.path = data.path, 
                                    plots.path = plots.path)

    merged.data <- full_join(outflank.data, 
                             pcadapt.data, 
                             by = c("Pos", "Chr", "unique"))

    write.table(merged.data, paste(data.path, "/final_merged_data.txt", sep = ""))
}

 wrapper(data.file = "data/large_data/genotypeMatrix_exclude_LM.rds",
         rand.snps = "data/thinned_snps/thinnedMatrixAndMetaData50000Window_exclude_LM.rds",
         subset.name = "exclude.LM")

# wrapper(data.file = "data/large_data/genotypeMatrix_unrelated.rds",
#         rand.snps = "data/thinned_snps/thinnedMatrixAndMetaData50000Window_unrelated.rds",
#         subset.name = "unrelated")
