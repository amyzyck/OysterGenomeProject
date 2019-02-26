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
source("OutFLANKv2.R")
library(pcadapt)
library(bigsnpr)

preprocess.fstByChr <- function(data,
                                metadata,
                                chr.index, 
                                rm.pops = NULL, 
                                select.pops = NULL) {
    # Load the genotype matrix that has been provided
    data <- readRDS(data)

    # Load metadata to check the population identifiers agaisnt that the user has provided.
    metadata <- read.csv(metadata)

    # Check if the user has provided pop ids that for exclusion and inclusion at the same time
    if (length(rm.pops) != 0 && length(select.pops) != 0) {
        stop("Error: Only one of the parameters 'rm.pops' and 'select.pops' should be provided by the user.")
    }

    # Selecting just one chromosome to decouple the data processing
    one.chr <- list(genotype = data$genotype[which(data$chromosome == chr.index),],
                    position = data$positions[which(data$chromosome == chr.index)],
                    chromosome = data$chromosome[which(data$chromosome == chr.index)],
                    sample.id = data$sample.id)
    # If the user has provided any populations to remove, select it and remove
    # the column associated with it from the genotype matrix
    if (length(rm.pops) != 0) {
        # Check if the correct data format is being used
        if (typeof(rm.pops) != "character") {
            stop("The parameter 'rm.pops' must be a character vector")
        }
        # Check if the populations specified are in the metadata
        if (!all(rm.pops %in% metadata$Pop.ID)) {
            stop("Error: at least one of the populations you have selected to remove is not in the data")
        }
        one.chr$genotype <- one.chr$genotype[, -which(one.chr$sample.id %in% rm.pops)]
        one.chr$sample.id <- data$sample.id[-which(one.chr$sample.id %in% rm.pops)]
        
    }

    # If populations have been selected to include only...
    if (length(select.pops) != 0) {
        # Check if the correct data format is being used
        if (typeof(select.pops) != "character") {
            stop("The parameter 'select.pops' must be a character vector")
        }
        # Check if the populations specified are in the metadata
        if (!all(select.pops %in% metadata$Pop.ID)) {
            stop("Error: at least one of the populations you have selected to keep is not in the data")
        }
        # Removing the specified populations from the full data
        one.chr$genotype <- one.chr$genotype[, which(one.chr$sample.id %in% select.pops)]
        one.chr$sample.id <- one.chr$sample.id[which(one.chr$sample.id %in% select.pops)]
    }
    # Throw an error if the number of samples does not correspond to the number
    # of columns provided.
    if (length(one.chr$sample.id) != ncol(one.chr$genotype)) {
        stop("Number of individuals in dataset and labels does not match")
    }
    # Return the processed data
    return(one.chr)
}

fstByChromosome <- function(data.file, 
                            chr.index = seq(1, 10),  
                            rm.pops = NULL, 
                            select.pops = NULL,
                            data.path = NULL, 
                            plots.path = NULL) {
    data.path <- file.path(data.path, "fstByChromosome")
    if (!dir.exists(data.path)) {
        dir.create(data.path)
    }
    
    for (i in chr.index) {
        preprocessedData <- preprocess.fstByChr(data.file, 
                                                chr.index = i, 
                                                rm.pops = rm.pops, 
                                                select.pops = select.pops)

        #keep <- apply(preprocessedData$genotype, 1, function(x) length(unique(x[!is.na(x)])) != 1)
        #preprocessedData$genotype <- preprocessedData$genotype[keep, ]
        #preprocessedData$positions <- preprocessedData$positions[keep]
        #preprocessedData$chromosome <- preprocessedData$chromosome[keep]

        fst.stat <- MakeDiploidFSTMat(t(preprocessedData$genotype), 
                                      locusNames = preprocessedData$position, 
                                      popNames = preprocessedData$sample.id)

        colnames(fst.stat)[1] <- "POS"

        fst.stat$CHR <- i 

        if (length(rm.pops) == 0 && length(select.pops) == 0) {
            colnames(fst.stat)[2:ncol(fst.stat)] <- paste("OutFLANK_0.2_PopSet_allpops", 
                                      colnames(fst.stat)[2:ncol(fst.stat)],
                                      sep="_")
            fst.stat$unique.id <- paste("chr_", i, "_pos_", fst.stat$POS, sep="")
            saveRDS(fst.stat, paste(data.path, 
                                    "/fstByChromosome_chr", i, 
                                    "_allpops.rds", sep=""))
        } else if (length(rm.pops) != 0) {
            colnames(fst.stat)[2:ncol(fst.stat)] <- paste("OutFLANK_0.2_PopSet",
                                      "excluding", paste(rm.pops, collapse="_"), 
                                      colnames(fst.stat)[2:ncol(fst.stat)],
                                      sep="_")
            fst.stat$unique.id <- paste("chr_", i, "_pos_", fst.stat$POS, sep="")
            saveRDS(fst.stat, paste(data.path, 
                                    "/fstByChromosome_chr", i, "_excludingpops_", 
                                    paste(rm.pops, collapse="_"), 
                                    ".rds", sep=""))
        } else {
            colnames(fst.stat)[2:ncol(fst.stat)] <- paste("OutFLANK_0.2_PopSet",
                                      "selecting", paste(select.pops, collapse="_"), 
                                      colnames(fst.stat)[2:ncol(fst.stat)],
                                      sep="_")
            fst.stat$unique.id <- paste("chr_", i, "_pos_", fst.stat$POS, sep="")
            saveRDS(fst.stat, paste(data.path, 
                                    "/fstByChromosome_chr", i, "_selectingpops_", 
                                    paste(select.pops, collapse="_"),
                                    ".rds", sep=""))
        }
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
                                     rm.pops = NULL,
                                     select.pops = NULL) {
    # Load data with subset of random independent SNPs
    data <- readRDS(rand.snps)
    # Load metadata to check the population IDs that are given by the user.
    metadata <- read.csv(metadata)
    # Check if the user has provided pop ids that for exclusion and inclusion at the same time
    if (length(rm.pops) != 0 && length(select.pops) != 0) {
        stop("Error: Only one of the parameters 'rm.pops' and 'select.pops' should be provided by the user.")
    }

    # Need to remove any populations that should be excluded from the neutral
    # paramaterization.
    if (length(rm.pops) != 0) {
        # Check if the correct data format is being used
        if (typeof(rm.pops) != "character") {
            stop("The parameter 'rm.pops' must be a character vector")
        }
        # Check if the populations specified are in the metadata
        if (!all(rm.pops %in% metadata$Pop.ID)) {
            stop("Error: at least one of the populations you have selected to remove is not in the data")
        }
        data$G <- data$G[, -which(data$sample.id %in% rm.pops)]
        data$sample.id <- data$sample.id[-which(data$sample.id %in% rm.pops)]
    }

    if (length(select.pops) != 0) {
        # Check if the correct data format is being used
        if (typeof(select.pops) != "character") {
            stop("The parameter 'rm.pops' must be a character vector")
        }
        # Check if the populations specified are in the metadata
        if (!all(select.pops %in% metadata$Pop.ID)) {
            stop("Error: at least one of the populations you have selected to remove is not in the data")
        }
        data$G <- data$G[, which(data$sample.id %in% select.pops)]
        data$sample.id <- data$sample.id[which(data$sample.id %in% select.pops)]
    }

    if (length(data$sample.id) != ncol(data$G)) {
        stop("Number of individuals in dataset and labels does not match")
    }

    return(data)
}

altered.outflank <- function(data, NumberOfSamples) {
    outflank.data <- OutFLANK.v2(data, 
                                 NumberOfSamples = NumberOfSamples, 
                                 qthreshold = 0.05, 
                                 Hmin = 0.1)
    return(outflank.data)
}

neutralParams <- function (rand.snps, 
                           rm.pops = NULL,
                           select.pops = NULL,
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
                                            rm.pops = rm.pops, 
                                            select.pops = select.pops)

    
    # subset.keep <- apply(subset.data$G, 1, function(x) length(unique(x[!is.na(x)])) != 1)
    # subset.data$G <- subset.data$G[subset.keep, ]
    # subset.data$positions <- subset.data$positions[subset.keep]
    # subset.data$chromosomes <- subset.data$chromosome[subset.keep]

    pops <- length(unique(subset.data$sample.id))

    fst.stat <- MakeDiploidFSTMat(t(subset.data$G), 
                                  locusNames = subset.data$positions, 
                                  popNames = subset.data$sample.id)

    outflank.data <- tryCatch(OutFLANK(fst.stat, 
                                       NumberOfSamples = pops, 
                                       qthreshold = 0.05, 
                                       Hmin = 0.1),
                              error=function(cond){
                                  print("Used altered.outflank")
                                  data <- altered.outflank(fst.stat, pops)
                                  return(data)
                              })

    # keep <- apply(fst.stat, 1, function(x) all(!is.na(x)))
    # keep2 <- complete.cases(fst.stat)
    # 
    #fst.stat <- fst.stat[keep, ]
    ## fst.stat <- fst.stat[which(fst.stat$meanAlleleFreq >= 0.05), ]
    # zeros <- which(rowSums(fst.stat[, 2:ncol(fst.stat)]) == 0)
    # zeros <- as.numeric(attr(zeros, "names"))

    if (length(rm.pops) != 0) {
        plot(fst.stat$He, fst.stat$FST)
        dev.copy(png, paste(plots.path, 
                            "/heteroVsFst_using_",
                            rand.file,
                            "_excluding_", 
                            paste(rm.pops, collapse="_"), 
                            ".png", sep=""))
        dev.off()

        plot(fst.stat$FST, fst.stat$FSTNoCorr)
        abline(0,1)
        dev.copy(png, paste(plots.path, 
                            "/fstNoCorrVsFst_using_",
                            rand.file,
                            "_excluding_", 
                            paste(rm.pops, collapse="_"), 
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
                            "/outflankResults_using_",
                            rand.file,
                            "_excluding_", 
                            paste(rm.pops, collapse="_"), 
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
                            "/outflankResultsRightTail_using",
                            rand.file,
                            "_excluding_", 
                            paste(rm.pops, collapse="_"), 
                            ".png", sep=""))
        dev.off()

        saveRDS(fst.stat, paste(data.path, 
                                "/fstMatrix_using",
                                rand.file,
                                "_excluding_", 
                                paste(rm.pops, collapse="_"), 
                                ".rds", sep=""))
        saveRDS(outflank.data, paste(data.path,
                                     "/outflank_using",
                                    rand.file,
                                    "_excluding_", 
                                     paste(rm.pops, collapse="_"), 
                                     ".rds", sep=""))
    } else if (length(select.pops) != 0) {
        plot(fst.stat$He, fst.stat$FST)
        dev.copy(png, paste(plots.path, 
                            "/heteroVsFst_using",
                            rand.file,
                            "_selecting_",
                            paste(select.pops, collapse="_"), 
                            ".png", sep=""))
        dev.off()

        plot(fst.stat$FST, fst.stat$FSTNoCorr)
        abline(0,1)
        dev.copy(png, paste(plots.path, 
                            "/fstNoCorrVsFst_using",
                            rand.file,
                            "_selecting_", 
                            paste(select.pops, collapse="_"), 
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
                            "/outflankResults_using",
                            rand.file,
                            "_selecting_",
                            paste(select.pops, collapse="_"), 
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
                            "/outflankResultsRightTail_using",
                            rand.file,
                            "_selecting_", 
                            paste(select.pops, collapse="_"), 
                            ".png", sep=""))
        dev.off()

        saveRDS(fst.stat, paste(data.path,
                                "/fstMatrix_using",
                                rand.file,
                                "_selecting_", 
                                paste(select.pops, collapse="_"), 
                                ".rds", sep=""))
        saveRDS(outflank.data, paste(data.path, 
                                     "/outflank_using",
                                    rand.file,
                                    "_selecting_", 
                                     paste(select.pops, collapse="_"), 
                                     ".rds", sep=""))
    } else {
        plot(fst.stat$He, fst.stat$FST)
        dev.copy(png, paste(plots.path, "heteroVsFst_allpops.png", sep="/"))
        dev.off()

        plot(fst.stat$FST, fst.stat$FSTNoCorr)
        abline(0,1)
        dev.copy(png, paste(plots.path, "fstNoCorrVsFst_allpops.png", sep="/"))
        dev.off()

        OutFLANKResultsPlotter(outflank.data, 
                               withOutliers = TRUE, 
                               NoCorr = TRUE, 
                               Hmin = 0.1, 
                               binwidth = 0.001, 
                               Zoom = FALSE, 
                               RightZoomFraction = 0.05, 
                               titletext = NULL)
        dev.copy(png, paste(plots.path, "outflankResults_allpops.png", sep="/"))
        dev.off()

        OutFLANKResultsPlotter(outflank.data, 
                       withOutliers = TRUE, 
                       NoCorr = TRUE, 
                       Hmin = 0.1, 
                       binwidth = 0.001, 
                       Zoom = TRUE, 
                       RightZoomFraction = 0.15, 
                       titletext = NULL)
        dev.copy(png, paste(plots.path, "outflankResultsRightTail_allpops.png", sep="/"))
        dev.off()

        saveRDS(fst.stat, paste(data.path,
                                paste("fstMatrix_using",
                                      rand.file,
                                      "allpops.rds",
                                      sep="_"),
                                sep="/"))
        saveRDS(outflank.data, paste(data.path,
                                     paste("outflank_using",
                                            rand.file,
                                            "allpops.rds",
                                            sep="_"),
                                sep="/"))
    }
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

outflank.outlierFinder <- function(data.file = NULL,
                                   null.dist = NULL,
                                   data.path = NULL,
                                   plots.path = NULL) {
    # Creating folder for outlier plots. If the folder does not exist, create
    # it
    null.file <- strsplit(null.dist, split=".rds")[[1]]
    plots.path <- file.path(plots.path, 
                            "outflank.outlierFinder", 
                            paste("neutralParams", null.file, sep="_"))
    if (!dir.exists(plots.path)) {
         dir.create(plots.path, recursive = TRUE)
    }
    # Creating folder for outlier data. If the folder does not exist, create
    # it
    results.path <- file.path(data.path, 
                              "outflank.outlierFinder", 
                              paste("neutralParams", null.file, sep="_"))
    if (!dir.exists(results.path)) {
         dir.create(results.path, recursive = TRUE)
    }

    # Load FST matrix from folder
    fst.stat <- readRDS(paste(data.path, "fstByChromosome", data.file, sep="/"))
    # Get the name of the full FST matrix file that is being used.
    data.file <- strsplit(data.file, split=".rds")[[1]]
    # Get the random quasi-independent SNPs in order to be used to establish
    # neutral parameters
    outflank.data <- readRDS(paste(data.path, "neutralParams", null.dist, sep="/"))
    # Find outliers using the neutral parameters
    P1 <- tryCatch(pOutlierFinderChiSqNoCorr(fst.stat, 
                                             Fstbar = outflank.data$FSTNoCorrbar, 
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
     
    png(paste(plots.path, 
              "/manhattanPlotUsing_", 
              data.file, 
              ".png", sep=""), 
        height = 720,
        width = 2048)
    par(mar = c(8,4,2,1))
    plot(P1$LocusName[P1$He>0.1], 
         P1$FST[P1$He>0.1],
         ylab="FST",
         xlab = "",
         ylim = c(-0.2, 1.0),
         main = paste("Manhattan Plot: Chromsome", chr), 
         col=rgb(0,0,0,0.2),
         xaxt="n")
    axis(side = 1,
         las = 2,
         at=seq(P1$LocusName[1], 
         P1$LocusName[length(P1$LocusName)], 
         by=1000000))
    mtext(text = "Position (BP)",
          side = 1,
          line = 6)
    points(P1$LocusName[outlier], P1$FST[outlier], col="magenta", pch=20)
    dev.off()

    tmp <- strsplit(data.path, split="/")[[1]][3]

    identifier <- strsplit(tmp, split="outlierAnalysis_")[[1]][2]

    colnames(P1)[1] <- "POS"
    P1$CHR <- chr

    colnames(P1)[2:ncol(P1)] <- paste("OutFLANK_0.2", 
                                      identifier, 
                                      colnames(P1)[2:ncol(P1)],
                                      sep="_")
    P1$unique.id <- paste("chr_", P1$CHR, "_pos_", P1$POS, sep="")

    saveRDS(P1, paste(results.path, "/outflank_outlierFinder_", data.file, ".rds", sep=""))
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

pcadaptAnalysis <- function (genotype.matrix,
                             metadata,
                             rand.snps,
                             chr.index = NULL,
                             rm.pops = NULL,
                             select.pops = NULL,
                             data.path = NULL,
                             plots.path = NULL) {
    # Loading the full genotype matrix with associated bp positions, chromosomal
    # positions and population id
    data <- readRDS(genotype.matrix)
    # Load metadata to check the population identifiers agaisnt that the user has provided.
    metadata <- read.csv(metadata)
    # Load the subset of random independent SNPs 
    subset.data <- readRDS(rand.snps)

    # Check if the user has provided pop ids that for exclusion and inclusion at the same time
    if (length(rm.pops) != 0 && length(select.pops) != 0) {
        stop("Error: Only one of the parameters 'rm.pops' and 'select.pops' should be provided by the user.")
    }
    
    # If populations have been selected to be removed...
    if (length(rm.pops) != 0) {
        # Check if the correct data format is being used
        if (typeof(rm.pops) != "character") {
            stop("The parameter 'rm.pops' must be a character vector")
        }
        # Check if the populations specified are in the metadata
        if (!all(rm.pops %in% metadata$Pop.ID)) {
            stop("Error: at least one of the populations you have selected to remove is not in the data")
        }
        # Removing the specified populations from the full data
        data$genotype <- data$genotype[, -which(data$sample.id %in% rm.pops)]
        data$sample.id <- data$sample.id[-which(data$sample.id %in% rm.pops)]
        # Removing the specified populations from the subset
        subset.data$G <- subset.data$G[, -which(subset.data$sample.id %in% rm.pops)]
        subset.data$sample.id <- subset.data$sample.id[-which(subset.data$sample.id %in% rm.pops)]
    }
    
    # If populations have been selected to include only...
    if (length(select.pops) != 0) {
        # Check if the correct data format is being used
        if (typeof(select.pops) != "character") {
            stop("The parameter 'select.pops' must be a character vector")
        }
        # Check if the populations specified are in the metadata
        if (!all(select.pops %in% metadata$Pop.ID)) {
            stop("Error: at least one of the populations you have selected to keep is not in the data")
        }
        # Removing the specified populations from the full data
        data$genotype <- data$genotype[, which(data$sample.id %in% select.pops)]
        data$sample.id <- data$sample.id[which(data$sample.id %in% select.pops)]
        # Removing the specified populations from the subset
        subset.data$G <- subset.data$G[, which(subset.data$sample.id %in% select.pops)]
        subset.data$sample.id <- subset.data$sample.id[which(subset.data$sample.id %in% select.pops)]
    }

    # Subsetting the full data by chromosome
    if (typeof(chr.index) != "NULL") {
        # Check to see if the chr.index specified is in the data
        if (!chr.index %in% unique(data$chromosome)) {
            stop("Error: The chromosome index which has been specified is not an index found in the data")
        }
        # Select only the SNPs associated with the specified chromosome
        data$genotype <- data$genotype[which(data$chromosome == chr.index), ]
        data$positions <- data$positions[which(data$chromosome == chr.index)]
        data$chromosome <- data$chromosome[which(data$chromosome == chr.index)]
    }

    # keep <- apply(data$genotype, 1, function(x) length(unique(x[!is.na(x)])) != 1)
    # data$genotype <- data$genotype[keep, ]
    # data$positions <- data$positions[keep]
    # data$chromosome <- data$chromosome[keep]
    # 
    # subset.keep <- apply(subset.data$G, 1, function(x) length(unique(x[!is.na(x)])) != 1)
    # subset.data$G <- subset.data$G[subset.keep, ]
    # subset.data$positions <- subset.data$positions[subset.keep]
    # subset.data$chromosomes <- subset.data$chromosome[subset.keep]

    # Getting the coded genotype matrix in order to run snp_autoSVD()
    G.coded <- add_code256(big_copy(t(data$genotype),
                                    type="raw"), 
                           code=bigsnpr:::CODE_012)
    # Getting the coded genotype matrix of the subset in order to run snp_autoSVD()
    subset.coded <- add_code256(big_copy(t(subset.data$G),
                                         type = "raw"), 
                                code = bigsnpr:::CODE_012)
    
    # appending the coded genotype matrix to the original data
    data$G.coded <- G.coded

    subset.data$subset.coded <- subset.coded
    
    # Removing the variables which will not be used any further
    rm(G.coded, subset.coded)

    newpc <- snp_autoSVD(G = subset.data$subset.coded,
                         infos.chr = subset.data$chromosomes,
                         infos.pos = subset.data$positions,
                         is.size.in.bp = TRUE)

    # pcadapt on all loci but with the pruned PCs
    full_stats_pcadapt_4.0.3_bigsnpr_0.8.2 <- snp_gc(snp_pcadapt(data$G.coded, U.row = newpc$u[, 1:10]))
    # get the negative log10 values of the results.
    pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p <- -predict(full_stats_pcadapt_4.0.3_bigsnpr_0.8.2, log10=T)
    # put the full results and the negative log10 values into a dataframe.
    pcadapt.results <- list("full_stats_pcadapt_4.0.3_bigsnpr_0.8.2" = full_stats_pcadapt_4.0.3_bigsnpr_0.8.2, 
                            "pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p" = pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p)
    
    if (length(rm.pops) != 0) {
        results.path <- file.path(data.path, 
                                  paste("/pcadaptAnalysis_excludingpops_", 
                                        paste(rm.pops, collapse = "_"), 
                                        sep = ""))
        if (!dir.exists(results.path)) {
            dir.create(results.path)
        }
        saveRDS(pcadapt.results, 
                paste(results.path, 
                      "/pcadaptResults_chr", 
                      chr.index, 
                      ".rds",
                      sep = ""))
    } else if (length(select.pops) != 0) {
        results.path <- file.path(data.path, 
                                  paste("/pcadaptAnalysis_excludingpops_", 
                                        paste(rm.pops, collapse = "_"), 
                                        sep = ""))
        if (!dir.exists(results.path)) {
            dir.create(results.path)
        }
        saveRDS(pcadapt.results, 
                paste(results.path, 
                      "/pcadaptAnalysis_chr", 
                      chr.index, 
                      ".rds", 
                      sep = ""))
    } else {
        results.path <- file.path(data.path, "/pcadaptAnalysis_allpops")
        if (!dir.exists(results.path)) {
            dir.create(results.path)
        }
        saveRDS(pcadapt.results, 
                paste(results.path, 
                      "/pcadaptResults_chr", 
                      chr.index, 
                      ".rds",
                      sep = ""))
    }
}

# Below is a wrapper function that will allow you to run all of the analyses
# that you would like while just pressing 'enter'. If you would like to run the
# analyses for all populations, leave the rm.pops and select.pops parameters as
# NULL. If you would like to run analyses where populations are removed or only
# certain populations are selected, you can concatenate those parameters into
# those parameters, respectively. However, do not put NULL into your
# concatenated parameters (i.e rm.pops = list(c("LM"), NULL, c("CL")) 

wrapper <- function (data.file = "data/large_data/genotypeMatrix.rds",
                     metadata = "data/modified_samplemetadata.csv",
                     rand.snps = "data/thinned_snps/thinnedMatrixAndMetaData1e+05Window.rds",
                     chr.index = seq(1, 10), 
                     rm.pops = NULL,
                     select.pops = NULL) {
    # Create folder structure for the analyses and their outputs depending upon
    # whether they remove or select certain populations
    
        if (length(select.pops) != 0) {
            data.path <- file.path("large_outputs", 
                                     "outlierAnalysis", 
                                     paste("outlierAnalysis_selectingpops_", 
                                            paste(select.pops, collapse="_"), sep=""),
                                     "data")
            plots.path <- file.path("large_outputs", 
                                     "outlierAnalysis", 
                                     paste("outlierAnalysis_selectingpops_", 
                                            paste(select.pops, collapse="_"), sep=""),
                                     "plots")                                 
        } else if (length(rm.pops) != 0) {
            data.path <- file.path("large_outputs", 
                                     "outlierAnalysis", 
                                     paste("outlierAnalysis_excludingpops_", 
                                            paste(rm.pops, collapse="_"), sep=""),
                                     "data")
            plots.path <- file.path("large_outputs", 
                                     "outlierAnalysis", 
                                     paste("outlierAnalysis_excludingpops_", 
                                            paste(rm.pops, collapse="_"), sep=""),
                                     "plots")
        } else {
            data.path <- file.path("large_outputs", 
                                   "outlierAnalysis", 
                                   "outlierAnalysis_allpops",
                                   "data")
            plots.path <- file.path("large_outputs", 
                                    "outlierAnalysis", 
                                    "outlierAnalysis_allpops",
                                    "plots")
        }
        if (!dir.exists(data.path)) {
            dir.create(data.path, recursive = TRUE)
        }
        if (!dir.exists(plots.path)) {
            dir.create(plots.path, recursive = TRUE)
        }

        fstByChromosome(data.file = data.file, 
                        chr.index = chr.index,  
                        rm.pops = rm.pops, 
                        select.pops = select.pops,
                        data.path = data.path, 
                        plots.path = plots.path)

        neutralParams(rand.snps = rand.snps, 
                      rm.pops = rm.pops, 
                      select.pops = select.pops,
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

        for (i in files) {
            outflank.outlierFinder(data.file = i, 
                                   null.dist = null.dist,
                                   data.path = data.path,
                                   plots.path = plots.path)
        }

#    if (length(rm.pops) != 0 | length(select.pops) != 0) {
#        for (i in chr.index) {
#                pcadaptAnalysis(chr.index = i, rm.pops = j)
#        }
#
#        for (i in chr.index) {
#            for (j in select.pops) {
#                pcadaptAnalysis(chr.index = i, rm.pops = j)
#            }
#        }
#    }
#    for (i in chr.index) {
#        pcadaptAnalysis(chr.index = i)
#    }    
}

wrapper()
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds")

wrapper(rm.pops="LM")
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        rm.pops="LM")

wrapper(select.pops = c("NEH", "UMFS"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("NEH", "UMFS"))

wrapper(select.pops = c("CL", "SL"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("CL", "SL"))

wrapper(select.pops = c("CLP", "HC_VA"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("CLP", "HC_VA"))

wrapper(select.pops = c("HP", "CS"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("HP", "CS"))

wrapper(select.pops = c("HG", "NEH"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("HG", "NEH"))

wrapper(select.pops = c("NG", "NEH"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("NG", "NEH"))

wrapper(select.pops = c("CS", "NEH"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("CS", "NEH"))

wrapper(select.pops = c("HC", "DEBY"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("HC", "DEBY"))

wrapper(select.pops = c("HC_VA", "DEBY"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("HC_VA", "DEBY"))

wrapper(select.pops = c("CL", "LOLA"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("CL", "LOLA"))

wrapper(select.pops = c("CLP", "LOLA"))
wrapper(rand.snps = "50KRandomSNPs1e+05Window.rds", 
        select.pops = c("CLP", "LOLA"))