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

library(OutFLANK)

preprocess <- function(data.file = "50KRandomSNPs1e+05Window.rds", 
                       exclude.pops = "") {
    data <- readRDS(data.file)

    # Need to remove any populations that should be excluded from the neutral
    # paramaterization.
    if (exclude.pops != "") {
        for (i in exclude.pops) {
            to.delete <- which(data$sample.id == i)
            data$G <- data$G[,-to.delete]
            data$sample.id <- data$sample.id[-to.delete]
        }
    } else {
        
    }

    if (length(data$sample.id) != ncol(data$G)) {
        stop("Number of individuals in dataset and labels does not match")
    }

    return(data)
}

getNeutralParams <- function (data.file = "50KRandomSNPs1e+05Window.rds", 
                              exclude.pops = "") {
    folder.path <- paste("neutralParamPlotsUsing", strsplit(data.file, split=".rds")[[1]], sep="")
    if (!dir.exists(folder.path)) {
        dir.create(folder.path)
    }
    data.subset <- preprocess(data.file, exclude.pops)

    # getting the number of populations in order to inform OUTFLANK
    pops <- length(unique(data.subset$sample.id))

    fst.stat <- MakeDiploidFSTMat(t(data.subset$G), 
                                  locusNames = data.subset$positions, 
                                  popNames = data.subset$sample.id)

    outflank.data <- OutFLANK(fst.stat, 
                              NumberOfSamples = pops, 
                              qthreshold = 0.05, 
                              Hmin = 0.1)

    if (exclude.pops != "") {

        plot(fst.stat$He, fst.stat$FST)
        dev.copy(png, paste(folder.path, "/heteroVsFst_minus", exclude.pops, ".png", sep=""))
        dev.off()

        plot(fst.stat$FST, fst.stat$FSTNoCorr)
        abline(0,1)
        dev.copy(png, paste(folder.path, "/fstNoCorrVsFst_minus", exclude.pops, ".png", sep=""))
        dev.off()

        OutFLANKResultsPlotter(outflank.data, 
                               withOutliers = TRUE, 
                               NoCorr = TRUE, 
                               Hmin = 0.1, 
                               binwidth = 0.001, 
                               Zoom = FALSE, 
                               RightZoomFraction = 0.05, 
                               titletext = NULL)
        dev.copy(png, paste(folder.path, "/outflankResults_minus", exclude.pops, ".png", sep=""))
        dev.off()

        OutFLANKResultsPlotter(outflank.data, 
                       withOutliers = TRUE, 
                       NoCorr = TRUE, 
                       Hmin = 0.1, 
                       binwidth = 0.001, 
                       Zoom = TRUE, 
                       RightZoomFraction = 0.15, 
                       titletext = NULL)
        dev.copy(png, paste(folder.path, "/outflankResultsRightTail_minus", exclude.pops, ".png", sep=""))
        dev.off()

        saveRDS(fst.stat, paste("fstMatrix_50KSubset_minus", exclude.pops, ".rds", sep=""))
        saveRDS(outflank.data, paste("outflank_50KSubset_minus", exclude.pops, ".rds", sep=""))
    } else {
        plot(fst.stat$He, fst.stat$FST)
        dev.copy(png, paste(folder.path, "/heteroVsFst_allpops.png", sep=""))
        dev.off()

        plot(fst.stat$FST, fst.stat$FSTNoCorr)
        abline(0,1)
        dev.copy(png, paste(folder.path, "/fstNoCorrVsFst_allpops.png", sep=""))
        dev.off()

        OutFLANKResultsPlotter(outflank.data, 
                               withOutliers = TRUE, 
                               NoCorr = TRUE, 
                               Hmin = 0.1, 
                               binwidth = 0.001, 
                               Zoom = FALSE, 
                               RightZoomFraction = 0.05, 
                               titletext = NULL)
        dev.copy(png, paste(folder.path, "/outflankResults_allpops.png", sep=""))
        dev.off()

        OutFLANKResultsPlotter(outflank.data, 
                       withOutliers = TRUE, 
                       NoCorr = TRUE, 
                       Hmin = 0.1, 
                       binwidth = 0.001, 
                       Zoom = TRUE, 
                       RightZoomFraction = 0.15, 
                       titletext = NULL)
        dev.copy(png, paste(folder.path, "/outflankResultsRightTail_allpops.png", sep=""))
        dev.off()

        saveRDS(fst.stat, "fstMatrix_50KSubset_allpops.rds")
        saveRDS(outflank.data, "outflank_50KSubset_allpops.rds")
    }
}

wrapper <- function(pop.list = list("LM", "")) {
    for (i in pop.list) {
        getNeutralParams(exclude.pops = i)
    }
}

wrapper()