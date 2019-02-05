###############################################################################
#
# File      : fstMatrixByChromosome.R 
# History   : 01/30/2019  Created by K Bodie Weedop (KBW)
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

preprocess <- function(data = "genotypeMatrix.rds", 
                       chr.index = 1, 
                       exclude.pops = "") {
    # Load the genotype matrix that has been provided
    data <- readRDS(data)

    # Selecting just one chromosome to decouple the data processing
    one.chr <- list(genotype = data$genotype[which(data$chromosome == chr.index),],
                    position = data$positions[which(data$chromosome == chr.index)],
                    chromosome = data$chromosome[which(data$chromosome == chr.index)],
                    sample.id = data$sample.id)
    # If the user has provided any populations to remove, select it and remove
    # the column associated with it from the genotype matrix
    if (exclude.pops != "") {
        for (i in exclude.pops) {
            to.delete <- which(data$sample.id == i)
            one.chr$genotype <- one.chr$genotype[,-to.delete]
            one.chr$sample.id <- data$sample.id[-to.delete]
        }
    } else {
        # Pass if no populations have been provided
    }
    # Throw an error if the number of samples does not correspond to the number
    # of columns provided.
    if (length(one.chr$sample.id) != ncol(one.chr$genotype)) {
        stop("Number of individuals in dataset and labels does not match")
    }
    # Return the processed data
    return(one.chr)
}

wrapper <- function(data = "genotypeMatrix.rds", 
                    chr.index = seq(1,10),  
                    delete.pops = "") {
    folder.path <- "fstMatricesByChr"
    if (!dir.exists(folder.path)) {
        dir.create(folder.path)
    }
    for (i in chr.index) {
        preprocessedData <- preprocess(data, chr.index = i, exclude.pops = delete.pops)
        fst.stat <- MakeDiploidFSTMat(t(preprocessedData$genotype), 
                                      locusNames = preprocessedData$position, 
                                      popNames = preprocessedData$sample.id)
        if (delete.pops == "") {
            saveRDS(fst.stat, paste(folder.path, "/fstMatrix_chr", i, "_allpops.rds", sep=""))
        } else if (delete.pops != "") {
           saveRDS(fst.stat, paste(folder.path, "/fstMatrix_chr", i, "_minus", delete.pops, ".rds", sep=""))
        }
    }
}

