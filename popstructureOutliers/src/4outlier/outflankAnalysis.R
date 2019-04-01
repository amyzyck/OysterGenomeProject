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

library(OutFLANK)  # outflank package
library(ggplot2)

runAnalysis <- function(data.file = "fstMatrix_chr1_minusLM.rds",
                        null.dist = "outflank_50KSubset_minusLM.rds") {
     # Creating folder for outlier plots. If the folder does not exist, create
     # it
     folder.path <- "outlierPlots"
     if (!dir.exists(folder.path)) {
          dir.create(folder.path)
     }
     # Load FST matrix from folder
     fst.stat <- readRDS(paste("fstMatricesByChr/", data.file, sep=""))
     # Get the name of the full FST matrix file that is being used.
     data.file <- strsplit(data.file, split=".rds")[[1]]
     # Get the random quasi-independent SNPs in order to be used to establish
     # neutral parameters
     outflank.data <- readRDS(null.dist)
     # Find outliers using the neutral parameters
     P1 <- pOutlierFinderChiSqNoCorr(fst.stat, 
                                     Fstbar = outflank.data$FSTNoCorrbar, 
                                     dfInferred = outflank.data$dfInferred, 
                                     qthreshold = 0.00001, 
                                     Hmin = 0.1)

     outlier <- P1$OutlierFlag == TRUE

     # png(paste(folder.path, "/outlierPlotUsing_", data.file , ".png", sep=""), 
     #     height = 512,
     #     width = 2048)
     # plot(P1$He, P1$FST, pch = 19, col = rgb(0,0,0,0.1))
     # points(P1$He[outlier], P1$FST[outlier], col = "blue")
     # dev.off()

     # hist(P1$pvaluesRightTail)
     
     png(paste(folder.path, "/manhattanPlotUsing_", data.file , ".png", sep=""), 
         height = 512,
         width = 2048)
     plot(P1$LocusName[P1$He>0.1], 
          P1$FST[P1$He>0.1],
          xlab="Position", 
          ylab="FST", 
          col=rgb(0,0,0,0.2),
          xaxt="n")
     axis(side = 1, at=seq(P1$LocusName[1], P1$LocusName[length(P1$LocusName)], by=100000))
     points(P1$LocusName[outlier], P1$FST[outlier], col="magenta", pch=20)
     dev.off()
}

wrapper <- function(data.folder="fstMatricesByChr") {
     # Check to see if the folder which has been provided has been created.
     if (!dir.exists(data.folder)) {
          stop("The folder holding the data files does not exist")
     # Check to see if there are any files in the folder provided. If there are
     # none, return with an error     
     } else if (length(list.files(path=data.folder)) == 0) {
          stop("There are no data files in the folder that you have provided")
     } else {
     # Run the function above across the files which the user has provided.
     }
     files <- list.files(path=data.folder)
     for (i in files) {
          runAnalysis(data.file = i)
     }
}
