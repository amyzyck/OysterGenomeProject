###############################################################################
#
# File      : populationStructureScript.R 
# History   : 10/10/2018  Created by K Bodie Weedop (KBW)
#           : 10/11/2018 KBW - continued the analysis using OutFLANK
#           : 12/11/2018 KBW - consolidated code into a few functions to 
#                              streamline the process of data preparation.
#           : 12/17/2018 KBW - Wrapper added to run full analysis and 
#                              subsetting of data with a shell command.
#           : 02/26/2018 KBW - Refactor to provide the options to remove or 
#                              select for certain populations
#
###############################################################################

###############################################################################
#
# This script is the full script going from the original VCF file
# Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz) to a
# genotype matrix, combining this matrix with metadata (loci positions,
# chromosomal locations and available data on individuals).
#
# After converting the matrix, the script will go through the process of
# performing these various analyses:
#   - PCA on full genotype matrix (SNPs thinned and long range LD regions 
#     removed)
#   - Subset data to get 20K and 50K random SNPs from the preprocessed data
#   - Write all R objects and txt files for saving the data.
#
###############################################################################


library(vcfR)
library(pcadapt)
library(bigsnpr)

# Data preprocessing

# First thing to be done is converting the .vcf file to a genotype matrix which
# can be done using this function. However, this takes a while if you have a
# large .vcf file so my suggestion is to run this function using a shell script
# of some kind and have it run in the background on some server while you do
# other things. 

# If you do run it on a server, you will come back to a saved R object named
# 'genotypeMatrix.rds'. This is the .vcf file converted to a genotype matrix
# with sites where all individuals are heterozygotes removed.

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

getGenotypeMatrix <- function(vcf.file = NULL,
                              metadata.file = NULL,
                              subset.name = NULL) {
    
    metadata <- read.csv(metadata.file, stringsAsFactors = FALSE)

    if (!is.null(subset.name)) {
        if (!subset.name %in% names(metadata)) {
            stop("Error: the subset that you provided is not specified in the metadata.")
        }
    }

    # read the VCF file using the 'vcfR' package
    vcf <- read.vcfR(vcf.file)
    # Grab the metadata from the VCF file
    geno <- vcf@gt[,-c(1)]
    # Get the positions of the SNPs. Positions in bp.
    POS <- getPOS(vcf)
    assign("POS", POS, envir = .GlobalEnv)
    # Get chromosome information.
    CHR <- getCHROM(vcf)
    # Factoring the different chromosomes so they can be encoded as integer values.
    factor.chromosome <- as.integer(factor(CHR))

    assign("CHR", sort(unique(CHR)), envir = .GlobalEnv)

    # Matrix to hold binary values...
    G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

    # Changing the values to either 0, 1, 2
    G[grepl("0/0", geno, fixed = TRUE)] <- 0
    G[grepl("0|0", geno, fixed = TRUE)] <- 0
    G[grepl("0/1", geno, fixed = TRUE)] <- 1
    G[grepl("0|1", geno, fixed = TRUE)] <- 1
    G[grepl("1/0", geno, fixed = TRUE)] <- 1
    G[grepl("1|0", geno, fixed = TRUE)] <- 1
    G[grepl("1/1", geno, fixed = TRUE)] <- 2
    G[grepl("1|1", geno, fixed = TRUE)] <- 2

    # Next we will map the metadata onto the genotype matrix and combine them both
    # into a single object.

    if (!is.null(subset.name)) {
        # Select the indices that have been specified with a 1 in the metadata
        G <- G[, which(metadata[, subset.name] == 1)]
        Pop.ID <- metadata$Pop.ID[which(metadata[, subset.name] == 1)]
        Sample.ID <- metadata$Sample.ID[which(metadata[, subset.name] == 1)]
        # Name the file according to the subset that the user has provided.
        file.name <- paste("genotypeMatrix_", 
                            gsub("\\.", "_", subset.name), 
                            ".rds", sep="")
    } else {
        # Name the file in case that no subset.name is provided
        file.name <- "genotypeMatrix.rds"
    
        # Get population and sample IDs. May be reset if subset is provided.
        Pop.ID <- metadata$Pop.ID
        Sample.ID <- metadata$Sample.ID
    }

    # Remove sites where all individuals are heterozygotes.
    allFixedSites <- getFixedSites(G)

    allData <- list(G = G[-allFixedSites, ],
                    Pos = POS[-allFixedSites],
                    Chr = factor.chromosome[-allFixedSites],
                    Pop.ID = Pop.ID,
                    Sample.ID = Sample.ID)
    
    # End of the conversion from VCF to genotype matrix.
    # Save allData for the sake of time when doing additional analysis

    data.path <- file.path("data", "large_data")

    saveRDS(allData, paste(data.path, file.name, sep="/"))

    return(paste(data.path, file.name, sep="/"))
}

# Performing SNP thinning, removal of long-range LD regions and PCA using the
# 'bigsnpr' package

getThinnedSnps <- function(data = NULL,
                           subset.name = NULL,
                           pca.loadings = 13, 
                           window = 5000, 
                           nb.cores = 4){
    # Load the full data set
    allData <- readRDS(data)

    # Generate the path for the thinned SNPs if has yet to be created.
    data.path <- file.path("data", "thinned_snps")
    if (!dir.exists(data.path)) {
        dir.create(data.path)
    }

    # keep <- apply(allData$G, 1, function(x) length(unique(x[!is.na(x)])) != 1)
    # new.data$G <- allData$G[keep, ]
    # new.data$Pos <- allData$Pos[keep]
    # new.data$Chr <- allData$Chr[keep]

    # Code the genotype matrix so that it can be used in snp_autoSVD()
    G_FBM <- add_code256(big_copy(t(allData$G), type="raw"), code = bigsnpr:::CODE_012)
    
    # Run snp_autoSVD on the full dataset
    fullDatasetSVD <- snp_autoSVD(G_FBM,
                                  infos.chr = allData$Chr,
                                  infos.pos = allData$Pos,
                                  k = pca.loadings,
                                  size = window,
                                  is.size.in.bp = TRUE,
                                  ncores = nb.cores)
    
    
    if (!is.null(subset.name)) {
        file.name <- paste("thinnedMatrixAndMetaData",
                            window, 
                            "Window_", 
                            gsub("\\.", "_", subset.name), 
                            ".rds",
                            sep="")
        svd.name <- paste("fullDatasetSVDStats", 
                          window, 
                          "Window_", 
                          gsub("\\.", "_", subset.name), 
                          ".rds", 
                          sep = "")
    } else {
        file.name <- paste("thinnedMatrixAndMetaData", 
                            window, 
                            "Window_allpops.rds", 
                            sep="")
        svd.name <- paste("fullDatasetSVDStats", 
                          window, 
                          "Window_allpops.rds", 
                          sep = "")
    }
    
    # Save the SVD stats in a .rds file
    saveRDS(fullDatasetSVD, paste(data.path, svd.name, sep="/"))

    # Get the SNPs and their associated positions, chromosome locations and Pop.ID 
    thinnedMatrixAndMetaData <- list(G = allData$G[attr(fullDatasetSVD, "subset"), ],
                                     Pos = allData$Pos[attr(fullDatasetSVD, "subset")],
                                     Chr = allData$Chr[attr(fullDatasetSVD, "subset")],
                                     Pop.ID = allData$Pop.ID,
                                     Sample.ID = allData$Sample.ID)

    # Save the R object containin the data in a .rds file. 
    saveRDS(thinnedMatrixAndMetaData, paste(data.path, file.name, sep="/"))

    # Get the indices of the thinned SNPs and save them in a txt file
    write(paste(CHR[thinnedMatrixAndMetaData$Chr],
                thinnedMatrixAndMetaData$Pos, sep='    '),
                paste(data.path, gsub(".rds", ".txt", file.name), sep="/"))

    # Return the path and file name to use in next function
    return(paste(data.path, file.name, sep="/"))
}

###############################################################################

getSubsets <- function(data = NULL, 
                       window = NULL,
                       subset.name = NULL){

    postAnalysisData <- readRDS(data)

    data.path <- "data/thinned_snps/"

    if (!is.null(subset.name)) {
        name.50K <- paste("50KRandomSNPs",
                          window, 
                          "Window_", 
                          gsub("\\.", "_", subset.name), 
                          ".rds",
                          sep="")
        name.20K <- paste("20KRandomSNPs", 
                          window, 
                          "Window_", 
                          gsub("\\.", "_", subset.name), 
                          ".rds", 
                          sep = "")
    } else {
        name.50K <- paste("50KRandomSNPs", 
                            window, 
                            "Window_allpops.rds", 
                            sep="")
        name.20K <- paste("20KRandomSNPs", 
                          window, 
                          "Window_allpops.rds", 
                          sep = "")
    }

    # get subset (20K) of the thinned SNPs

    # Set seed and generate 20K random samples of loci and subset the dataset.
    set.seed(20)
    index <- 0
    random20 <- c(rep(0, 20000))
    for (i in unique(postAnalysisData$Chr)) {
        tmp <- floor(length(which(postAnalysisData$Chr == i)) * 0.1)
        random20[index + 1:tmp] <- sample(which(postAnalysisData$Chr == i), tmp)
        index <- index + tmp
    }

    while(length(random20) > 20000){
        t <- length(random20) - 20000
        chr <- 1
        index <- 1
        tmp <- floor(length(which(postAnalysisData$Chr == chr)) * 0.1)
        random20 <- random20[-index]
        chr =+ 1
        if (chr == 11) {
            chr <- 1
        }
        index =+ tmp
    }

    random20Subset <- list(G = postAnalysisData$G[random20, ],
                           Pos = postAnalysisData$Pos[random20],
                           Chr = postAnalysisData$Chr[random20],
                           Pop.ID = postAnalysisData$Pop.ID,
                           Sample.ID = postAnalysisData$Sample.ID)
    # Write chromosome and positions of the thinned randomly selected SNPs to txt 
    # file
    write(paste(CHR[random20Subset$Chr],
                random20Subset$Pos, sep='    '),
                paste(data.path, gsub(".rds", ".txt", name.20K), sep=""))

    # Write genotype matrix, positions and corresponding chromosome to R object file
    saveRDS(random20Subset, paste(data.path, name.20K, sep=""))

    
    set.seed(50)
    index <- 0
    random50 <- c(rep(0, 50000))
    for (i in unique(postAnalysisData$Chr)) {
        tmp <- floor(length(which(postAnalysisData$Chr == i)) * 0.24)
        random50[index + 1:tmp] <- sample(which(postAnalysisData$Chr == i), tmp)
        index <- index + tmp
    }

    while(length(random50) > 50000){
        t <- length(random50) - 50000
        chr <- 1
        index <- 1
        tmp <- floor(length(which(postAnalysisData$Chr == chr)) * 0.24)
        random50 <- random50[-index]
        chr =+ 1
        if (chr == 11) {
            chr <- 1
        }
        index =+ tmp
    }

    random50Subset <- list(G = postAnalysisData$G[random50, ],
                           Pos = postAnalysisData$Pos[random50],
                           Chr = postAnalysisData$Chr[random50],
                           Pop.ID = postAnalysisData$Pop.ID,
                           Sample.ID = postAnalysisData$Sample.ID)

    # Write chromosome and positions of the thinned randomly selected SNPs to txt 
    # file
    write(paste(CHR[random50Subset$Chr],
                random50Subset$Pos, sep='    '),
                paste(data.path, gsub(".rds", ".txt", name.50K), sep=""))

    # Write genotype matrix, positions and corresponding chromosome to R object file
    saveRDS(random50Subset, paste(data.path, name.50K, sep = ""))
}

###############################################################################

getDataWrapper <- function(vcf.file = "data/large_data/Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz",
                           metadata.file = "data/modified_samplemetadata.csv",
                           subset.name = subset.name,
                           window = 5000,
                           pca.loadings = 13,
                           nb.cores = 4){
    # Run the preprocessing function. This function will return the path to the
    # processed data to feed into the thinning function.
    data.path <- getGenotypeMatrix(vcf.file = vcf.file,
                                   metadata.file = metadata.file,
                                   subset.name = subset.name)
    # Run the thinning function on the genotype matrix data. The function will
    # also return the path to the data.
    thinned.path <- getThinnedSnps(data = data.path,
                                   subset.name = subset.name,
                                   pca.loadings = pca.loadings, 
                                   window = window, 
                                   nb.cores = nb.cores)

    getSubsets(data = thinned.path, 
               window = window,
               subset.name = subset.name)
}

getDataWrapper(subset.name = "exclude.LM", window = 5000)
getDataWrapper(subset.name = "exclude.LM", window = 50000)
getDataWrapper(subset.name = "unrelated", window = 5000)
getDataWrapper(subset.name = "unrelated", window = 50000)