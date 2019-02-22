###############################################################################
#
# File      : populationStructureScript.R 
# History   : 10/10/2018  Created by K Bodie Weedop (KBW)
#           : 10/11/2018 KBW - continued the analysis using OutFLANK
#           : 12/11/2018 KBW - consolidated code into a few functions to 
#                              streamline the process of data preparation.
#           : 12/17/2018 KBW - Wrapper added to run full analysis and 
#                              subsetting of data with a shell command. 
#
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
library(ggplot2)

# Data preprocessing

# First thing to be done is converting the .vcf file to a genotype matrix which
# can be done using this function. However, this takes a while if you have a
# large .vcf file so my suggestion is to run this function using a shell script
# of some kind and have it run in the background on some server while you do
# other things. 

# If you do run it on a server, you will come back to a saved R object named
# 'genotypeMatrix.rds'. This is the .vcf file converted to a genotype matrix
# with sites where all individuals are heterozygotes removed.

getFixedHeterozygoteSites <- function(genotype.matrix) {
    # Convert to format (file-backed matrix; 'FBM') that Prive uses
    # in his walk-through.
    G_FBM <- add_code256(big_copy(t(genotype.matrix), type="raw"),
                         code = bigsnpr:::CODE_012)

    # Search all of the SNP sites where all individuals are heterozygotes (have
    # a value of 1)
    return(which(big_counts(G_FBM)[2,] == 90))
}

getProcessedData <- function(vcf.file = "Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz",
                             metadata.file = "modifiedColors_SampleMetaData.csv") {
    # read the VCF file using the 'vcfR' package
    vcf <- read.vcfR(vcf.file)
    # Grab the metadata from the VCF file
    geno <- vcf@gt[,-c(1)]
    # Get the positions of the SNPs. Positions in bp.
    POS <- getPOS(vcf)
    assign("POS", POS, envir = .GlobalEnv)
    # Get chromosome information.
    CHR <- getCHROM(vcf)
    assign("CHR", CHR, envir = .GlobalEnv)
    # Factoring the different chromosomes so they can be encoded as integer values.
    factor.chromosome <- as.integer(factor(CHR))

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

    # Remove sites where all individuals are heterozygotes.
    allFixedHeteroSites <- getFixedHeterozygoteSites(G)

    allData <- list(genotype = G[-allFixedHeteroSites, 1:90],
                    positions = POS[-allFixedHeteroSites],
                    chromosome = factor.chromosome[-allFixedHeteroSites])
    
    # End of the conversion from VCF to genotype matrix.
    # Save allData for the sake of time when doing additional analysis
    saveRDS(allData, 'genotypeMatrix.rds')

    # Next we will map the metadata onto the genotype matrix and combine them both
    # into a single object.
    metadata <- read.csv(metadata.file)
    genotypeMatrixAndMetadata <- list(genotype = allData$genotype,
                                      fam = data.frame(pop.id = metadata$custom.id,
                                                       plot.id = metadata$plot.id),
                                      map = data.frame(positions = allData$positions,
                                                       chromosome = allData$chromosome))

    saveRDS(genotypeMatrixAndMetadata, 'genotypeMatrixAndMetadata.rds')
}

# Performing SNP thinning, removal of long-range LD regions and PCA using the
# 'bigsnpr' package

analysisFullDataset <- function(data = "genotypeMatrixAndMetadata.rds",
                                pca.loadings = 13, 
                                window = 100000, 
                                nb.cores = 4){
    allData <- readRDS(data)

    G_FBM <- add_code256(big_copy(t(allData$genotype), type="raw"),
                         code = bigsnpr:::CODE_012)
    
    # Run snp_autoSVD on the full dataset
    fullDatasetSVD <- snp_autoSVD(G_FBM,
                                  infos.chr = allData$map$chromosome,
                                  infos.pos = allData$map$positions,
                                  k = pca.loadings,
                                  size = window,
                                  is.size.in.bp = TRUE,
                                  ncores = nb.cores)

    saveRDS(fullDatasetSVD, paste("fullDatasetSVD", window, "Window.rds", sep = ""))

    thinnedMatrixAndMetaData <- list(G = allData$genotype[attr(fullDatasetSVD, "subset"), ],
                                     positions = allData$map$positions[attr(fullDatasetSVD, "subset")],
                                     chromosome = allData$map$chromosome[attr(fullDatasetSVD, "subset")])

    saveRDS(thinnedMatrixAndMetaData, 
            paste("thinnedMatrixAndMetaData", window, "Window.rds", sep = ""))
}

# INSERT CODE FOR THE PCA PLOTS...

###############################################################################

getSubsets <- function(data = "thinnedMatrixAndMetaData1e+05Window.rds"){

    window <- strsplit(data, "Data", fixed = TRUE)[[1]][2]

    postAnalysisData <- readRDS(data)

    # get subset (20K) of the thinned SNPs

    # Set seed and generate 20K random samples of loci and subset the dataset.
    set.seed(20)
    index <- 0
    random20 <- c(rep(0, 20000))
    for (i in unique(postAnalysisData$chromosome)) {
        tmp <- floor(length(which(postAnalysisData$chromosome == i)) * 0.1)
        random20[index + 1:tmp] <- sample(which(postAnalysisData$chromosome == i), tmp)
        index <- index + tmp
    }

    while(length(random20) > 20000){
        t <- length(random20) - 20000
        chr <- 1
        index <- 1
        tmp <- floor(length(which(postAnalysisData$chromosome == chr)) * 0.1)
        random20 <- random20[-index]
        chr =+ 1
        if (chr == 11) {
            chr <- 1
        }
        index =+ tmp
    }

    random20Subset <- list(G = postAnalysisData$G[random20, 1:90],
                           positions = postAnalysisData$positions[random20],
                           chromosomes = postAnalysisData$chromosome[random20])
    # Write chromosome and positions of the thinned randomly selected SNPs to txt 
    # file
    write(paste(CHR[random20Subset$chromosome],
                random20Subset$positions, sep='    '),
                paste("20KRandomSNPs", strsplit(window, ".rds")[[1]][1], ".txt", sep = ""))

    # Write genotype matrix, positions and corresponding chromosome to R object file
    saveRDS(random20Subset, paste("20KRandomSNPs", window, sep = ""))

    
    set.seed(50)
    index <- 0
    random50 <- c(rep(0, 50000))
    for (i in unique(postAnalysisData$chromosome)) {
        tmp <- floor(length(which(postAnalysisData$chromosome == i)) * 0.24)
        random50[index + 1:tmp] <- sample(which(postAnalysisData$chromosome == i), tmp)
        index <- index + tmp
    }

    while(length(random50) > 50000){
        t <- length(random50) - 50000
        chr <- 1
        index <- 1
        tmp <- floor(length(which(postAnalysisData$chromosome == chr)) * 0.24)
        random50 <- random50[-index]
        chr =+ 1
        if (chr == 11) {
            chr <- 1
        }
        index =+ tmp
    }

    random50Subset <- list(G = postAnalysisData$G[random50, 1:90],
                           positions = postAnalysisData$positions[random50],
                           chromosomes = postAnalysisData$chromosome[random50])

    # Write chromosome and positions of the thinned randomly selected SNPs to txt 
    # file
    write(paste(CHR[random50Subset$chromosome],
                random50Subset$positions, sep='    '),
                paste("50KRandomSNPs", strsplit(window, ".rds")[[1]][1], ".txt", sep = ""))

    # Write genotype matrix, positions and corresponding chromosome to R object file
    saveRDS(random50Subset, paste("50KRandomSNPs", window, sep = ""))

}

###############################################################################

getDataWrapper <- function(){
    getProcessedData(vcf.file = "Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz",
                     metadata.file = "modified_samplemetadata.csv")

    analysisFullDataset(data = "genotypeMatrixAndMetadata.rds",
                        pca.loadings = 13, 
                        window = 100000, 
                        nb.cores = 4)

    getSubsets(data = "thinnedMatrixAndMetaData1e+05Window.rds")
}

getDataWrapper()