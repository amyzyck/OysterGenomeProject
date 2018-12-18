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
                                                       color = metadata$color,
                                                       plot.id = metadata$plot.id),
                                      map = data.frame(positions = allData$positions,
                                                       chromosome = allData$chromosome))

    saveRDS(genotypeMatrixAndMetadata, 'genotypeMatrixAndMetadata.rds')
}

# Performing SNP thinning, removal of long-range LD regions and PCA using the
# 'bigsnpr' package

analysisFullDataset <- function(data = "genotypeMatrixAndMetadata.rds",
                                pca.loadings = 13, 
                                window = 10000, 
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

    saveRDS(fullDatasetSVD, 'fullDatasetSVD.rds')

    thinnedMatrixAndMetaData <- list(G = allData$genotype[attr(fullDatasetSVD, "subset"), ],
                                 positions = allData$map$positions[attr(fullDatasetSVD, "subset")],
                                 chromosome = allData$map$chromosome[attr(fullDatasetSVD, "subset")])

    saveRDS(thinnedMatrixAndMetaData, "thinnedMatrixAndMetaData.rds")
}

# INSERT CODE FOR THE PCA PLOTS...

###############################################################################

getSubsets <- function(data = "thinnedMatrixAndMetaData.rds"){

    # get subset (20K) of the thinned SNPs

    postAnalysisData <- readRDS(data)

    # Set seed and generate 20K random samples of loci and subset the dataset.
    set.seed(20)
    random20 <- sample(1:nrow(postAnalysisData$G), 20000)

    random20Subset <- list(G = postAnalysisData$G[random20, 1:90],
                            positions = postAnalysisData$positions[random20],
                            chromosomes = postAnalysisData$chromosome[random20])
    # Write chromosome and positions of the thinned randomly selected SNPs to txt 
    # file
    write(paste(CHR[random20Subset$chromosome],
                random20Subset$positions, sep='    '),
                "20KRandomSNPs.txt")

    set.seed(50)
    random50 <- sample(1:nrow(postAnalysisData$G), 50000)

    random50Subset <- list(G = postAnalysisData$G[random50, 1:90],
                           positions = postAnalysisData$positions[random50],
                           chromosomes = postAnalysisData$chromosome[random50])

    # Write chromosome and positions of the thinned randomly selected SNPs to txt 
    # file
    write(paste(CHR[random50Subset$chromosome],
                random50Subset$positions, sep='    '),
                "50KRandomSNPs.txt")
}

###############################################################################

getDataWrapper <- function(){
    getProcessedData(vcf.file = "Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz",
                     metadata.file = "modifiedColors_SampleMetaData.csv")

    analysisFullDataset(data = "genotypeMatrixAndMetadata.rds",
                        pca.loadings = 13, 
                        window = 10000, 
                        nb.cores = 4)

    getSubsets(data = "thinnedMatrixAndMetaData.rds")
}

getDataWrapper()