###############################################################################
#
# File      : populationStructureScript.R 
# History   : 10/10/2018  Created by K Bodie Weedop (KBW)
#
###############################################################################
#
# This script is the full script going from the original VCF file
# SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz) to a
# genotype matrix, combining this matrix with metadata (loci positions,
# chromosomal locations and available data on individuals).
#
# After converting the matrix, the script will go through the process of
# performing a PCA on the full genotype matrix using the 'pcadapt' package.
#
# This script will also go through the process of using the 'bigsnpr' package to
# thin SNPs, remove long-range LD regions and perform a PCA on the data. 
#
###############################################################################


library(vcfR)
library(pcadapt)
library(bigsnpr)
library(ggplot2)

# read the VCF file using the 'vcfR' package
vcf <- read.vcfR('SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz')

# Load in metadata on individuals
metadata <- read.csv('modified_SampleMetaData.csv')

# Grab the metadata from the VCF file
geno <- vcf@gt[,-c(1)]
# Get the positions of the SNPs. Positions in bp.
POS <- getPOS(vcf)
# Get chromosome information.
CHR <- getCHROM(vcf)

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

# Factoring the different chromosomes so they can be encoded as integer values.
factor.chromosome <- factor(CHR)
allData <- list(G = G,
                positions = POS,
                chromosome = as.integer(factor.chromosome))

# Save allData for the sake of time when doing additional analysis
saveRDS(allData, 'matrixAndMetadata.rds')

# End of the conversion from VCF to genotype matrix.

###############################################################################

# Performing PCA on full genotype matrix using 'pcadapt' package.

# Given the format of the data, use the 'pcadapt' type to read the genotype
# matrix into the format necessary for other package functions
genoMatPcadapt <- read.pcadapt(allData$G, type='pcadapt')

fullPcadapt <- pcadapt(genoMatPcadapt, K=13)

# Save the PCA object for any future analyses.
saveRDS(fullPcadapt, 'pcadaptOnFullGenotypeMatrix.rds')

###############################################################################

# Performing SNP thinning, removal of long-range LD regions and PCA using the
# 'bigsnpr' package

# Convert to format (file-backed matrix; 'FBM') that Prive uses
# in his walk-through.
G_FBM <- add_code256(big_copy(t(allData$G), type="raw"),
                     code = bigsnpr:::CODE_012)

# Search all of the SNP sites where all individuals are heterozygotes (have a 
# value of 1)
allHeteroSites <- which(big_counts(G_FBM)[2,] == 90)

# Reset the variable G_FBM after taking out the sites where all individuals are
# heterozygotes.
G_FBM <- add_code256(big_copy(t(allData$G[-allHeteroSites, 1:90]), type="raw"),
                     code = bigsnpr:::CODE_012)

fbmAllData <- list(genotype = G_FBM,
                   fam = data.frame(pop.id = metadata$custom.id,
                                    color = metadata$color,
                                    plot.id = metadata$plot.id),
                   map = data.frame(positions = allData$positions[-allHeteroSites],
                                    chromosome = allData$chromosome[-allHeteroSites]))

# Run snp_autoSVD on the full dataset with 13 PCs 
fullDatasetSVD13 <- snp_autoSVD(G_FBM,
                                infos.chr = fbmAllData$map$chromosome,
                                infos.pos = fbmAllData$map$positions,
                                k = 13,
                                ncores=4)

saveRDS(fullDatasetSVD13, 'fullDatasetSVD13.rds')

# INSERT CODE FOR THE PCA PLOTS...

###############################################################################

# Perform analysis on subset (10K) of the thinned SNPs

thinnedMatrixAndMetaData <- list(G = allData$G[attr(fullDatasetSVD13, "subset"), ],
                                 positions = allData$positions[attr(fullDatasetSVD13, "subset")],
                                 chromosome = allData$chromosome[attr(fullDatasetSVD13, "subset")])

# Set seed and generate 10K random samples of loci to include in the analysis
# and subset the dataset.
set.seed(10)
random10 <- sample(1:nrow(thinnedMatrixAndMetaData$G), 10000)

random10AllData <- list(G = thinnedMatrixAndMetaData$G[random10, 1:90],
                           positions = thinnedMatrixAndMetaData$positions[random10],
                           chromosomes = thinnedMatrixAndMetaData$chromosome[random10])
# Write chromosome and positions of the thinned randomly selected SNPs to txt 
# file
write(paste(CHR[random10AllData$chromosome],
              random10AllData$positions, sep='    '),
              "10KRandomSNPs13PCs.txt")
# Get the subsetted data into the correct format.
subset10_G_FBM <- add_code256(big_copy(t(random10AllData$G), type="raw"),
                              code = bigsnpr:::CODE_012)
fbm.random10 <- list(genotype = subset10_G_FBM,
                       fam = data.frame(sample.id=metadata$Sample.ID,
                                        pop.id=metadata$Pop.ID,
                                        sex = metadata$Sex),
                       map = data.frame(positions = random10AllData$positions,
                                        chromosome = random10AllData$chromosome))
# Perform the clumping and check for long-range LD regions
subset10_fullDatasetSVD13 <- snp_autoSVD(subset10_G_FBM,
                                         infos.chr = random10AllData$chromosomes,
                                         infos.pos = random10AllData$positions,
                                         k = 13,
                                         ncores = 4)


saveRDS(subset10_fullDatasetSVD13, '10kRandomSubsetThinned13PCs.rds')

###############################################################################

# Perform analysis on subset (50K) of the thinned SNPs


set.seed(50)
random50 <- sample(1:nrow(thinnedMatrixAndMetaData$G), 50000)

random50_allData <- list(G = thinnedMatrixAndMetaData$G[random50, 1:90],
                         positions = thinnedMatrixAndMetaData$positions[random50],
                         chromosomes = thinnedMatrixAndMetaData$chromosome[random50])

subset50_G_FBM <- add_code256(big_copy(t(random50_allData$G), type="raw"),
                              code = bigsnpr:::CODE_012)

fbm.random50 <- list(genotype = subset50_G_FBM,
                     fam = data.frame(sample.id=metadata$Sample.ID,
                                      pop.id=metadata$Pop.ID,
                                      sex = metadata$Sex),
                     map = data.frame(positions = random50_allData$positions,
                                      chromosome = random50_allData$chromosome))

subset50_fullDatasetSVD13 <- snp_autoSVD(subset50_G_FBM,
                                         infos.chr = random50_allData$chromosome,
                                         infos.pos = random50_allData$positions,
                                         k = 13,
                                         ncores = 4)

saveRDS(subset50_fullDatasetSVD13, '50kRandomSubsetThinned13PCs.rds')