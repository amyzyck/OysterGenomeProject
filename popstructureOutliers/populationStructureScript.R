library(vcfR)
library(bigsnpr)
library(ggplot2)

setwd('~/popstructureOutliers/')
vcf <- read.vcfR('SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz')
# Load in metadata
metadata <- read.csv('modified_SampleMetaData.csv')


geno <- vcf@gt[,-c(1)]
# Get the positions of the SNPs. Positions in bp.
positions <- getPOS(vcf)
# Get chromosome information.
chromosome <- getCHROM(vcf)

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
factor.chromosome <- factor(chromosome)
allData <- list(G = G,
                positions = positions,
                chromosome = as.integer(factor.chromosome))
# Save allData for the sake of time if necessary
saveRDS(allData, 'matrixAndMetadata.rds')

# Convert to format (file-backed matrix; 'FBM') that Prive uses
# in his walk-through.
G_FBM <- add_code256(big_copy(t(allData$G),
                                type="raw"), code = bigsnpr:::CODE_012)

# Search all of the SNP sites where all individuals are heterozygotes (have a 
# value of 1)
heteros <- which(big_counts(G_FBM)[2,] == 90)

# Reset the variable G_FBM after taking out the sites where all individuals are
# heterozygotes.
G_FBM <- add_code256(big_copy(t(allData$G[-heteros, 1:90]), type="raw"),
                     code = bigsnpr:::CODE_012)

fbm.allData <- list(genotype = G_FBM,
                    fam = data.frame(pop.id = metadata$custom.id,
                                     color = metadata$color,
                                     plot.id = metadata$plot.id),
                    map = data.frame(positions = allData$positions[-heteros],
                                     chromosome = allData$chromosome[-heteros]))

fullDatasetSVD13 <- snp_autoSVD(G_FBM,
                              infos.chr = allData$chromosome,
                              infos.pos = allData$positions,
                              k = 13,
                              ncores=4)

# INSERT CODE FOR THE PCA PLOTS...

# Same analysis but just with a random sample of 10,000 Loci.

# I need to first remove the sites where all individuals are heterozygotes so
# they are not included in the random sample.
noheteros_allData <- list(G = allData$G[-heteros,],
                          positions = allData$positions[-heteros],
                          chromosomes = allData$positions[-heteros])
# Set seed and generate 10K random samples of loci to include in the analysis 
# and subset the dataset.
set.seed(10)
random10 <- sample(1:nrow(noheteros_allData$G), 10000)
random10_allData <- list(G = noheteros_allData$G[random10, 1:90],
                           positions = noheteros_allData$positions[random10],
                           chromosomes = noheteros_allData$positions[random10])

# Get the subsetted data into the correct format.
subset10_G_FBM <- add_code256(big_copy(t(random10_allData$G), type="raw"),
                              code = bigsnpr:::CODE_012)

# Perform the clumping and check for long-range LD regions
> subset10_svd13 <- snp_autoSVD(subset10_G_FBM,
                                infos.chr = random10_allData$chromosomes,
                                infos.pos = random10_allData$positions,
                                k = 13,
                                ncores = 4)

# INSERT CODE FOR THE PCA PLOTS...

# Same analysis but just with a random sample of 50,000 Loci.

set.seed(50)
random50 <- sample(1:nrow(noheteros_allData$G), 50000)
random50_allData <- list(G = noheteros_allData$G[random50, 1:90],
                         positions = noheteros_allData$positions[random50],
                         chromosomes = noheteros_allData$positions[random50])

subset50_G_FBM <- add_code256(big_copy(t(random50_allData$G), type="raw"),
                              code = bigsnpr:::CODE_012)

subset50_svd13 <- snp_autoSVD(subset50_G_FBM,
                              infos.chr = random50_allData$chromosomes,
                              infos.pos = random50_allData$positions,
                              k = 13,
                              ncores = 4)

# INSERT CODE FOR PLOTS HERE...
