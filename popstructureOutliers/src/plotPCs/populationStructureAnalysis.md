# Detecting Population Structure with PCA Using SNP Data

## Remove Tri-allelic Sites

Load the original data from KITT server.

I made a separate file on the server and copied the files into that folder to
not alter the original data mistakenly.

```R
> setwd('~/popstructureOutliers/')
> library(vcfR)
> vcf <- read.vcfR('SNP.TRSdp5g95FnDNAmaf05.vcf.gz')
Scanning file to determine attributes.
File attributes:
  meta lines: 63
  header_line: 64
  variant count: 7030071
  column count: 99
Meta line 63 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 7030071
  Character matrix gt cols: 99
  skip: 0
  nrows: 7030071
  row_num: 0
Processed variant: 7030071
All variants processed
```

Remove 'FORMAT' column...

```R
> geno <- vcf@gt[,-c(1)]
```

Matrix to hold binary values...

```R
> G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
```

Changing the values to either 0, 1, 2...

```R
> G[grepl("0/0", geno, fixed = TRUE)] <- 0
> G[grepl("0|0", geno, fixed = TRUE)] <- 0

> G[grepl("0/1", geno, fixed = TRUE)] <- 1
> G[grepl("0|1", geno, fixed = TRUE)] <- 1
> G[grepl("1/0", geno, fixed = TRUE)] <- 1
> G[grepl("1|0", geno, fixed = TRUE)] <- 1

> G[grepl("1/1", geno, fixed = TRUE)] <- 2
> G[grepl("1|1", geno, fixed = TRUE)] <- 2
```

There are some values are still NAs after doing these steps. Let's take a look
at some of the values that are not affected by these last conversion steps...

```R
> geno[11,1]
CL_1
2|2:9:0,0,0,8,0,0,0,0,1:0:0:0,0,8,0,0,0,0,1:0,0,328,0,0,0,0,15:-24.8931,
-24.8931,-24.8931,-24.8931,-24.8931,-24.8931,-2.40824,-2.40824,-2.40824,0,
-24.8931,-24.8931,-24.8931,-2.40824,-24.8931,-24.8931,-24.8931,-24.8931,
-2.40824,-24.8931,-24.8931,-24.8931,-24.8931,-24.8931,-2.40824,-24.8931,
-24.8931,-24.8931,-24.8931,-24.8931,-24.8931,-2.40824,-24.8931,-24.8931,
-24.8931,-24.8931,-24.1335,-24.1335,-24.1335,-1.50338,-24.1335,-24.1335,
-24.1335,-24.1335,-23.8325

> geno[204,1]
CL_1
./.:3:3,0,0,0,0:3:110:0,0,0,0:0,0,0,0:0,-0.90309,-9.52694,-0.90309,-9.52694,
-9.52694,-0.90309,-9.52694,-9.52694,-9.52694,-0.90309,-9.52694,-9.52694,
-9.52694,-9.52694

> geno[205,1]
CL_1
./.:3:3,0,0,0,0:3:110:0,0,0,0:0,0,0,0:0,-0.90309,-9.52694,-0.90309,-9.52694,
-9.52694,-0.90309,-9.52694,-9.52694,-9.52694,-0.90309,-9.52694,-9.52694,
-9.52694,-9.52694
```

We can take a look at how many occurrences of these tri-allelic sites there are
using zgrep in the terminal.

```shell
[cvsnps@KITT VCF_files]$ zgrep -c "|2" SNP.TRSdp5g95FnDNAmaf05.vcf.gz
243375
[cvsnps@KITT VCF_files]$ zgrep -c "2|" SNP.TRSdp5g95FnDNAmaf05.vcf.gz
206620
[cvsnps@KITT VCF_files]$ zgrep -c "2/" SNP.TRSdp5g95FnDNAmaf05.vcf.gz
209893
[cvsnps@KITT VCF_files]$ zgrep -c "/2" SNP.TRSdp5g95FnDNAmaf05.vcf.gz
206225
```

So, for the time being we decided to remove sites which were tri-allelic. We
could do this using the terminal or VCFtools has a feature where you can
actually filter sites which have above a certain amount of alleles.

VCFtools manual:

```
--min-alleles <integer>
--max-alleles <integer>

Include only sites with a number of alleles greater than or equal to the "--min-alleles" value and less than or equal to the "--max-alleles" value. One of these options may be used without the other. 
For example, to include only bi-allelic sites, one could use:

vcftools --vcf file1.vcf --min-alleles 2 --max-alleles 2
```

```shell
[cvsnps@KITT VCF_files]$ vcftools --gzvcf SNP.TRSdp5g95FnDNAmaf05.vcf.gz --min-alleles 0 --max-alleles 2 --out SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2
```

Just to make sure that this command is doing what I think it is doing I want to
zgrep once again those sites patterns that we did earlier

```shell
[cvsnps@KITT VCF_files]$ zgrep -c "2|" SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz
0
[cvsnps@KITT VCF_files]$ zgrep -c "|2" SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz
0
[cvsnps@KITT VCF_files]$ zgrep -c "/2" SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz
0
[cvsnps@KITT VCF_files]$ zgrep -c "2/" SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz
0
```

From the looks of it, the option built in to VCFtools does what we need it to
but now we can try it out in R

Load the data just like before using vcfR...

```R
> setwd('~/popstructureOutliers/')
> library(vcfR)
> vcf <- read.vcfR('SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz')
```

Remove 'FORMAT' column...

```R
> geno <- vcf@gt[,-c(1)]
```

Matrix to hold binary values...

```R
> G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
```

Changing the values to either 0, 1, 2

```R
> G[grepl("0/0", geno, fixed = TRUE)] <- 0
> G[grepl("0|0", geno, fixed = TRUE)] <- 0
> G[grepl("0/1", geno, fixed = TRUE)] <- 1
> G[grepl("0|1", geno, fixed = TRUE)] <- 1
> G[grepl("1/0", geno, fixed = TRUE)] <- 1
> G[grepl("1|0", geno, fixed = TRUE)] <- 1
> G[grepl("1/1", geno, fixed = TRUE)] <- 2
> G[grepl("1|1", geno, fixed = TRUE)] <- 2
> saveRDS(G, "GMatrix_TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.rds")
```

```R
> dim(G)
[1] 6782970      90
```

All sites and individuals are accounted for from the original .vcf file.

```R
> head(G[1:10,1:10])

     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    0    0    0    0    0    0    0    0    0     0
[2,]    0    0    0    1    0    0    1    1    1     1
[3,]    0    0    0    1    0    0    1    1    1     1
[4,]    0    0    0    1    0    0    2    1    1     1
[5,]    0    0    0    0    0    0    0    0    0     0
[6,]    0    0    0    0    0    0    0    0    0     0
```

```R
> rem = c(which(rowSums(G)==0), which(rowSums(G-2)==0)) ##fixed loci
> rem
integer(0)

> position[rem]
integer(0)
```

No fixed loci?

## 2. Filter For Sites That Have 0% Missing

VCFtools also provides a feature for doing this sort of site filtering. 

VCFtools manual:

```
--max-missing <float>

Exclude sites on the basis of the proportion of missing data (defined to be
between 0 and 1, where 0 allows sites that are completely missing and 1
indicates no missing data allowed).
```

We can do it on the file that we filtered for the tri-allelic site and output
the result in another file.

```shell
[cvsnps@KITT popstructureOutliers] vcftools --gzvcf SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out SNP.TRSdp5g95FnDNAmaf05_min0_max2_noMissing --stdout | pigz -c > SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz

VCFtools - 0.1.15
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --gzvcf SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz
        --recode-INFO-all
        --max-missing 1
        --out SNP.TRSdp5g95FnDNAmaf05_min0_max2_noMissing
        --recode
        --stdout

Using zlib version: 1.2.7
After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 3211528 out of a possible 6782970 Sites
Run Time = 1656.00 seconds
```

Read the data into R using vcfR package

```R
> vcf <- read.vcfR('SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz')
```

Remove 'FORMAT' column...

```R
> geno <- vcf@gt[,-c(1)]
# Get the positions of the SNPs. Positions in bp.
> positions <- getPOS(vcf)
# Get chromosome information.
> chromosome <- getCHROM(vcf)
```

Matrix to hold binary values...

```R
> G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
```

Changing the values to either 0, 1, 2

```R
> G[grepl("0/0", geno, fixed = TRUE)] <- 0
> G[grepl("0|0", geno, fixed = TRUE)] <- 0
> G[grepl("0/1", geno, fixed = TRUE)] <- 1
> G[grepl("0|1", geno, fixed = TRUE)] <- 1
> G[grepl("1/0", geno, fixed = TRUE)] <- 1
> G[grepl("1|0", geno, fixed = TRUE)] <- 1
> G[grepl("1/1", geno, fixed = TRUE)] <- 2
> G[grepl("1|1", geno, fixed = TRUE)] <- 2
> saveRDS(G, "GMatrix_TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.rds")
```

Combine the matrix, positions and chromosome (encoded as integer values)
information all into one data structure.

```R
# Factoring the different chromosomes so they can be encoded as integer values.
> factor.chromosome <- factor(chromosome)
> allData <- list(G = G,
                  positions = positions,
                  chromosome = as.integer(factor.chromosome))
saveRDS(allData, 'matrixAndMetadata.rds')
```

## Thin SNPs for LD Based on Recommendation of (Privé et al. 2018) Using 
## 'bigsnpR' Package

Now we need it to be in the format (file-backed matrix; 'FBM') that Prive uses
in his walk-through.

```R
> allData <- readRDS('matrixAndMetadata.rds')
> G_FBM <- add_code256(big_copy(t(allData$G),
                                type="raw"), code = bigsnpr:::CODE_012)
# Puts it in the raw format and stores likelihood genotype probability
> dim(G_FBM)
  [1]      90 3211528
```

Now I can add in the metadata with this data structure...

```R
> metadata <- read.csv('SampleMetaData.csv')
> fbm.allData <- list(genotype = G_FBM,
                      fam = data.frame(pop.id = metadata$custom.id,
                                       color = metadata$color,
                                       plot.id = metadata$plot.id),
                      map = data.frame(positions = allData$positions,
                                       chromosome = allData$chromosome))
```

### When pruning only

Well, Prive is using clumping on the Minor Allele Frequencies (MAF) instead of
pruning.

```R
> svd0 <- snp_autoSVD(G_FBM,
                      infos.chr = allData$chromosome,
                      infos.pos = allData$positions,
                      ncores=4)
Phase of clumping at r2 > 0.2.. keep 480417 SNPs.

Iteration 1:
Computing SVD..
Error in na.fail.default(data) : missing values in object
```

Okay so let's check for missing data or NAs in the data.

```R
big_counts(G_FBM, ind.col=1:10)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
0      67   43    3   83   83   51   79   61   76    60
1      23   47   32    4    4   35   11   21   13    22
2       0    0   55    3    3    4    0    8    1     8
<NA>    0    0    0    0    0    0    0    0    0     0

# Okay but we have WAY more columns than just 10...

> missing.data <- big_counts(G_FBM)
> str(missing.data)
 int [1:4, 1:3211528] 67 23 0 0 43 47 0 0 3 32 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:4] "0" "1" "2" NA
  ..$ : NULL
> sum(missing.data[4,])
[1] 0
```

So no NAs in the character matrix. How about in the chromosomes or the
positions

```R
> any(is.na(allData$chromosome))
[1] FALSE
> any(is.na(allData$positions))
[1] FALSE
```

Okay so let's try to run snp_autoSVD on a subset of the data and see what
happens.

```R
> subset_G_FBM <- add_code256(big_copy(t(allData$G[1:10000, 1:90]),
                                        type="raw"), code = bigsnpr:::CODE_012)
> svd0 <- snp_autoSVD(subset_G_FBM, 
                      infos.chr = allData$chromosome[1:10000],
                      infos.pos = allData$positions[1:10000],
                      ncores=4)
Phase of clumping at r2 > 0.2.. keep 1311 SNPs.

Iteration 1:
Computing SVD..

Converged!
```

It works on a subset but that subset didn't have any LD regions. What if we
increase the size of the subset?

```R
> subset_G_FBM <- add_code256(big_copy(t(allData$G[1:105000, 1:90]),
                                        type="raw"), code = bigsnpr:::CODE_012)
> svd0 <- snp_autoSVD(subset_G_FBM,
                      infos.chr = allData$chromosome[1:105000],
                      infos.pos = allData$positions[1:100000],
                      ncores=4)
Phase of clumping at r2 > 0.2.. keep 15188 SNPs.

Iteration 1:
Computing SVD..

Converged!
```

The function still works but still no LD regions. However, I may have found at
least one spot that is causing some issues.

```R
> subset_G_FBM <- add_code256(big_copy(t(allData$G[1:106000, 1:90]),
                                        type="raw"), code = bigsnpr:::CODE_012)
> svd0 <- snp_autoSVD(subset_G_FBM,
                      infos.chr = allData$chromosome[1:106000],
                      infos.pos = allData$positions[1:106000],
                      ncores=4)
Phase of clumping at r2 > 0.2.. keep 15430 SNPs.

Iteration 1:
Computing SVD..
Error in na.fail.default(data) : missing values in object
```

Somewhere between the positions 105k and 106k is causing the issue here. What 
about other spots?

```R
subset_G_FBM <- add_code256(big_copy(t(allData$G[200000:808296, 1:90]),
                                     type="raw"), code = bigsnpr:::CODE_012)
svd0 <- snp_autoSVD(subset_G_FBM,
                    infos.chr = allData$chromosome[200000:808296],
                    infos.pos = allData$positions[200000:808296],
                    ncores=4)
Phase of clumping at r2 > 0.2.. keep 98376 SNPs.

Iteration 1:
Computing SVD..
1 long-range LD regions were detected..

Iteration 2:
Computing SVD..

Converged!
```

Works here and there is even a LD region that was detected but...

```R
> subset_G_FBM <- add_code256(big_copy(t(allData$G[200000:808300, 1:90]),
                                        type="raw"), code = bigsnpr:::CODE_012)
> svd0 <- snp_autoSVD(subset_G_FBM,
                      infos.chr = allData$chromosome[200000:808300],
                      infos.pos = allData$positions[200000:808300],
                      ncores=4)
Phase of clumping at r2 > 0.2.. keep 98377 SNPs.

Iteration 1:
Computing SVD..
Error in na.fail.default(data) : missing values in object
```

So here somewhere between 808295 and 808300 is what is causing the issue in
this subset but when we look at the data...

```R
> G_FBM[,808295:808300]
      [,1] [,2] [,3] [,4] [,5] [,6]
 [1,]    0    1    0    0    2    0
 [2,]    0    1    0    1    0    0
 [3,]    1    1    0    0    0    0
 [4,]    2    1    0    0    0    0
 [5,]    2    1    0    0    1    0
 [6,]    1    1    0    0    2    0
 [7,]    0    1    0    0    0    0
 [8,]    1    1    0    0    1    1
 [9,]    1    1    0    0    1    0
[10,]    1    1    0    0    1    0
[11,]    1    1    0    0    1    0
[12,]    2    1    0    0    2    0
[13,]    2    1    0    0    2    0
[14,]    0    1    0    0    0    0
[15,]    2    1    0    0    1    1
[16,]    2    1    0    0    2    0
[17,]    1    1    0    0    1    0
[18,]    0    1    0    0    0    0
[19,]    0    1    0    0    0    0
[20,]    1    1    0    0    1    0
[21,]    0    1    0    0    0    0
[22,]    0    1    0    0    0    0
[23,]    2    1    0    0    2    0
[24,]    1    1    0    0    2    1
[25,]    0    1    0    0    0    0
[26,]    0    1    0    0    0    0
[27,]    0    1    0    0    0    0
[28,]    0    1    0    0    0    0
[29,]    1    1    0    0    1    0
[30,]    0    1    0    0    0    0
[31,]    0    1    0    1    0    0
[32,]    0    1    0    0    1    0
[33,]    0    1    0    0    0    0
[34,]    0    1    0    0    0    0
[35,]    0    1    0    0    1    1
[36,]    1    1    0    0    0    0
[37,]    1    1    0    0    0    0
[38,]    1    1    0    0    0    0
[39,]    1    1    0    0    0    0
[40,]    1    1    0    0    1    0
[41,]    0    1    0    0    0    0
[42,]    0    1    0    0    0    0
[43,]    1    1    0    1    1    0
[44,]    1    1    0    0    1    0
[45,]    1    1    0    0    1    0
[46,]    0    1    2    0    0    0
[47,]    0    1    2    0    0    0
[48,]    0    1    2    0    0    0
[49,]    0    1    1    1    0    0
[50,]    0    1    2    0    0    0
[51,]    2    1    0    0    0    0
[52,]    2    1    0    0    1    0
[53,]    1    1    0    0    0    0
[54,]    1    1    0    0    0    0
[55,]    2    1    0    0    2    0
[56,]    2    1    0    0    0    0
[57,]    1    1    0    0    0    0
[58,]    1    1    0    0    1    1
[59,]    0    1    0    0    0    0
[60,]    1    1    0    0    1    1
[61,]    1    1    0    0    0    0
[62,]    0    1    0    0    0    0
[63,]    1    1    0    0    1    1
[64,]    1    1    0    0    0    0
[65,]    1    1    0    0    1    1
[66,]    1    1    0    0    1    1
[67,]    2    1    0    0    2    0
[68,]    2    1    0    0    2    0
[69,]    2    1    0    0    2    0
[70,]    1    1    0    0    1    0
[71,]    1    1    0    0    0    0
[72,]    0    1    0    0    0    0
[73,]    0    1    0    1    0    0
[74,]    1    1    0    1    0    0
[75,]    0    1    0    0    1    0
[76,]    0    1    0    0    0    0
[77,]    1    1    0    1    0    0
[78,]    1    1    0    0    1    0
[79,]    0    1    0    0    0    0
[80,]    1    1    0    0    0    0
[81,]    1    1    0    0    1    0
[82,]    2    1    0    0    1    0
[83,]    2    1    0    0    0    0
[84,]    1    1    0    0    1    0
[85,]    0    1    0    1    0    0
[86,]    0    1    0    1    1    0
[87,]    0    1    0    0    0    0
[88,]    1    1    0    0    2    0
[89,]    0    1    0    1    1    0
[90,]    0    1    0    0    1    0
```

You can see here that each individual at SNP 808296 is a heterozygote. We need
to filter for these. We can find them using ```big_counts()``` and then filter
them.

```R
# Search all of the SNP sites where all individuals are heterozygotes (have a 
# value of 1)
> heteros <- which(big_counts(G_FBM)[2,] == 90)
 [1]  105422  153456  153457  153458  153459  808296 1943627 2459357 2470340
[10] 2499552 2499553
```

Now let's see how it runs without these sites...

```R
> subset_G_FBM <- add_code256(big_copy(t(allData$G[-heteros, 1:90]),
                                        type="raw"), code = bigsnpr:::CODE_012)
> svd0 <- snp_autoSVD(subset_G_FBM,
                      infos.chr = allData$chromosome[-heteros],
                      infos.pos = allData$positions[-heteros],
                      k = 13,
                      ncores=4)
Phase of clumping at r2 > 0.2.. keep 480406 SNPs.

Iteration 1:
Computing SVD..
11 long-range LD regions were detected..

Iteration 2:
Computing SVD..
0 long-range LD regions were detected..

Iteration 3:
Computing SVD..
0 long-range LD regions were detected..

Iteration 4:
Computing SVD..

Converged!
```

IT WORKED! So, let's see some plots

## Output Thinned SNP Matrix for Future Analysis and a File That Lists and Identifies Them

I will write the locations of the SNPs that remained in the character matrix
after removing SNPs where all individuals were heterozygotes and SNPs that were
pruned. I will write to another file all of the locations of the SNPs that were
removed (as described above) into another file.

```R
write(paste(allData$chromosome[attr(svd0, "subset")],
            allData$positions[attr(svd0, "subset")], sep='_'),
      "allLociLocationsAfterThinning13PCs.txt")

all_heteros_prunedSNPs <- which(!c(1:ncol(G_FBM)) %in% attr(svd0, "subset"))
write(all_heterozygotes_prunedSNPs, "allHeterozygoteAndPrunedSNPs.txt")
```

## Conduct Principal Components of population structure (was also going to compare to mitochondria)

```R
# Scree plot
> plot(svd0$d)
# Let's increase the amount of PCs and see what happens to this plot
> svd20 <- snp_autoSVD(subset_G_FBM,
                      infos.chr = allData$chromosome[-heteros],
                      infos.pos = allData$positions[-heteros],
                      k = 20,
                      ncores = 4)
Phase of clumping at r2 > 0.2.. keep 480406 SNPs.

Iteration 1:
Computing SVD..
16 long-range LD regions were detected..

Iteration 2:
Computing SVD..
0 long-range LD regions were detected..

Iteration 3:
Computing SVD..

Converged!
#Hmm more LD regions were detected...
> plot(svd20$d)
> summary(svd20)
# From this I would say that we need to make the cutoff for PCs at 13. Correct
# me if this isn't right...
> svd13 <- snp_autoSVD(subset_G_FBM,
                      infos.chr = allData$chromosome[-heteros],
                      infos.pos = allData$positions[-heteros],
                      k = 13,
                      ncores = 4)
Phase of clumping at r2 > 0.2.. keep 480406 SNPs.

Iteration 1:
Computing SVD..
11 long-range LD regions were detected..

Iteration 2:
Computing SVD..
0 long-range LD regions were detected..

Iteration 3:
Computing SVD..
0 long-range LD regions were detected..

Iteration 4:
Computing SVD..

Converged!

> summary(svd13)
> plot(svd13)
> dev.copy(png, "screePlot13PCs.png")
> dev.off()
```

Now let's see some plots of the principle components

```R
> plot(svd0, type = "scores") +
  aes(shape = fbm.allData$fam$pop.id, color = fbm.allData$fam$pop.id) +
  geom_point(size=4) +
  scale_shape_manual(values=1:nlevels(fbm.allData$fam$pop.id)) +
  scale_color_manual(values = unique(fbm.allData$fam$color)) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "13PCsPC1-2.png")
> dev.off()

> plot(svd0, type = "scores", scores=3:4) +
  aes(shape = fbm.allData$fam$pop.id, color = fbm.allData$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.allData$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "13PCsPC3-4.png")
> dev.off()

> plot(svd0, type = "scores", scores=c(1, 3)) +
  aes(shape = fbm.allData$fam$pop.id, color = fbm.allData$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.allData$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "13PCsPC1-3.png")
> dev.off()

> plot(svd0, type = "scores", scores=c(1, 4)) +
  aes(shape = fbm.allData$fam$pop.id, color = fbm.allData$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.allData$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "13PCsPC1-4.png")
> dev.off()

> plot(svd0, type = "scores", scores=c(2, 3)) +
  aes(shape = fbm.allData$fam$pop.id, color = fbm.allData$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.allData$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "13PCsPC2-3.png")
> dev.off()

> plot(svd0, type = "scores", scores=c(2, 4)) +
  aes(shape = fbm.allData$fam$pop.id, color = fbm.allData$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.allData$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "13PCsPC2-4.png")
> dev.off()

> plot(svd0, type = "loadings", loadings = 1:4, coeff = 0.6)
```

Set of random snps (10000, 50000). I will do set.seed(233) so that the same samples may be found when doing this same analysis later.

```R
# Load the data and metadata...
> allData <- readRDS('matrixAndMetadata.rds')
> metadata <- read.csv('SampleMetaData.csv')
# I need to first remove the sites where all individuals are heterozygotes so
# they are not included in the random sample.
> G_FBM <- add_code256(big_copy(t(allData$G), type="raw"),
                                code = bigsnpr:::CODE012)
> heteros <- which(big_counts(G_FBM)[2,] == 90)
> noheteros_allData <- list(G = allData$G[-heteros, ],
                            positions = allData$positions[-heteros],
                            chromosomes = allData$positions[-heteros])
# Set seed and generate 10K random samples of loci to include in the analysis 
# and subset the dataset.
> set.seed(10)
> random10 <- sample(1:nrow(noheteros_allData$G), 10000)
> random10_allData <- list(G = noheteros_allData$G[random10, 1:90],
                           positions = noheteros_allData$positions[random10],
                           chromosomes = noheteros_allData$positions[random10])
# Get the subsetted data into the correct format.
> subset10_G_FBM <- add_code256(big_copy(t(random10_allData$G), type="raw"),
                                           code = bigsnpr:::CODE_012)
> fbm.random10 <- list(genotype = subset10_G_FBM,
                       fam = data.frame(sample.id=metadata$Sample.ID,
                                        pop.id=metadata$Pop.ID,
                                        sex = metadata$Sex),
                       map = data.frame(positions = random10_allData$positions,
                                        chromosome = random10_allData$chromosome))
# Perform the clumping and check for long-range LD regions
> subset10_svd13 <- snp_autoSVD(subset10_G_FBM,
                                infos.chr = random10_allData$chromosomes,
                                infos.pos = random10_allData$positions,
                                k = 13,
                                ncores = 4)
Phase of clumping at r2 > 0.2.. keep 10000 SNPs.

Iteration 1:
Computing SVD..
20 long-range LD regions were detected..

Iteration 2:
Computing SVD..
8 long-range LD regions were detected..

Iteration 3:
Computing SVD..

Converged!
> plot(subset10_svd13)
> subset10_svd4 <- snp_autoSVD(subset10_G_FBM,
                               infos.chr = random10_allData$chromosomes,
                               infos.pos = random10_allData$positions,
                               k = 4,
                               ncores = 4)
Phase of clumping at r2 > 0.2.. keep 10000 SNPs.

Iteration 1:
Computing SVD..
11 long-range LD regions were detected..

Iteration 2:
Computing SVD..
3 long-range LD regions were detected..

Iteration 3:
Computing SVD..
2 long-range LD regions were detected..

Iteration 4:
Computing SVD..

Converged!
```

Now we can take a look at some of the plots

```R
> plot(subset10_svd4, type = "scores") +
    aes(shape = fbm.random10$fam$pop.id, color = fbm.random10$fam$pop.id) +
    scale_shape_manual(values = 1:nlevels(fbm.random10$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset10_svd4_PC1-2.png")
> dev.off()

> plot(subset10_svd4, type = "scores", scores=3:4) +
    aes(shape = fbm.random10$fam$pop.id, color = fbm.random10$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random10$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset10_svd4_PC3-4.png")
> dev.off()

> plot(subset10_svd4, type = "scores", scores=c(1, 3)) +
    aes(shape = fbm.random10$fam$pop.id, color = fbm.random10$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random10$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset10_svd4_PC1-3.png")
> dev.off()

> plot(subset10_svd4, type = "scores", scores=c(1, 4)) +
  aes(shape = fbm.random10$fam$pop.id, color = fbm.random10$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.random10$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "subset10_svd4_PC1-4.png")
> dev.off()

> plot(subset10_svd4, type = "scores", scores=c(2, 3)) +
    aes(shape = fbm.random10$fam$pop.id, color = fbm.random10$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random10$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset10_svd4_PC2-3.png")
> dev.off()

> plot(subset10_svd4, type = "scores", scores=c(2, 4)) +
    aes(shape = fbm.random10$fam$pop.id, color = fbm.random10$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random10$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset10_svd4_PC2-4.png")
> dev.off()

```

Now let's do the same for 50K random samples.

```R
> set.seed(50)
> random50 <- sample(1:nrow(noheteros_allData$G), 50000)
> random50_allData <- list(G = noheteros_allData$G[random50, 1:90],
                           positions = noheteros_allData$positions[random50],
                           chromosomes = noheteros_allData$positions[random50])
> subset50_G_FBM <- add_code256(big_copy(t(random50_allData$G), type="raw"),
                                           code = bigsnpr:::CODE_012)
> fbm.random50 <- list(genotype = subset50_G_FBM,
                       fam = data.frame(sample.id=metadata$Sample.ID,
                                        pop.id=metadata$Pop.ID,
                                        sex = metadata$Sex),
                       map = data.frame(positions = random50_allData$positions,
                                        chromosome = random50_allData$chromosome))
> subset50_svd13 <- snp_autoSVD(subset50_G_FBM,
                                infos.chr = random50_allData$chromosomes,
                                infos.pos = random50_allData$positions,
                                k = 13,
                                ncores = 4)
Phase of clumping at r2 > 0.2.. keep 50000 SNPs.

Iteration 1:
Computing SVD..
100 long-range LD regions were detected..

Iteration 2:
Computing SVD..
30 long-range LD regions were detected..

Iteration 3:
Computing SVD..
5 long-range LD regions were detected..

Iteration 4:
Computing SVD..
1 long-range LD regions were detected..

Iteration 5:
Computing SVD..

Converged!
# Take a look at the scree plot
> plot(subset50_svd13)
# Let's do the same analysis but set the cutoff at 4 PCs
> subset50_svd4 <- snp_autoSVD(subset50_G_FBM,
                               infos.chr = random50_allData$chromosomes,
                               infos.pos = random50_allData$positions,
                               k = 4,
                               ncores = 4)
Phase of clumping at r2 > 0.2.. keep 50000 SNPs.

Iteration 1:
Computing SVD..
58 long-range LD regions were detected..

Iteration 2:
Computing SVD..
21 long-range LD regions were detected..

Iteration 3:
Computing SVD..
0 long-range LD regions were detected..

Iteration 4:
Computing SVD..

Converged!

```

Plot the results and save your plots

```R
> plot(subset50_svd4, type = "scores") +
    aes(shape = fbm.random50$fam$pop.id, color = fbm.random50$fam$pop.id) +
    scale_shape_manual(values = 1:nlevels(fbm.random50$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset50_svd4_PC1-2.png")
> dev.off()

> plot(subset50_svd4, type = "scores", scores=3:4) +
    aes(shape = fbm.random50$fam$pop.id, color = fbm.random50$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random50$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset50_svd4_PC3-4.png")
> dev.off()

> plot(subset50_svd4, type = "scores", scores=c(1, 3)) +
    aes(shape = fbm.random50$fam$pop.id, color = fbm.random50$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random50$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset50_svd4_PC1-3.png")
> dev.off()

> plot(subset50_svd4, type = "scores", scores=c(1, 4)) +
  aes(shape = fbm.random50$fam$pop.id, color = fbm.random50$fam$pop.id) +
  scale_shape_manual(values=1:nlevels(fbm.random50$fam$pop.id)) +
  geom_point(size=4) +
  labs(shape = "Population", color = "Population")
> dev.copy(png, "subset50_svd4_PC1-4.png")
> dev.off()

> plot(subset50_svd4, type = "scores", scores=c(2, 3)) +
    aes(shape = fbm.random50$fam$pop.id, color = fbm.random50$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random50$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset50_svd4_PC2-3.png")
> dev.off()

> plot(subset50_svd4, type = "scores", scores=c(2, 4)) +
    aes(shape = fbm.random50$fam$pop.id, color = fbm.random50$fam$pop.id) +
    scale_shape_manual(values=1:nlevels(fbm.random50$fam$pop.id)) +
    geom_point(size=4) +
    labs(shape = "Population", color = "Population")
> dev.copy(png, "subset50_svd4_PC2-4.png")
> dev.off()
```

## Use thinned SNPs to calculate null distribution for OutFLANK and PCAdapt,
## then conduct outlier tests on all SNPs. (If you are planning any other 
## methods, let’s make sure to use the same set of thinned SNPs)

## Once I get environmental data, will also calculate SNP-environment associations for methods that aren’t sensitive to recombination rate variation (see manuscript I shared with you)
