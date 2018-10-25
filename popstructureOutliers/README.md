### We’re going to work with

* SNP.TRSdp5g95FnDNAmaf05.vcf.gz
  * This file has been filtered for MAF > 0.05, and missing data < 5%. The 5% missing data is for a locus.  For this data set, I didn't filter any individuals for missing data because our sample size was so small already.  The bioinformatics page has a graph of the missingness per individual (before MAF and missingness filtering) https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome_files/figure-markdown_github/unnamed-chunk-14-1.svg
Also of note is that this VCF contains all the individuals, including the aquaculture and inbred lines.
* SNP.TRSdp5g75mtDNA.recode.vcf (mitochondria)
  * Filtered SNPs from mitochondrial genome
* INDELs.TRSdp5g75FnDNA.vcf.gz (indels)
  * Filtered INDel variants with missing data of less than 25% allowed

### Our plan:
* LD analysis
  * on original vcf file
  * using vcftools
  * Populationstructure/LD_decay (check folder)
  * LD decay
    * use this to help determine the window size for SNP thinning 
    * 25, 50, 500, 1000, 2500, 5000, 10000, 30000, 50000 (for bigger sizes 100 bp window)
  * LD Manhattan plot
    * 4600-5000
    * around 30000
    
* SNP thinning for population structure
  * Population_structure/SNP_thinning
  * Filtering in vcf tools
    * remove triallelic sites
    * filter for sites that have 0% missing data
    * filter for sites that have all heterozygotes
  * New VCF file name: 
  * thin SNPs for LD based on recommendation of (Privé et al. 2018) using bigsnpR package.
  * use clumping instead of pruning, and thin for long-range LD not captured by clumping
  * thin on a window size where the correlation among SNPs < 0.1 or maybe 0.05
    * based on prelim data that is 5000-XXX SNPs
    * write a revised snp_autoSVD function to accoutn for window size in bp (not size in adjacent SNPs)
  * output thinned SNP matrix for future analysis and a file that lists which SNPs they are
    * in addition, want a set of ~ 50,000 SNPs and 10,000 SNPs to compare for population structure
    that are a subset of the initial thinned set. These may be used for SNP chip.
  
* population structure
  * Principal Components output from snp_autoSVD function
  * compare and mak
  
* Genome scans
  * Use thinned SNPs to calculate null distribution for OutFLANK and PCAdapt, then conduct outlier tests on all SNPs. (If you are planning any other methods, let’s make sure to use the same set of thinned SNPs)
  * Once I get environmental data, will also calculate SNP-environment associations for methods that aren’t sensitive to recombination rate variation (see manuscript I shared with you)


### Triallelic sites
Example of a triallelic site:
```
zgrep "0|2" SNP.TRSdp5g95FnDNAmaf05.vcf.gz | head -n1
```
https://gist.github.com/inutano/f0a2f5c219ab4920c5b5

If a call cannot be made for a sample at a given locus
'.' should be specified for each missing allele in the GT field
e.g. './.' for a diploid genotype and '.' for haploid genotype

### Files in KITT server directory (~/popstructureOutliers) being used in this directory and descriptions of those files.
* The following were used to create vcf files of SNPs after thinning and subsets
  of SNPs:
  * __10KRandomSNPs13PCs.txt__ - Chromosome location and positions (tab separated)
    of 10K randomly selected SNPs from the full dataset of thinned SNPs.
  * __10KThinnedRandomSNPs13PCs.txt__ - Chromosome location and positions (tab
    separated) of a subset of thinned SNPs from the 10K randomly selected SNPs.
  * __50KRandomSNPs13PCs.txt__ - Chromosome location and positions (tab separated)
    of 50K randomly selected SNPs from the full dataset of thinned SNPs.
  * __50KThinnedRandomSNPs13PCs.txt__ - Chromosome location and positions (tab
    separated) of a subset of thinned SNPs from the 10K randomly selected SNPs.
  * __allLociLocationsAfterThinning13PCs_test.txt__ - Chromosome location and
    positions (tab separated) of a subset of thinned SNPs from the 10K randomly
    selected SNPs.
* __INDELs.TRSdp5g75FnDNA.vcf.gz__ - Copy of file in the VCF_files directory.
* __matrixAndMetadata.rds__ - Genotype matrix, positions and chromosome of
  positions generated from the
  SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz vcf file
* __README.md__
* __SampleMetaData.csv__ - SampleMetaData.md file found on github repo converted
  to .csv format.
* __SNP.TRSdp5g95FnDNAmaf05_10KRandomSNPs13PCs.log__ - Output from VCFtools for
  the corresponding file SNP.TRSdp5g95FnDNAmaf05_10KRandomSNPs13PCs.vcf.gz
* __SNP.TRSdp5g95FnDNAmaf05_10KRandomSNPs13PCs.vcf.gz__ - vcf file containing
  the subset of 10K SNPs.
* __SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs.log__ - Output from
  VCFtools for the corresponding file
  SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs.vcf.gz
* __SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs.vcf.gz__ - vcf file
  containing the SNPs which remained after thinning the subset of 10K SNPs.
* __SNP.TRSdp5g95FnDNAmaf05_50KRandomSNPs13PCs.log__ - Output from VCFtools for
  the corresponding file SNP.TRSdp5g95FnDNAmaf05_50KRandomSNPs13PCs.vcf.gz
* __SNP.TRSdp5g95FnDNAmaf05_50KRandomSNPs13PCs.vcf.gz__ - vcf file containing
  the subset of 50K SNPs.
* __SNP.TRSdp5g95FnDNAmaf05_50KThinnedRandomSNPs13PCs.log__ - Output from
  VCFtools for the corresponding file
  SNP.TRSdp5g95FnDNAmaf05_50KThinnedRandomSNPs13PCs.vcf.gz
* __SNP.TRSdp5g95FnDNAmaf05_50KThinnedRandomSNPs13PCs.vcf.gz__ - vcf file
  containing the SNPs which remained after thinning the subset of 50K SNPs.
* __SNP.TRSdp5g95FnDNAmaf05_AfterThinning13PCs.vcf.gz__ - vcf file
  containing the SNPs which remained after thinning the full set of SNPs.
* __SNP.TRSdp5g95FnDNAmaf05allLociLocationsAfterThinning.log__ - Output from
  VCFtools for the corresponding file
  SNP.TRSdp5g95FnDNAmaf05_AfterThinning13PCs.vcf.gz
* __SNP.TRSdp5g95FnDNAmaf05_min0_max2.log__ - Output from VCFtools when
  filtering for out the tri-allelic sites.
* __SNP.TRSdp5g95FnDNAmaf05_min0_max2_noMissing.log__ - Output from VCFtools
  when filtering for out the tri-allelic sites and any sites which had any
  missing data.
* __SNP.TRSdp5g95FnDNAmaf05_min3_max3.log__ - Output from VCFtools when filtering for only the tri-allelic sites.
* __SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2_noMissingData.vcf.gz__ -
  'SNP.TRSdp5g95FnDNAmaf05.vcf.gz' filtered further to exclude any tri-allelic
  sites as well as any sites that had any missing data.
* __SNP.TRSdp5g95FnDNAmaf05_min-allele0_max-allele2.vcf.gz__ -
  'SNP.TRSdp5g95FnDNAmaf05.vcf.gz' filtered further to exclude any tri-allelic
  sites.
* __SNP.TRSdp5g95FnDNAmaf05_min-allele3_max-allele3.vcf.gz__ - vcf file containing only the tri-allelic sites.
* __SNP.TRSdp5g95FnDNAmaf05.vcf.gz__ - Copy of file in the VCF_files directory.
