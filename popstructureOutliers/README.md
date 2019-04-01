### We’re going to work with

* SNP VCF File: __Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz__
  * This file can be found on Jon Puritz's KITT server.
  * This filtered SNPs with a minor allele frequency> 0.05 and no missing data, 2 alleles  The bioinformatics page has a graph of the missingness per individual (before MAF and missingness filtering) https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome_files/figure-markdown_github/unnamed-chunk-14-1.svg
Also of note is that this VCF contains all the individuals, including the aquaculture and inbred lines.
* SNP.TRSdp5g75mtDNA.recode.vcf (mitochondria)
  * Filtered SNPs from mitochondrial genome
* INDELs.TRSdp5g75FnDNA.vcf.gz (indels)
  * Filtered INDel variants with missing data of less than 25% allowed

### Our plan:
* LD analysis
  * on original vcf file
  * Use vcftools to calculate LD using several different window sizes
  * popstructureOutliers/1LD_analysis (check folder)
  * LD decay
    * use this to help determine the window size for SNP thinning 
    * 25, 50, 500, 1000, 2500, 5000, 10000, 30000, 50000 (for bigger sizes 100 bp window)
  * LD Manhattan plot
    * 4600-5000
    * around 30000
* SNP thinning for population structure
  * folder: popstructureOutliers/2thinning/
  * These steps were done by JP after the full genotyping was completed:
    * Filtering in vcf tools
      * remove triallelic sites
      * filter for sites that have 0% missing data
    * The VCF file that JP provided after doing these filtering steps is __Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz__
  * The following steps are done using the populationStructureScipt.R:
    * Thin SNPs for LD based on recommendation of (Privé et al. 2018) using bigsnpr::snp_autoSVD()
      * Uses clumping instead of pruning, and thin for long-range LD not captured by clumping
    * Thin on a window size where the correlation among SNPs < 0.1 or maybe 0.05
      * based on the LD decay analysis that window is a bp window of 100K.
  * Output thinned SNP matrix for future analysis and a file that lists which SNPs they are
    * The thinned SNP matrix can be found in __popstructureOutliers/data/thinned_data/__ folder
  * In addition, want a set of ~ 50,000 SNPs and 20,000 SNPs to compare for population structure that are a subset of the initial thinned set. These may be used for SNP chip.
    * The randome subsets of 50K and 20K can also be found in the __popstructureOutliers/data/thinned_data/__
  
* population structure
  * Principal Components output from snp_autoSVD function
  * compare and mak
  
* Genome scans
  * Use thinned SNPs to calculate null distribution for OutFLANK and PCAdapt, then conduct outlier tests on all SNPs. (If you are planning any other methods, let’s make sure to use the same set of thinned SNPs)
    * Thinned SNP sets are located in the __popstructureOutliers/data/thinned_snps/__ folder
    * The set of thinned SNPs should be the set located in the __thinnedMatrixAndMetaData1e+05Window.rds__ file
      * There might be multiple files with this naming scheme given that we have run some analyses with some populations excluded but the file name should indicate this.
  * Once I get environmental data, will also calculate SNP-environment associations for methods that aren’t sensitive to recombination rate variation (see manuscript I shared with you)


### Triallelic sites
* We have yet to determine what we will do with the triallelic sites. For now, they have been filtered out of the data for analysis.
* Any suggestions on how to deal with these sites would be appreciated.

# Files on KITT server and a description of their contents

* __Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz__ - Copy of file in the VCF_files directory.
* __exome__ - Folder containing ...
* __genotypeMatrix.rds__ - R object containing the conversion of Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz to a genotype matrix, chromosome labels positions of all loci minus the fixed heterozygote sites.
* __genotypeMatrixAndMetadata.rds__ - R object containing the data from genotypeMatrix.rds along with metadata from modifiedColors_samplemetadata.csv
* __The following files contain LD calculations with using different physical distance windows which can be seen in the file name. These were created using ldAnalysis.sh__
    1. __geno_ld_window_200-250.geno.ld__
    2. __geno_ld_window_4500-5000.geno.ld__
    3. __geno_ld_window_450-500.geno.ld__
    4. __geno_ld_window_49500-50000.geno.l__
    5. __geno_ld_window_50-50.geno.ld__
    6. __geno_ld_window_9500-10000.geno.ld__
    7. __geno_ld_window_99500-100000.geno.ld__
    8. __geno_ld_window_499500-500000.geno.ld__
* __INDELs.TRSdp5g75FnDNA.vcf.gz__ - Copy of file in the VCF_files directory.
* __ldAnalysis.sh__ - Bash script that will get the LD calculations needed for the LD analyses
* __modified_samplemetadata.csv__ - metadata that has been modified in order to include the order in which the samples are found in the VCF file
* __populationStructureScript.R__ - Script to perform thinning and getting 50K and 20K
* __README.md__ - This file : D
* __SampleMetaData.csv__ - SampleMetaData.md file found on github repo converted
  to .csv format.
  
  
  ## Epigenetics
HP: I am including here a link to the CpG o/e file, as well as the description of the file generation.

http://gannet.fish.washington.edu/Atumefaciens/20190225_cpg_oe/ID_CpG_labelled_all

https://github.com/hputnam/EastOyEpi/blob/master/methods.md
