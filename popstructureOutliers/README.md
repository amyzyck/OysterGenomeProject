### We’re going to work with

* Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz
  * This filtered SNPs with a minor allele frequency> 0.05 and no missing data, 2 alleles  The bioinformatics page has a graph of the missingness per individual (before MAF and missingness filtering) https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome_files/figure-markdown_github/unnamed-chunk-14-1.svg
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
* __ldAnalysis.sh__
* __modifiedColors_samplemetadata.csv__
* __nohup.out__
* __populationStructureScript.R__
* __README.md__
* __SampleMetaData.csv__ - SampleMetaData.md file found on github repo converted
  to .csv format.
* __SNP.TRSdp5g95FnDNAmaf05.vcf.gz__ - Copy of file in the VCF_files directory.
