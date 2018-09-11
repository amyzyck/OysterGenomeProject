### We’re going to work with

* SNP.TRSdp5g95FnDNAmaf05.vcf.gz
  * This file has been filtered for MAF > 0.05, and missing data < 5%.
* SNP.TRSdp5g75mtDNA.recode.vcf (mitochondria)
  * Filtered SNPs from mitochondrial genome
* INDELs.TRSdp5g75FnDNA.vcf.gz (indels)
  * Filtered INDel variants with missing data of less than 25% allowed

### Our plan:
1) thin SNPs for LD based on recommendation of (Privé et al. 2018) using bigsnpR package.
2) Conduct Principal Components of population structure (was also going to compare to mitochondria)
3) Use thinned SNPs to calculate null distribution for OutFLANK and PCAdapt, then conduct outlier tests on all SNPs. (If you are planning any other methods, let’s make sure to use the same set of thinned SNPs)
4) Once I get environmental data, will also calculate SNP-environment associations for methods that aren’t sensitive to recombination rate variation (see manuscript I shared with you)
