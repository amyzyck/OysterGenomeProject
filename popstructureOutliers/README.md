### We’re going to work with

* SNP.TRSdp5g95FnDNAmaf05.vcf.gz
  * This file has been filtered for MAF > 0.05, and missing data < 5%. The 5% missing data is for a locus.  For this data set, I didn't filter any individuals for missing data because our sample size was so small already.  The bioinformatics page has a graph of the missingness per individual (before MAF and missingness filtering) https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome_files/figure-markdown_github/unnamed-chunk-14-1.svg
Also of note is that this VCF contains all the individuals, including the aquaculture and inbred lines.
* SNP.TRSdp5g75mtDNA.recode.vcf (mitochondria)
  * Filtered SNPs from mitochondrial genome
* INDELs.TRSdp5g75FnDNA.vcf.gz (indels)
  * Filtered INDel variants with missing data of less than 25% allowed

### Our plan:
* remove triallelic sites
* filter for sites that have 0% missing 
* thin SNPs for LD based on recommendation of (Privé et al. 2018) using bigsnpR package.
* output thinned SNP matrix for future analysis and a file that lists which SNPs they are
* Conduct Principal Components of population structure (was also going to compare to mitochondria)
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
