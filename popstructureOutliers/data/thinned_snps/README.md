# Description of File Names

In this folder are the sets of thinned SNPs that were produced after running the thinning (via pruning) procedures on full genotype matrices. Additionally, subsets of these thinned SNPs are also here and are named in such a way to identify them apart from the full set of thinned SNPs. I have attempted to make it easy for someone to find the file(s) that they are looking for even though they may not be familiar with how they were created. This is done by using a naming scheme across all the files which I describe below.

## Naming scheme

The data within each file should be described by the file name in this order separated using an underscore:

1. Type (or size) of thinned SNPs
   * Full sets: `thinnedMatrixAndMetaData`
   * Subsets of full set of thinned SNPs:
     * `20KRandomSNPs`
     * `50KRandomSNPs`
2. Size of window used for thinning procedure
   * We used a window size of 5000 bp and 50,000 bp. Those values should appear in the file names.
3. Subset of samples used.
   * Several subsets of the total samples were developed for comparisons. These subset names should appear at the end of the file name. The subset name is the same as what should be in the `modified_samplemetadata.csv` file.
