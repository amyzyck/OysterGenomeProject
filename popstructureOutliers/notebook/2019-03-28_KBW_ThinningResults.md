# 2019-03-28 - KBW - Thinning results

After the LD decay analysis, window sizes of 5000 and 50000 were suggested for the thinning steps. I went forward with both to see how thinning performed and the amount of SNPs that were left afterwards.

## Methods

To perform thinning, I used the [bigsnpr package](https://github.com/privefl/bigsnpr/tree/master/R). You can take a look at the code on github [here](https://github.com/jpuritz/OysterGenomeProject/blob/master/popstructureOutliers/src/2thinning/populationStructureScript.R).

### Thinning with subsets

I have made it to where the thinning can be performed while only using specific samples that have been specified using the metadata. The metadata file named `modified_samplemetadata.csv` has been edited to include columns with either a 0 or 1 for specific samples that should be included (1) or excluded (0) from the data prior to any analyses such as the thinning procedures.

The two columns that I used when performing thinning this round was `exclude.LM` and `unrelated`. The first (`exclude.LM`) has 0 for all Laguna Madre individuals (5) and 1 for all other individuals (85). The second (`unrelated`) has 0 for one individual within a pair that were highly related (relatedness > 0.5). This included all LM samples and 9 other individuals from the inbred and selection lines leaving a total of 76 individuals.

## Results

I have pushed the thinned sets of SNPs to github and they can be found [here](https://github.com/jpuritz/OysterGenomeProject/tree/master/popstructureOutliers/data/thinned_snps).

The naming scheme for the full set of thinned SNPs for a particular window and subset is this:
`thinnedMatrixAndMetaData#####Window_name_of_subset_used.rds`

Here is a table of describing the thinned SNPs, windows used, and subsets used.

| File name                                      | Number of SNPs|
|------------------------------------------------|---------------|
| thinnedMatrixAndMetaData5000Window_exclude_LM  |     334011    |
| thinnedMatrixAndMetaData50000Window_exclude_LM |     193255    |
| thinnedMatrixAndMetaData5000Window_unrelated   |     341044    |
| thinnedMatrixAndMetaData50000Window_unrelated  |     184326    |

You will also see that there are subsets from the thinned set of SNPs in the file as well. There should be a subset of 20K and 50K SNPs for each set of thinned SNPs. These subsets include 20K and 50K SNPs that were selected randomly from each chromosome at equal proportions. The naming scheme is similar to the full set of thinned SNPs:
`##KRandomSNPs####Window_name_of_subset_used`
