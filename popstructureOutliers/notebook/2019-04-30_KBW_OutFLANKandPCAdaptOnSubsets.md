# 2019-05-01 - KBW - Running OutFLANK on Wild and Selection Subsets

## Neutral Parameters Not Fitting a Chi-Distribution

After running OutFLANK on the Wild and Selection subset, I realized the neutral parameters were not fitting the chi distribution that is necessary for the analysis. This is likely due to the removal of individuals from populations that do not fall into the Wild or Selection subsets from the thinned SNP set.

## Solution

In order to resolve this issue, we decided to run the thinning steps again after the removal of populations outside the respective subset.

### Thinning After Removal of Populations

Running the thinning procedures to get a thinned SNP set for the Wild and Selection subsets does not work as well as one might think. Doing so eliminates nearly all SNPs from the dataset only leaving < 10K of the original ~4 million and also does not provide a null distribution that is any better than the fit when just removing populations outside the subsets from the thinned SNP set for unrelated subset. I tried using a base pair window size of 5K and 50K during the thinning steps and did not change this dramatic decrease.

<p float="left">
    <img src="../data/large_outputs/outlierAnalysis/outlierAnalysis_atlantic_wild_subset/plots/neutralParams/outflankResults_atlantic_wild_subset.png" width="400">
    <img src="../data/large_outputs/outlierAnalysis/outlierAnalysis_atlantic_selection_subset/plots/neutralParams/outflankResults_atlantic_selection_subset.png" width="400">
</p>
<div>
    <b>Null Distributions Established Using Thinned SNPs of Specific subsets</b>
    <p>Right: Atlantic Wild</p>
    <p>Left: Atlantic Selection</p>
</div>

### Going Forward

I would suggest that we either just deal with the mediocre fit of the null distribution that we had before. If this is not acceptable, we may not be able to use OutFLANK or PCAdapt.