# 2019-05-01 - KBW - Running OutFLANK on Wild and Selection Subsets

## Neutral Parameters Not Fitting a Chi-Distribution

After running OutFLANK on the Wild and Selection subset, I realized the neutral parameters were not fitting the chi distribution that is necessary for the analysis. This is likely due to the removal of individuals from populations that do not fall into the Wild or Selection subsets in the thinned SNP set.

## Solution

In order to resolve this issue, we decided to run the thinning steps again after the removal of populations outside the subset.