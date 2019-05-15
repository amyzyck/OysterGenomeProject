# 2019-05-09 - KBW - Combined Outlier Plots

I have combined the OutFLANK and PCAdapt outlier plots into a 2 x 1 matrix plot for each chromosome. Each plot was then combined into a single pdf file. I have uploaded the individual plots for each chromosome and the pdf file to github [here](https://github.com/jpuritz/OysterGenomeProject/tree/master/popstructureOutliers/figures/4outlier).

## PCAdapt Outliers and Transformation of Negative Log10 P-values

Through following the PCAdapt [documentation](https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html) , I realized that I had not written the code in `outlierAnalysis.R` to include a column in the final data that marks SNPs as outliers using PCAdapt. Also, the final data that we have right now has the negative log10 p-value rather than the raw p-values. When finding outliers using PCAdapt you need the raw p-values. For the time being and to create the current plots I have just transformed the negative log10 p-values back using the inverse: `10^-(negativeLog10p)`.

When finding the outliers, I used the Benjamini-Hochberg procedure which is a moderately conservative method. This can be changed easily if a different procedure is desirable.

## PCAdapt vs. OutFLANK Outliers

There is a sizeable difference between the amount of outliers that are flagged by the outlier detection methods. This is readily apparent in the plots using the exclude_LM subset. I am interested in the difference in these outliers. I will look into plotting the outliers from both detection methods onto the same plot with different colored points denoting outliers flagged by PCAdapt, OutFLANK, and those flagged by both. I feel it would be interesting to see the discrepancies between both methods as well as those SNPs that both methods agree are outliers.
