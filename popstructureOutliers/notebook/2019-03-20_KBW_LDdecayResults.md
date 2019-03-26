# 2019-03-21 - KBW - LD Decay Results

## LD Decay Analysis

In order to make decisions throughout the thinning process, we needed to revisit the LD decay analysis.

## Methods

### Calculating and summarizing LD for each chromosome

We used `vcftools` for the LD calculations. Across all LD decay analyses, we used eight different base pair window sizes from 50 b.p. to 500K b.p. (50, 200-250, 450-500, 4500-5000, 9500-10000, 49500-50000, 99500-100000, 499500-500000). An example shell command can be seen here:

```shell
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --geno-r2 --ld-window-bp-min 200 --ld-window-bp 250 --out ldAnalysisData/geno_ld_window_200-250
```

Lastly, to summarize and plot LD for each chromosome, I took a mean of the LD values for each chromosome. Each point that you see on the plots below represent that mean LD value. To see the variance in LD at each window size for each chromosome, error bars are shown which one standard deviation from the mean.

### Individuals and Subsets of Individuals

We discussed removing individuals from the Laguna Madre (LM) population as they are highly (some might say 'abnormally') related with one another and most likely affecting all other analyses. Therefore, to start, we ran two different LD decay analyses (1) including all populations and (2) all populations excluding LM samples.

We also wanted to see how LD decay is affected when only considering specific subsets of the data as well. We analyzed two specific subsets, the first being (1) Atlantic wild populations (Maine_SM, DelBay_HC, DelBay_CS, ChesBay_CLP, ChesBay_HC_VA) and individuals from (2) Atlantic selection lines (UMFS, NEH, DEB, LOLA).

#### Related Individuals

One caveat with the subset of individuals from Atlantic selection lines is that there were individuals that were highly related (relatedness > 0.5). These individuals (NEH_1 & LOLA_2) were removed from the LD decay analysis of Atlantic selection lines.

## Results

### All populations and Excluding Laguna Madre

The effect that they have on the LD decay is fairly obvious. Here is the LD decay when LM samples are included...

![](https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/popstructureOutliers/figures/1LD_analysis/ldDecayPlot_allpops.png)
__../figures/1LD_analysis/LdDecayPlot_allpops.png__

Then, compare that to same analysis when LM samples are excluded...

![](https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/popstructureOutliers/figures/1LD_analysis/ldDecayPlot_excluding_LM.png)
__../figures/1LD_analysis/LdDecayPlot_excluding_LM.png__

The LD decay is much more pronounced when excluding LM individuals. This will lead us to using a much smaller window when performing the SNP thinning. However, we wanted to see how LD decay was behaving when considering certain subsets of the populations as well.

### Atlantic Wild and Selection subsets

Here are the plots from those analyses, respectively:

![](https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/popstructureOutliers/figures/1LD_analysis/ldDecayPlot_excluding_selection.png)
__../figures/1LD_analysis/LdDecayPlot_excluding_selection.png__

![](https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/popstructureOutliers/figures/1LD_analysis/ldDecayPlot_excluding_wild.png)
__../figures/1LD_analysis/LdDecayPlot_excluding_wild.png__

The comparison between the wild and selection plots might be what you would expect given the relatedness amongst the selection lines. The LD decay does not go quite as low with selection lines

One thing that sticks out immediately to me is the change in LD decay on chromosome NC_035785.1 from when only excluding LM and when considering only the wild populations.

## Follow up

Using these analyses, we need to decide upon a window size for the thinning steps. When thinning initially, we decided on a window size where the LD decay had reached ~0.10.

Now that we are removing LM, we might consider using a window around 5K or 10K. However, we could be conservative and use a window of 50K.

Let me know what you think!
