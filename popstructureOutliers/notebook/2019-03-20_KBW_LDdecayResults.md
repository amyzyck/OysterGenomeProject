# 2019-03-21 - KBW - LD Decay Results

## LD Decay Analysis

In order to make decisions throughout the thinning process, we needed to revisit the LD decay analysis. We discussed removing individuals from the Laguna Madre (LM) population as they are highly (some might say 'abnormally') related with one another and most likely affecting all other analyses.

The affect that they have on the LD decay is fairly obvious. Here is the LD decay when LM samples are included...

![](../figures/1LD_analysis/LdDecayPlot_allpops.png)
__../figures/1LD_analysis/LdDecayPlot_allpops.png__

Then, compare that to same analysis when LM samples are excluded...

![](../figures/1LD_analysis/LdDecayPlot_excluding_LM.png)
__../figures/1LD_analysis/LdDecayPlot_excluding_LM.png__

The LD decay is much more pronounced when excluding LM individuals. This will lead us to using a much smaller window when performing the SNP thinning. However, we wanted to see how LD decay was behaving when considering certain subsets of the populations as well.

### Atlantic Wild and Selection subsets

We also wanted to see how LD decay is affected when only considering...
* Atlantic wild populations(Maine_SM, DelBay_HC, DelBay_CS, ChesBay_CLP, ChesBay_HC_VA)
* and then, the Atlantic selection lines

Here are the plots from those analyses, respectively:

![](../figures/1LD_analysis/LdDecayPlot_excluding_selection.png)
__../figures/1LD_analysis/LdDecayPlot_excluding_selection.png__

![](../figures/1LD_analysis/LdDecayPlot_excluding_wild.png)
__../figures/1LD_analysis/LdDecayPlot_excluding_wild.png__

The comparison between the wild and selection plots might be what you would expect given the relatedness amongst the selection lines. The LD decay does not go quite as low with selection lines

One thing that sticks out immediately to me is the change in LD decay on chromosome NC_035785.1 from when only excluding LM and when considering only the wild populations.

## Follow up

Using these analyses, we need to decide upon a window size for the thinning steps. When thinning initially, we decided on a window size where the LD decay had reached ~0.10.

Now that we are removing LM, we might consider using a window around 5K or 10K. However, we could be conservative and use a window of 50K. 

Let me know what you think!
