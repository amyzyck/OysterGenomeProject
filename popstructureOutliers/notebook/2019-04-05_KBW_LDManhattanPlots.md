# 2019-04-05 - KBW - LD Manhattan Plots

I have changed the script for plotting the LD manhattan plots to (1) run as a sequence of functions (to be able to plot different LD window sizes and subsets easier), (2) visualize the breakpoints and misassemblies, (3) visualize the thinned SNPs on the chromosomes where there are breakpoints.

1. I am lazy and want the ability to run different window sizes and subsets quickly. You can see the new code to generate the plots at `src/1LD_analysis/ld_on_manhtn.R`
2. There are misassemblies that can be visualized on chromosomes NC_035784.1, NC_035785.1, and NC_035788.1. The plots of these chromosomes have a semi-transparent, green rectangle where these are located.
3. You will also see some semi-transparent grey vertical lines on the plots with the misassemblies. Each line denotes the location of a thinned SNP that is on that chromosome. We did this in hopes of visualizing the proportion of thinned SNPs in areas of the misassemblies and areas we are more confident about.

You can take a look at the plots [here](https://github.com/jpuritz/OysterGenomeProject/tree/master/popstructureOutliers/figures/1LD_analysis)