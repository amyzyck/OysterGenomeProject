CpGOE Environment Association Analysis
================
adowneywall
5/22/2019

# How to Use LFMMCPoe.R

OverviewLFMMCPoe.R is

## Intended Usage

Start in the popstructureOutliers directory.

`Rscript src/6envi_assoc/LFMMwild.R`

This will start the analysis. The script will produce periodic updates
on its progress, and will automatically save plots and results tables to
the appropriate directoris. Here is an example of successful output for
the analysis of association with the “Lat” variable:

    ----------------------------------------------------
    STARTING ANALYSIS FOR: Lat
    ----------------------------------------------------
    [1] "Building lfmm ridge model"
    [1] "Running lfmm test"
    [1] "lfmm.test.ridge$gif" "2.89063197471484"   
    [1] "Calculating spearman's correlation"
    [1] "Plotting LF 1 and 2"
    [1] "Creating data frame"
    [1] "saving dataframe"

## Outputs

1)  An RDS file containing the subset of genotype data being used (only
    wild populations). This is stored in
    “data/genotypeMatrix\_selecting\_Wild.rds”, and I would recommend
    that if you have to restart the analysis or run it multiple times
    you should modfiy the code to read in this object instead of redoing
    the subsetting.
2)  Before LFMM is run, a scree plot is produced in
    “figures/6envi\_assoc” called “PCA\_scree\_plot.png”
3)  (Optional) a Manhattan plot of the p values from the LFMM ridge test
    is created in “figures/6envi\_assoc” called
    "LFMM\_ridge\_0.0\_\[variable name\]\_pvalues\_plot.png"
4)  (Optional) a plot of the populations along the first 2 latent
    factors is created to show clustering. This is saved in
    “figures/6envi\_assoc” and called "LFMM\_ridge\_0.0\_\[variable
    name\]\_LF\_plot.png"
5)  (Optional) a plot of Spearmanns Rho vs -log10 p values for LFMM is
    created in “figures/envi\_assoc” and called
    "Spearmanns\_vs\_LFMM\_\[variable name\]\_plot.png"
6)  A table in “data/envi\_assoc\_results/” called "\[variable
    name\]\_assoc\_results.txt" This table has 5 columns: “Pos”, “Chr”,
    “Unique” (a combination of Chr and Pos), Spearmanns Rho for the
    locus, and LFMM ridge -log10(p) for the locus.
