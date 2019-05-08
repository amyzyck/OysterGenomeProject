# How to Use LFMMwild.R

LFMMwild.R is a rather lengthy script with many inputs and outputs, but it is designed so that it is trivial to adjust what environmental variable the association analysis uses. 

## Intended Usage

Start in the popstructureOutliers directory. 

`Rscript src/6envi_assoc/LFMMwild.R`

This will start the analysis. The script will produce periodic updates on its progress, and will automatically save plots and results tables to the appropriate directoris. Here is an example of successful output for the analysis of association with the "Lat" variable:

```
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
```

## Outputs

1) An RDS file containing the subset of genotype data being used (only wild populations). This is stored in "data/genotypeMatrix_selecting_Wild.rds", and I would recommend that if you have to restart the analysis or run it multiple times you should modfiy the code to read in this object instead of redoing the subsetting.
2) Before LFMM is run, a scree plot is produced in "figures/6envi_assoc" called "PCA_scree_plot.png"
3) (Optional) a Manhattan plot of the p values from the LFMM ridge test is created in "figures/6envi_assoc" called "LFMM_ridge_0.0_[variable name]_pvalues_plot.png"
4) (Optional) a plot of the populations along the first 2 latent factors is created to show clustering. This is saved in "figures/6envi_assoc" and called "LFMM_ridge_0.0_[variable name]_LF_plot.png"
5) (Optional) a plot of Spearmanns Rho vs -log10 p values for LFMM is created in "figures/envi_assoc" and called "Spearmanns_vs_LFMM_[variable name]_plot.png"
6) A table in "data/envi_assoc_results/" called "[variable name]_assoc_results.txt" This table has 5 columns: "Pos", "Chr", "Unique" (a combination of Chr and Pos), Spearmanns Rho for the locus, and LFMM ridge -log10(p) for the locus. 


## Inputs

The script at some point uses the following files, stored on Github. If their names or locations change, the script will have to be adjusted.

1) data/modified_samplemetadata.csv
2) data/environment/full_sample_metadata_4_20_19_ER.csv
3) data/PopPlotting_COLORS.csv

It also needs to read in a genotype matrix, which is too large to be stored on GitHub. The default is to read in: 

data/large_data/genotypeMatrix.rds

And then subset it. Remember, you will have to add this file locally after pulling the repository. A more efficient method if you do not need to change the way the matrix was subset is to comment out the line where "subsetGenoData" is applied and replace it with a line to read in the already subsetted genotype matrix. EX:

```
# wild <- subsetGenoData(all_data, all_metadata)
>wild <- readRDS("data/genotypeMatrix_selecting_Wild.rds")
```

You can find this pre-subsetted genotype matrix [on google drive](https://drive.google.com/open?id=14RkcdKIHvP0uX9moHsrpik3DRHwmPKAt)

## Modifying the script

If you want to change the variables that are being analyzed, simply modify the vector at lines 201-202

```
## Calc statistics and generate plots

envi_variables <- c("Lat","Long", "Temp_C", "Mean_Annual_Temp_Celsius", "Max_temperature_Celsius", "Min_temperature_Celsius",
                    "Mean_Annual_Salinity_ppt", "dd_0", "dd_30")
```

If you just want the table and not the plots, change the value of 'plots' in the call to calcEnviLFMMandSpRho() from TRUE to FALSE

```
out_table <- calcEnviLFMMandSpRho(envi_var = var, pop_object = wild, metadata = all_metadata, plots = TRUE)
```

If you want to adjust how the plots look, take a look at this block of code in calcEnviLFMMandSpRho() from lines 133-156

```
 if (plots){
    ### LF Pplot
    print("Plotting LF 1 and 2")
    png(paste("figures/6envi_assoc/LFMM_ridge_0.0", envi_var, "LF_plot.png", sep = "_"))
    plot(lfmm.ridge$U[,1], lfmm.ridge$U[,2], col = metadata[which(metadata$wild_for_assoc == 1),]$color, pch = 19, 
         main = paste("LFMM Ridge", envi_var,"Association"), xlab = "LF1", ylab = "LF2")
    text(lfmm.ridge$U[,1], lfmm.ridge$U[,2] + 20, labels = pop_object$Pop.ID, cex = 0.6)
    dev.off()
    rm(lfmm.test.ridge)
    rm(lfmm.ridge)

    # LFMM p-value plot
    LFMM_ridge_0.0_log10p <- -log10(as.numeric(p.values.ridge))
    rm(p.values.ridge)
    png(paste("figures/6envi_assoc/LFMM_ridge_0.0", envi_var, "pvalues_plot.png", sep = "_"))
    plot(pop_object$Pos, LFMM_ridge_0.0_log10p, main = "LFMM Ridge P-values")
    dev.off()
    
    # Spearmann's vs LFMM
    png(paste("figures/6envi_assoc/Spearmanns_vs_LFMM", envi_var, "plot.png", sep = "_"))
    plot(absSpCor, LFMM_ridge_0.0_log10p, main = "Spearmann's Rho vs LFMM Ridge")
    abline(lm(LFMM_ridge_0.0_log10p ~ as.vector(absSpCor)), col = "red")
    dev.off()
  }
```



