### Libraries
library(qvalue)
library(lfmm)
library(psych)
library(data.table)
library(dplyr)
library(factoextra)
library(ape)
library(tidyverse)
library(R.utils)
library(ncdf4)

### Setting Directory and Specifying Files  
# Set working directory 
print("Reading in data...")
setwd("/home/downeyam/Github/OysterGenomeProject/popstructureOutliers/")


# CpGOA File (gene, exon, or CDS)
cpOE_dataTable <- fread("http://gannet.fish.washington.edu/Atumefaciens/20190225_cpg_oe/ID_CpG_labelled_all",header = TRUE)
# Gff file for gene_id
#dat <- read.gff("data/GCF_002022765.2_C_virginica-3.0_genomic.genes.gff3",GFF3 = TRUE)
dat <- read.gff("data/GCF_002022765.2_C_virginica-3.0_genomic.CDS.gff3",GFF3 = TRUE)
dat$cpgdat.ID <- paste0(dat$seqid,":",dat$start-1,"-",dat$end)
head(dat)
## NOTE: The cpOE dataTable and gff need to reference the same gene features (i.e. genes or exons)

# Provides a brief sanity check to confirm that the cpgOE matrix is only matching to a single feature within the .gff3
# data file, and will also through a warning if the features in the matrix poorly matches features in the .gff3, indicating
# that an incorrect .gff3 file is being used.  
tab <- tabulate(match(cpOE_dataTable$ID,dat$cpgdat.ID))
if(max(tab)==1){
  print("No duplicate matches found in .gff file, proceeding with script....")
  if(sum(tab)/length(tab)<.8){
    print("Poor matching rate (<80%), confirm you are using the correct .gff3 file for the particular cpgOE matrix")
  }
}else{print("Exiting Script....features in cpgOE matrix mapping to multiple .gff features. Confirm you are using matching CpGOE and .gff files 
            (i.e. they are for the same genetic feature)");
  exit()
}

# Metadata
metadata       <- read.csv("data/modified_samplemetadata.csv", stringsAsFactors = FALSE, header = TRUE)
envi_metadata  <- read.csv("data/environment/full_sample_metadata_4_20_19_ER.csv", stringsAsFactors = FALSE, header = TRUE)
plot_metadata  <- read.csv("data/PopPlotting_COLORS.csv", stringsAsFactors = FALSE, header = TRUE)

########################################################################
#
# File: LFMMCgGOE.R
# History: 2019/05/21
#
# This script is used to perform LFMM ridge analysis  for association 
# between genotype and environmental factors in in CpGOE data at the gene level in wild type populations. 
# It starts with  a CpGOE matrix ...
# populationStructureScript.R and subsets it to only include wild type
# populations
#
#########################################################################

#////////////////////////////////////////////////////////////////////////
# combineMetadata is a function to combine the metedata from 
# modified_samplemetadata.csv with environmental data from
# full_sample_metadata_4_20_19_ER.csv
#-----------------------------------------------------------------------
# inputs: 
#  - metadata      -- modfied_samplemetadata.csv
#  - envi_metadata -- full_sample_metadata_4_20_19_ER.csv
# outputs:
#  - a dataframe that combines the data from the 2 files
#
#////////////////////////////////////////////////////////////////////////

combineMetadata <- function(metadata, envi_metadata, plot_colors){
  # make the Wild.Sel column match between the two metadata files
  envi_metadata$Wild.Sel[which(envi_metadata$Wild.Sel == "inbred")] <- "I"
  
  # combine the two metadata csv files
  common_cols <- intersect(names(envi_metadata), names(metadata))
  comb_metadata <-  merge(metadata, 
                          envi_metadata[which(envi_metadata$Sample.ID %in% metadata$Sample.ID), 
                                        which(! names(envi_metadata) %in% common_cols[-1])],   # index 1 is Sample.ID, keep that to merge 
                          by = "Sample.ID", all = T)
  comb_metadata <- merge(comb_metadata, plot_colors, by = "custom.id")
  comb_metadata <- comb_metadata[order(comb_metadata$vcf_order), ]
  return(comb_metadata)
}

#//////////////////////////////////////////////////////////////////////////////////
#
# calcEnviLFMMandSpRho does the bulk of the analysis. It calculates LFMM 
# p-values and Spearmann's Rho for the association between genotype and the 
# given environmental variable. Optionally, it produces plots to visualize 
# the results. It always saves the results in a dataframe with descriptive
# column names. 
# ---------------------------------------------------------------------------------
# Inputs: 
#        - envi_var   -- the name of the environmental variable to consider. Must
#                        match the name of one of the columns in metadata
#        - pop_object -- a list object of the type created by subsetGenoData,
#                        contains genotypes and pop IDs
#        - metadata   -- the metadata file that contains environmental variables
#                        and pop IDs
#        - plots      -- a boolean: TRUE means create plots, FALSE means only create
#                        the dataframe        
# Outputs: 
#       - a dataframe with columns for spearmann's rho values,  LFMM values, Chr,
#         Pos, and "Unique" (an identifier created from Chr and Pos)
#       - optional : plots that have been saved in figures/6envi_assoc
#///////////////////////////////////////////////////////////////////////////////////
out_table <- calcEnviLFMMandSpRho(envi_var = var[1], pop_object = cpOE_wild, metadata = wild_metadata,gene_ID = dat_reduce, plots = TRUE)
print("saving dataframe")
if (!dir.exists("data/cpgoe_enviAssoc_results")){
  dir.create("data/cpgoe_enviAssoc_results")
}

calcEnviLFMMandSpRho <- function(envi_var, pop_object, metadata,gene_ID, plots){
  # scale genotype matrix (REMOVED THE SCALING when running with CpOE matrix)
  #scaled.genotype <- scaled(as.matrix(t(pop_object)))
  cpOE <- as.matrix(t(pop_object))
  # create a temperature matrix
  envi        <- metadata[envi_var]
  #envi_matrix <- matrix(data = envi, nrow = length(envi), ncol = 1)
  scaled.envi <- scale(envi)
  
  # build lfmm ridge model
  print("Building lfmm ridge model")
  lfmm.ridge <- lfmm::lfmm_ridge(Y = cpOE, X = scaled.envi, K = 4, lambda = 1e-4)
  # perform association testing
  print("Running lfmm test")
  lfmm.test.ridge <- lfmm::lfmm_test(Y = cpOE, X = scaled.envi, lfmm = lfmm.ridge, calibrate = "gif")
  # get p values
  p.values.ridge <- lfmm.test.ridge$calibrated.pvalue
  
  print(c("lfmm.test.ridge$gif", lfmm.test.ridge$gif))
  
  
  # spearman's correlation
  print("Calculating spearman's correlation")
  corTest  <- corr.test(data.frame(scaled.envi), data.frame(cpOE), method = "spearman", ci = FALSE, adjust = "none")
  absSpCor     <- abs(corTest$r)
  absSpCor[is.na(absSpCor)] <- 0
  spCor_log10p <- -log10(corTest$p)
  spCor_log10p[is.na(spCor_log10p)] <- 0
  
  if (plots){
    ### LF Pplot
    print("Plotting LF 1 and 2")
    png(paste("figures/7cpgoe_enviAssoc/LFMM_ridge_0.0", envi_var, "LF_plot.png", sep = "_"))
    plot(lfmm.ridge$U[,1], lfmm.ridge$U[,2], col = metadata$color, pch = 19,
        main = paste("LFMM Ridge", envi_var,"Association"), xlab = "LF1", ylab = "LF2")
    text(lfmm.ridge$U[,1], lfmm.ridge$U[,2] + 20, labels = metadata$Pop, cex = 0.6)
    dev.off()
    rm(lfmm.test.ridge)
    rm(lfmm.ridge)
    
    # LFMM p-value plot
    LFMM_ridge_0.0_log10p <- -log10(as.numeric(p.values.ridge))
    LFMM_ridge_0.0_log10p[is.na(LFMM_ridge_0.0_log10p)] <- 0
    rm(p.values.ridge)
    png(paste("figures/7cpgoe_enviAssoc/LFMM_ridge_0.0", envi_var, "pvalues_plot.png", sep = "_"))
    #plot(row.names(pop_object),LFMM_ridge_0.0_log10p, main = paste(envi_var, "LFMM Ridge P-values"),)
    dev.off()
    
    # Spearmann's vs LFMM
    png(paste("figures/7cpgoe_enviAssoc/Spearmanns_vs_LFMM", envi_var, "plot.png", sep = "_"))
    #plot(absSpCor, LFMM_ridge_0.0_log10p, main = "Spearmann's Rho vs LFMM Ridge")
    #abline(lm(LFMM_ridge_0.0_log10p ~ as.vector(absSpCor)), col = "red")
    dev.off()
  }
  
  ### Save the results
  print("Creating data frame")
  unique_ID <- row.names(pop_object)#sprintf("%02d_%09d", pop_object$Chr, pop_object$Pos)
  split_ID <- unlist(strsplit(unique_ID,":"))
  chr_ID <- split_ID[c(TRUE,FALSE)]
  pos_ID <- split_ID[c(FALSE,TRUE)]
  pos_split_ID <- unlist(strsplit(pos_ID,"-"))
  start_ID <- pos_split_ID[c(TRUE,FALSE)]
  end_ID <- pos_split_ID[c(FALSE,TRUE)]
  
  # make data into a matrix
  stat_matrix <- matrix(c(chr_ID,start_ID,end_ID,unique_ID, 
                          LFMM_ridge_0.0_log10p, absSpCor, spCor_log10p), ncol = 7, nrow = length(unique_ID))
  # matrix to dataframe
  out_table <- as.data.frame(stat_matrix)
  
  # rename columns
  colnames(out_table) <- c("Chromosome","Start_Pos","End_Pos","cpgdat.ID", 
                           paste("LFMM_ridge_0.0_log10p_cpgoe", envi_var, "wild_for_assoc",sep = "_"),
                           paste("Spearmanns_Rho_ABS_cpgoe", envi_var, "wild_for_assoc", sep = "_"),
                           paste("Spearmanns_Rho_log10p_cpgoe", envi_var, "wild_for_assoc", sep = "_"))
  
  # Convert back into numeric or character
  out_table$cpgdat.ID <- as.character(out_table$cpgdat.ID)
  
  # Join gene ids and LOC
  out_table_v2<- left_join(out_table,gene_ID, by="cpgdat.ID")
  out_table_v2 <- out_table_v2[,c(9,4,1,2,3,8,5,6,7)]
  return(out_table_v2)
}

#####################################################################################

#### Read in and process the data

# Read in environmental data and other variables of interest
print("Processing metadata")

#### recode envi_metadata pressure variables
# low = 1, medium = 2, high = 3
envi_metadata$Dermo_pressure <- factor(envi_metadata$Dermo_pressure, levels(factor(envi_metadata$Dermo_pressure))[c(2,3,1)])
envi_metadata$Dermo_pressure <- as.numeric(envi_metadata$Dermo_pressure)
# none = 0, low = 1, sporadic = 2, high = 3
envi_metadata$MSX_pressure <- factor(envi_metadata$MSX_pressure, levels(factor(envi_metadata$MSX_pressure))[c(3,2,4,1)])
envi_metadata$MSX_pressure <- as.numeric(envi_metadata$MSX_pressure) - 1

all_metadata <- combineMetadata(metadata = metadata, envi_metadata = envi_metadata, plot_colors = plot_metadata)
wild_metadata <- all_metadata[all_metadata$wild_for_assoc==1,]
wild_ID <- all_metadata$Sample.ID[all_metadata$wild_for_assoc==1]


  #### Processing in cpOE matrix
print("Processing CpGOA matrix....")
colnames(cpOE_dataTable)[colnames(cpOE_dataTable) == "HC_VA"] <- c(paste0("HC_VA_",c(1:6)))
cpOE <- data.frame(cpOE_dataTable,row.names = 1)

cpOE_wild <- subset(cpOE,select=c(wild_ID))

#### Read in Gene Identifiers
gene_id <- unlist(strsplit(dat$attributes,";"))
gene_id <- gene_id[c(TRUE,FALSE)]
final_gene_id <- substr(gene_id,9,30)
dat$gene_id <- final_gene_id
dat  %>% select(cpgdat.ID,strand,gene_id) -> dat_reduce

### Make a scree plot to identify the K value to use
print("plotting PCA scree plot")
png("figures/7cpgoe_enviAssoc/PCA_scree_plot.png")
PCA <- prcomp(t(cpOE_wild))
plot(PCA, type = "l", main = "PCA Scree Plot")
dev.off()

print("plotting PCA biplot")
png("figures/7cpgoe_enviAssoc/PCA_biplot_plot.png")
groups <- as.factor(wild_metadata$Region)
# This requires library(factoextra)
fviz_pca_ind(PCA,
             col.ind = groups, # color by groups
             palette = c("red","green","blue","orange"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
dev.off()

## Calc statistics and generate plots

envi_variables <- c("Max_temperature_Celsius", "Min_temperature_Celsius",
                    "Mean_Annual_Salinity_ppt", "dd_0", "dd_30", "Lat",
                    "Long", "Temp_C", "Mean_Annual_Temperature_Celsius",
                    "Dermo_pressure", "MSX_pressure")
  
for (i in 1:length(envi_variables)){
  var <- envi_variables[i]
  message("----------------------------------------------------")
  message(paste0("STARTING ANALYSIS FOR: ", var))
  message("----------------------------------------------------")
  out_table <- calcEnviLFMMandSpRho(envi_var = var, pop_object = cpOE_wild, metadata = wild_metadata,gene_ID = dat_reduce, plots = TRUE)
  print("saving dataframe")
  if (!dir.exists("data/cpgoe_enviAssoc_results")){
    dir.create("data/cpgoe_enviAssoc_results")
  } 
  if(i == 1){
    running_table <- out_table
  }
  if(i > 1){
    temp <- out_table[,c(2,7:9)]
    running_table <- left_join(running_table,temp, by="cpgdat.ID")
  }
  write.table(out_table, file = paste0("data/cpgoe_enviAssoc_results/", var, "_assoc_results.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
}  
write.table(running_table, file = paste0("data/cpgoe_enviAssoc_results/all_assoc_results.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
#Also writting this table as .RData to compress is it slightly
saveRDS(running_table,file = paste0("data/cpgoe_enviAssoc_results/all_assoc_results.rds"))
# Read this version of the the full association table in with the readRDS() function.
