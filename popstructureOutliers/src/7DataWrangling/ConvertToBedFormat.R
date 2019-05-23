# CHR - Specific names for chromosomes that correspond to index in data
CHR <- c("NC_035780.1", "NC_035781.1", "NC_035782.1", "NC_035783.1", 
         "NC_035784.1", "NC_035785.1", "NC_035786.1", "NC_035787.1", 
         "NC_035788.1", "NC_035789.1")

# data - Load data to be converted to additional formats
data <- read.table("data/large_outputs/AllOutlier_WildForAssocEnviAssoc_MergedData_Lotterhos.txt", 
                   stringsAsFactors = FALSE,
                   header = TRUE)

#########################
# Convert to BED format #
#########################

# bed.tmp - initial dataframe with the required columns for the .bed format.
# From IGV website - 
#
# Zero-based index: Start and end positions are identified using a zero-based 
#                   index. The end position is excluded. For example, setting 
#                   start-end to 1-2 describes exactly one base, the second 
#                   base in the sequence.

bed.tmp <- data.frame("chrom" = CHR[data$Chr], 
                      "chromStart" = (data$Pos - 1),
                      "chromEnd" = data$Pos,
                      stringsAsFactors = FALSE)

# bed - final data frame that will be saved in the gwas format

# Columns from bed.tmp will be binded with the columns in loaded data that 
# have yet to be included.

bed <- cbind(bed.tmp, 
             data[, !names(data) %in% c("Chr", 
                                        "Pos")])

# Save the file to bed format.

write.table(bed, "AllOutlier_WildForAssocEnviAssoc_MergedData_Lotterhos.bed", sep = " ")