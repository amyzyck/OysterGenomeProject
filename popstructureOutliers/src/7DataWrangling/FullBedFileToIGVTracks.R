library(colorspace)

bed.tmp <- read.table("data/large_outputs/UnrelatedOutlier_WildForAssocEnviAssoc_MergedData_Lotterhos.bed",
                      sep = " ",
                      header = TRUE,
                      stringsAsFactors = TRUE )

#############################################
# Separate results to individual .bed files #
#############################################

for (i in 1:length(names(bed.tmp[, 5:ncol(bed.tmp)]))) {
    colors <- rainbow_hcl(length(names(bed.tmp[, 5:ncol(bed.tmp)])))
    tmp <- paste("C_virginica-3.0", names(bed.tmp[, 5:ncol(bed.tmp)])[i], sep = "_")
    track.name <- paste('track name=', '"', tmp, '"', sep = "")
    description <- paste('description=', '"', tmp, '"', sep = "")
    color <- paste("color=", '"', colors[i], '"', sep = "")

    tmp.txt <- paste(tmp, ".bed", sep = "")

    cat(paste(track.name, description, color, 'itemRgb="On" #gffTags', sep = " "), "\n", file = tmp.txt)
    write.table(bed.tmp[, c(1:3, i+4)], 
                file = tmp.txt,
                sep = "\t", 
                quote = FALSE, 
                append = TRUE,
                row.names = FALSE, 
                col.names = FALSE)
}