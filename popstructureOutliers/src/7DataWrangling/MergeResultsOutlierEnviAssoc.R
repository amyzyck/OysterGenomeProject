##############################
# Combining all LFMM Results #
##############################

merge.results <- function () {

    analyses.path <- "data/large_outputs/envi_assoc/"
    if (!dir.exists(analyses.path)){
        stop("You have no outlier analyses to merge.")
    } else {
        analyses.files <- paste(analyses.path, 
                                list.files(analyses.path), 
                                sep = "")
    }

    data <- NULL
    for (i in 1:length(analyses.files)) {
        if (is.null(data)) {
            data <- read.table(analyses.files[i], 
                               header = TRUE,
                               stringsAsFactors = FALSE)
            } else {
                tmp <- read.table(analyses.files[i], 
                                  header = TRUE, 
                                  stringsAsFactors = FALSE)
                data <- full_join(data, tmp, by = c("Pos", "Chr", "Unique"))
                rm(tmp)
            }
        }

    colnames(data)[which(names(data) == "Unique")] <- "unique"

    write.table(data, "data/large_outputs/envi_assoc/LFMM_merged_data_Lotterhos.txt")
}

outlier <- read.table("data/large_outputs/outlier_analysis_merged_data_Lotterhos.txt", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)
env <- read.table("data/large_outputs/envi_assoc/LFMM_merged_data_Lotterhos.txt", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

combined <- full_join(outlier, env, by = c("Pos", "Chr", "unique"))
