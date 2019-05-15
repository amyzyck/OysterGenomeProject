###################################################################################################
#
# File    : createCoordForEEMs.R
# History : 10/11/2018 - Created by Kevin Freeman (KF)
#
###################################################################################################
#
# This script is used to turn the metadata file (modifed_samplemetadata) into a .coord input file
# for EEMs. It adds approxLat and approxLong for samples that do not have coordinates. These 
# approx coordinates were manually obtained from google maps by searching the name of the 
# locations. These coordinates are converted to decimal degrees.
#
# The script also reads the .order file that indicates what order the samples are in in the vcf.
# It uses this and the metatdata to select only the wild samples and only write these to the .coord
# files. It then goes back to the diff matrix and rewrites it, removing all the non-wild samples.
#
###################################################################################################
library(measurements)

args <- commandArgs(trailingOnly = TRUE)

datap    <-'SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs' #args[1]
dataDir  <- '/media/kevin/TOSHIBA_EXT/oyster_genome_project/10KThinnedRandomSNPs13PCs' #"."

#if (length(args) != 1 || !file.exists(paste0(dataDir,"/",datap, ".diffs"))){
#  cat("Error: 1 argument expected.\n
#  Usage: createCoordForEEMS.R <data prefix>\n\nData prefix should be the name of the .diffs file without the .diffs extension
#  You may also want to check that you are running createCoordForEEMS.R in the directory that your data is in.")
#  quit()
#}


metadata <- read.csv(paste0(dataDir,"/modified_samplemetadata.csv"), header = TRUE, stringsAsFactors = FALSE)
order    <- read.csv(paste0(dataDir, "/", datap, ".order"), header = FALSE, sep = " ")
outfile  <- paste0(dataDir, "/", datap, ".coord")
diffs    <- paste0(dataDir, "/", datap, ".diffs")

## select only wild oyster samples
keep <- metadata[which (metadata$Wild.Sel == "W"), ]

## add indexes 
ind           <- as.list(1:length(order$V2))
names(ind)    <- order$V2
i             <- match("LM_1_pool", names(ind))
names(ind)[i] <- "LM_1"                    ### there's one inconsistency b/w the 2 files -- LM_1 is LM_1_pool in the order file
keep$vcfInd   <- unlist(ind[keep$Sample.ID])

## input latitudes and longitudes (approx.) for unlisted locations
lats <- list("39 14.01", "37 37.01", "29 14.57", "29 54.52")
longs <- list("-75 01.52", "-75 39.56", "-90 55.16", "-93 17.07")

locationNames <- c("Cape Shore", "Hummock Cove", "Sister Lake", "Calcasieu Lake")

names(lats)  <- locationNames
names(longs) <- locationNames 

appLats                             <- lats[keep$Pop]
appLats[sapply(appLats, is.null)]   <- NA
appLongs                            <- longs[keep$Pop]
appLongs[sapply(appLongs, is.null)] <- NA

## convert to decimal degrees
appLats   <- lapply(appLats, function(x) measurements::conv_unit(x, from = 'deg_dec_min', to = 'dec_deg'))
appLongs  <- lapply(appLongs, function(x) measurements::conv_unit(x, from = 'deg_dec_min', to = 'dec_deg'))

keep$ApproxLat  <- unlist(appLats, use.names = FALSE)
keep$ApproxLong <- unlist(appLongs, use.names = FALSE)

convertCoords <- function(x){
  if (grepl("'", x, fixed = TRUE)){
    x <- gsub("��", " ", x)
    x <- gsub("\\.", "", x)
    x <- gsub("' *", ".", x)
    x <- gsub("\" *[NW]*","",x)
    x <- gsub("  ", " ", x)
    x <- trimws(x)
    measurements::conv_unit(x, from = 'deg_dec_min', to = 'dec_deg')
  }
  else {
    x
  }
}

keep$Lat  <- lapply(keep$Lat, convertCoords)
keep$Long <- lapply(keep$Long, convertCoords)

keep$Long <- lapply(keep$Long, function(x) -abs(as.numeric(x)))

### print the lats and longs to a coord file
cmd <- paste("rm", outfile)
if(file.exists(outfile)){
  system(cmd, wait = TRUE)
}

for (i in (sort(keep$vcfInd))){
  sample <- keep[keep$vcfInd == i, ]
  print(i)
  print(sample$Sample.ID)
  if (is.na(sample$Lat)) {
    write(paste(sample$ApproxLong, sample$ApproxLat, sample$Sample.ID), file = outfile, append = TRUE)
  }
  else {
    write(paste(sample$Long, sample$Lat, sample$Sample.ID), file = outfile, append = TRUE)
  }
}

## edit the .diffs file -- it currently has info for every sample but we just want the wild samples
diffMatrix  <- read.table(diffs, header=FALSE, sep = " ")
if (dim(diffMatrix)[1] != dim(keep)[1]){
  diffMatrix  <- diffMatrix[, colSums(is.na(diffMatrix)) == 0]  # remove the 'NA' column that is created because of leading whitespace in file
  keepMatrix  <- diffMatrix[keep$vcfInd, keep$vcfInd]

  write.table(keepMatrix, file = diffs, row.names = FALSE, col.names = FALSE)
}