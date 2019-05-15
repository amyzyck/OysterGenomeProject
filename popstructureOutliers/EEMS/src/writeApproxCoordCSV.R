library(measurements)

args <- commandArgs(trailingOnly = TRUE)

datap    <- "SNP.TRSdp5g95FnDNAmaf05_50KThinnedRandomSNPs13PCs"
dataDir  <- "/media/kevin/TOSHIBA_EXT/oyster_genome_project/50KThinnedRandomSNPs13PCs/"

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

keep$Lat  <- unlist(lapply(keep$Lat, convertCoords))
keep$Long <- lapply(keep$Long, convertCoords)

keep$Long <- unlist(lapply(keep$Long, function(x) -abs(as.numeric(x))))

write.csv(keep, file = "/media/kevin/TOSHIBA_EXT/oyster_genome_project/modified_metadata_w_approx_coords.csv", quote = FALSE, 
          row.names = FALSE)