library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)
library(measurements)
library(sp)

datapath <- "/media/kevin/TOSHIBA_EXT/oyster_genome_project/10KThinnedRandomSNPs13PCs/"
metadata <- read.csv(paste0(datapath,"source_files/modified_samplemetadata.csv"), header = TRUE, stringsAsFactors = FALSE)
coords <- read.csv(paste0(datapath, "SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs.coord"), header = FALSE, sep = " ")
coords <- matrix(c(coords$V1, coords$V2), ncol = 2)
coords <- unique(coords)
coords[6,2] <- coords[6,2] + 0.25 ## offset hog island so that it is not on top of sherman marsh
coords[,2]  <- coords[,2] + .5
mcmcpath <- c(paste0(datapath, "eems_output/SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs-include_land-chain1"),
              paste0(datapath, "eems_output/SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs-include_land-chain2"),
              paste0(datapath, "eems_output/SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs-include_land-chain3"))
plotpath <- paste0(datapath, "plots/")

coords_merc <- sp::spTransform(SpatialPoints(coords, CRS(projection_none)), CRS(projection_mercator))
coords_merc <- coords_merc@coords

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

colors <- c("red", "green", "blue", "purple", "orange", "lightgreen", "black", "pink", "brown")

order         <- read.csv(paste0(datapath, "source_files/SNP.TRSdp5g95FnDNAmaf05_10KThinnedRandomSNPs13PCs.order"), header = FALSE, sep = " ")
keep          <- metadata[which (metadata$Wild.Sel == "W"), ]
ind           <- as.list(1:length(order$V2))
names(ind)    <- order$V2
i             <- match("LM_1_pool", names(ind))
names(ind)[i] <- "LM_1"                    ### there's one inconsistency b/w the 2 files -- LM_1 is LM_1_pool in the order file
keep$vcfInd   <- unlist(ind[keep$Sample.ID])
ordKeep       <- keep[order(keep$vcfInd), ]

labels <- unique(ordKeep$Pop)

map_world <- getMap()
map_NA <- map_world[which(map_world@data$continent == "North America"), ]
map_NA <- spTransform(map_NA, CRSobj = CRS(projection_mercator))

eems.plots(mcmcpath = mcmcpath, 
           plotpath = paste0(plotpath, "labels_land_500_demes"), 
           longlat = TRUE,
           projection.in = projection_none,
           projection.out = projection_mercator,
           add.map = TRUE,
           add.outline = TRUE,
           add.demes   = TRUE,
           m.plot.xy = { plot(map_NA, col = NA, add = TRUE);
             text(coords_merc, col = "green", labels = labels, font = 1, cex = 0.5); },
           col.outline = "black")
