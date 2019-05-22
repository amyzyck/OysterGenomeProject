#### Package Install ####
    PackagesExists <- require (c("png", "profvis"))
    library(grid)
    library(qvalue)
    
    if (!PackageExists) {
    install.packages (c("png", "profvis"))
    library (c("png", "profvis"))
    }
    #########################

add.breakpoints <- function (chr.index, ybottom, ytop, ytext) {
       
    breakpoints <- read.table("breakpoints_for_plots.txt")

    breakpoints <- breakpoints[which(breakpoints$Chr == chr.index), ]

    for (i in 1:nrow(breakpoints)) {
        rect(xleft=breakpoints$Low[i], 
             ybottom = ybottom, 
             xright = breakpoints$High[i],
             ytop = ytop,
             col="lightgrey",
             border="transparent")
        x_text <- floor((breakpoints$High[i]+breakpoints$Low[i])/2)
        text(x = x_text, y = ytext , paste("LG", breakpoints$LG[i], sep=" "))
    }
}

make.plots <- function (data = NULL, subset = NULL, chr.index = NULL) {
    
    #######################
    # Load data and Setup #
    #######################

    # data - subset of full data given plots are made chromosome at a time    
    data <- data[which(data$Chr == chr.index), ]
    
    # plots.path - path for plots to be saved given the subset specified
    plots.path <- file.path("figures", "4outlier", subset)
    
    if (!dir.exists(plots.path)) {
        dir.create(plots.path)
    }

    # misassembled - chromosome indices where misassemblies occurred and need 
    #                shown on plots
    misassembled <- c(5, 6, 9)

    # He - Heterozygosity calculated by OutFLANK given the subset
    He <- paste("OutFLANK_0.2_bigsnpr_0.8.2_He_", subset, sep = "")
    # FST - FST values calculated by OutFLANK given the subset
    FST <- paste("OutFLANK_0.2_bigsnpr_0.8.2_FST_", subset, sep = "")
    # outflank.outlier & pcadapt.outlier - outlier flags for the seperate methods.
    outflank.outlier <- paste("OutFLANK_0.2_bigsnpr_0.8.2_OutlierFlag_", 
                              subset, 
                              sep = "")
    pcadapt.outlier <- paste("pcadapt_4.0.3_bigsnpr_0.8.2_OutlierFlag_", 
                             subset, 
                             sep = "")
    # negativeLog10p - Negative log10 p-values from PCAdapt.
    negativeLog10p <- paste("pcadapt_4.0.3_bigsnpr_0.8.2_negativeLog10p_", 
                            subset, 
                            sep = "")
    
    # Some of the negativeLog10p values are very large. This makes the plots 
    # look very odd and squishes all other, smaller values to the bottom.
    data[, negativeLog10p][data[, negativeLog10p] >= 100] <- 100
    # plot_x - All loci positions that have a heterozygosity greater than 0.1
    plots.x <- data$Pos[data[, He] > 0.1]

    ################################
    # Plotting - Seperate Outliers #
    ################################

    # OutFLANK results
    png(paste(plots.path, 
              "/manhattan_plots_", 
              sprintf(paste("%0", 2, "d", sep=""), 
              chr.index), 
              ".png", 
              sep = ""), 
        height = 1440,
        width = 2048) 
    par(mfrow = c(2, 1), mar = c(8, 4, 2, 1))
    plot(plots.x, 
         data[, FST][data[, He] > 0.1],
         ylab = "FST",
         xlab = "",
         ylim = c(-0.2, 1.1),
         xlim = c(min(plots.x, na.rm = TRUE), max(plots.x, na.rm = TRUE)),
         main = paste("OutFLANK Manhattan Plot: Chromsome", chr.index), 
         col = rgb(0, 0, 0, 0.2),
         xaxt = "n")
    axis(side = 1,
         las = 2,
         at = seq(min(plots.x, na.rm = TRUE), 
                  max(plots.x, na.rm = TRUE), 
                  length.out = 50))
    mtext(text = "Position (BP)",
          side = 1,
          line = 6)
    # Add misassembly if chromosome has one
    if (chr.index %in% misassembled) {
        add.breakpoints(chr.index, ybottom = -0.25, ytop = 1.15, ytext = 1.1)
        points(data$Pos[data[ , He] > 0.1], data[, FST][data[ , He] > 0.1])
    }
    # Colored points for outliers
    points(data$Pos[which(data[, outflank.outlier] == TRUE)], 
           data[, FST][which(data[, outflank.outlier] == TRUE)], 
           col = "red", 
           pch = 20)

    # PCAdapt results
    plot(plots.x, 
         data[, negativeLog10p][data[, He] > 0.1],
         ylab = "Negative Log10 p-value",
         xlab = "",
         ylim = c(0, 100),
         xlim = c(min(plots.x, na.rm = TRUE), max(plots.x, na.rm = TRUE)),
         main = paste("PCAdapt Manhattan Plot: Chromsome", chr.index), 
         col = rgb(0, 0, 0, 0.2),
         xaxt = "n")
    axis(side = 1,
         las = 2,
         at = seq(min(plots.x, na.rm = TRUE), 
                  max(plots.x, na.rm = TRUE), 
                  length.out = 50))
    mtext(text = "Position (BP)",
          side = 1,
          line = 6)
    # Add misassembly if chromosome has one
    if (chr.index %in% misassembled) {
        add.breakpoints(chr.index, ybottom = 0, ytop = 103, ytext = 100)
        points(data$Pos[data[, He] > 0.1], 
               data[, negativeLog10p][data[ ,He] > 0.1])
    }
    # Add misassembly if chromosome has one
    points(data$Pos[which(data[, pcadapt.outlier] == TRUE)], 
           data[, negativeLog10p][which(data[, pcadapt.outlier] == TRUE)], 
           col = "blue", 
           pch = 20)

    dev.off()

    #################################
    # Plotting - Outlier Comparison #
    #################################

    png(paste(plots.path, 
              "/outlier_compare_chr_", 
              sprintf(paste("%0", 2, "d", sep=""), 
              chr.index), 
              ".png", 
              sep = ""), 
        height = 1440,
        width = 2048) 
    par(mfrow = c(2, 1), mar = c(8, 4, 2, 16), xpd = NA)
    # OutFLANK results
    plot(plots.x, 
         data[, FST][data[, He] > 0.1],
         ylab = "FST",
         xlab = "",
         ylim = c(-0.2, 1.1),
         xlim = c(min(plots.x, na.rm = TRUE), max(plots.x, na.rm = TRUE)),
         main = paste("OutFLANK Manhattan Plot: Chromsome", chr.index), 
         col = rgb(0, 0, 0, 0.2),
         xaxt = "n")
    axis(side = 1,
         las = 2,
         at = seq(min(plots.x, na.rm = TRUE), 
                  max(plots.x, na.rm = TRUE), 
                  length.out = 50))
    mtext(text = "Position (BP)",
          side = 1,
          line = 6)
    # Add misassembly if chromosome has one
    if (chr.index %in% misassembled) {
        add.breakpoints(chr.index, ybottom = -0.25, ytop = 1.15, ytext = 1.1)
        points(data$Pos[data[ , He] > 0.1], data[, FST][data[ , He] > 0.1])
    }
    # Colored points for outliers
    points(data$Pos[which(data[, outflank.outlier] == TRUE & data[, pcadapt.outlier] != TRUE)], 
           data[, FST][which(data[, outflank.outlier] == TRUE & data[, pcadapt.outlier] != TRUE)], 
           col = "red", 
           pch = 20)
    points(data$Pos[which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] != TRUE)], 
           data[, FST][which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] != TRUE)], 
           col = "blue", 
           pch = 20)
    points(data$Pos[which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] == TRUE)], 
           data[, FST][which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] == TRUE)], 
           col = "orange", 
           pch = 20)

    # PCAdapt results
    plot(plots.x, 
         data[, negativeLog10p][data[, He] > 0.1],
         ylab = "Negative Log10 p-value",
         xlab = "",
         ylim = c(0, 100),
         xlim = c(min(plots.x, na.rm = TRUE), max(plots.x, na.rm = TRUE)),
         main = paste("PCAdapt Manhattan Plot: Chromsome", chr.index), 
         col = rgb(0, 0, 0, 0.2),
         xaxt = "n")
    axis(side = 1,
         las = 2,
         at = seq(min(plots.x, na.rm = TRUE), 
                  max(plots.x, na.rm = TRUE), 
                  length.out = 50))
    mtext(text = "Position (BP)",
          side = 1,
          line = 6)
    # Add misassembly if chromosome has one
    if (chr.index %in% misassembled) {
        add.breakpoints(chr.index, ybottom = 0, ytop = 103, ytext = 100)
        points(data$Pos[data[, He] > 0.1], 
               data[, negativeLog10p][data[ ,He] > 0.1])
    }
    # Colored points for outliers
    points(data$Pos[which(data[, outflank.outlier] == TRUE & data[, pcadapt.outlier] != TRUE)], 
           data[, negativeLog10p][which(data[, outflank.outlier] == TRUE & data[, pcadapt.outlier] != TRUE)], 
           col = "red", 
           pch = 20)
    points(data$Pos[which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] != TRUE)], 
           data[, negativeLog10p][which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] != TRUE)], 
           col = "blue", 
           pch = 20)
    points(data$Pos[which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] == TRUE)], 
           data[, negativeLog10p][which(data[, pcadapt.outlier] == TRUE & data[, outflank.outlier] == TRUE)], 
           col = "orange", 
           pch = 20)

    # Add legend to top right, outside plot region
    legend("topright",
           inset = c(-0.12, -0.2),
           legend = c("OutFLANK","PCAdapt", "Both"),
           col = c("red", "blue", "orange"),
           pch = 20,
           title = "Outlier Detection Method",
           cex = 1.5,
           box.lty = 0)
    

    dev.off()

}

plots.wrapper <- function () {
    data <- read.table("AllOutlier_WildForAssocEnviAssoc_MergedData_Lotterhos.txt",
                       header = TRUE)

    subsets <- c("exclude_LM", 
                 "unrelated", 
                 "atlantic_selection_subset", 
                 "atlantic_wild_subset")
    
    for (i in subsets) {
        for (j in seq(1, 10)) {
            make.plots(data = data, subset = i, chr.index = j)
        }
    }
}

plots.wrapper()