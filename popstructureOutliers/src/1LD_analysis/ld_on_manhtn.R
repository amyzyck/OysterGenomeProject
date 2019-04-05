library(data.table)
library(tidyverse)
library(devtools)
library(fields)
library(dplyr)
library(ggplot2)
library(LaCroixColoR)
mycol <- lacroix_palette("Pamplemousse", n = 50, type = "continuous")

load.ld.data <- function (path = NULL, ld.window = NULL) {
  if (!dir.exists(path)) {
      stop("ERROR: The path that you have provided is not a directory")
  } else {
      ld.files <- paste(path, list.files(path), sep = "")
  }
  
  ld <- fread(ld.files[which(grepl(ld.window, ld.files, fixed=TRUE) == TRUE)], 
              sep="\t",
              header=TRUE,
              stringsAsFactors=FALSE)
  
  return(ld)
}

#######
### Calculate moving average
mov_ave <- function(x_pos, y_response, window_size_bp){
  df <- data.frame(x_pos, y_response)
  getmeanwindow <- function(j){
    lower <- (x_pos[j]-window_size_bp)
    upper <- (x_pos[j]+window_size_bp)
    y_sub <- df$y_response[which(x_pos < upper & x_pos > lower)]
      # this is 10x faster than using dplyr and >60x faster than 'for loop'
    return(mean(y_sub, na.rm=TRUE))
  }
  
  return(sapply(1:length(x_pos), getmeanwindow))
    # I also tried supplying a subset of indexes to the sapply function
    # and it was an order of magnitude slower!
}

manhattan.ld.plot <- function (path = NULL,
                               thinned.snps = NULL,
                               ld.window = NULL,
                               subset = NULL) {
  # Load LD calculations
  ld <- load.ld.data(path = path, ld.window = ld.window)
  
  # Load thinned SNPs to visualise where they are located on each chromosome
  thinned.snps <- readRDS(thinned.snps)

  # Get the unique values for chromosomes
  chrs <- levels(as.factor(ld$CHR))
  # Reset name R^2 to R.2
  names(ld)[5] <- "R.2"

  # for (i in 1:length(chrs)){
  #   ld_sub <- ld %>% filter(CHR==chrs[i])
  #   print(c(chrs[i], (max(ld_sub$POS1) - min(ld_sub$POS1))/10^6))
  # }

  #system.time(mov_ave(ld_sub$POS1[1:10000], ld_sub$R.2[1:10000],  window_size))

  # misassembled - misassembled chromosomes that need to have rectangles added.
  misassembled <- c(5, 6, 9)
  
  # breakpoints - data table with the breakpoint locations. Will be used to show
  # where these breakpoints are in the plots.
  breakpoints <- read.table("breakpoints_for_plots.txt")

  for (i in 1:length(chrs)){
    ld_sub <- ld %>% filter(CHR == chrs[i])
    ld_sub <- ld_sub[order(ld_sub$POS1),]
    cat(chrs[i], "Number pairs in LD calculation:", nrow(ld_sub), "\n")

    # 'snp.pos' - Positions of thinned SNPs for the chromosome that will be plotted.
    snp.pos <- thinned.snps$Pos[which(thinned.snps$Chr == i)]

    ### parameters for moving average calculation
    window_size <- 10000

    by = nrow(ld_sub)/100000 # will give 100000 points for moving ave
      # 100,000 points takes about 7 minutes
    which_ave <- round(seq(1, nrow(ld_sub), by), 0)
    every_bp <- ld_sub$POS1[which_ave[1:(length(which_ave)-1)]]-ld_sub$POS1[which_ave[2:length(which_ave)]]
    cat("bp gap for moving ave", abs(mean(every_bp)), "\n")

    ld_sub$movave <- NA

    ld_sub$movave[which_ave] <- mov_ave(ld_sub$POS1[which_ave], ld_sub$R.2[which_ave], window_size)

    d <- ggplot(ld_sub, aes(POS1, R.2))
    d2 <- ggplot(ld_sub, aes(POS1, movave))

    results.path <- file.path("figures/1LD_analysis/ld_on_manhtn", subset, ld.window)
    if (!dir.exists(results.path)) {
      dir.create(results.path, recursive = TRUE)
    }

    pdf(paste0(results.path, "/", chrs[i], "_LDplot_", ld.window, "_", subset, ".pdf"), width=20, height=5)
    layout(1)
    par(mar = c(3,4,3,1))
    #plot(ld_sub$POS1, ld_sub$R.2)
    d <- d + geom_hex(binwidth = c(200000, 0.05)) + 
             labs(x = paste(chrs[i], nrow(ld_sub)), y = paste("R^2 for", ld.window, sep=" ")) + 
             ylim(c(0, 1)) + 
             scale_fill_gradientn(colors = two.colors(100, 
                                                      start = "#56B1F7", 
                                                      end = "red", 
                                                      middle = "black"))
    if (i %in% misassembled) {
      brk.tmp <- breakpoints[which(breakpoints$Chr == i), ]
      #ld_sub$Low <- brk.tmp$Low
      #ld_sub$High <- brk.tmp$High
      for (j in 1:nrow(brk.tmp)) {
        d <- d + annotate("rect",
                          xmin = brk.tmp$Low[j], 
                          xmax = brk.tmp$High[j], 
                          ymin = -Inf, 
                          ymax = Inf,
                          fill = "green", 
                          alpha = 0.3)
      }
    }
    d <- d + geom_vline(data = snp.pos, color = "lightgrey")
    print(d)

    # first element of bin width is bp

    #plot(ld_sub$POS1, ld_sub$movave, pch=19, col=rgb(0,0,1,0.1),
    #     main = paste(i, chrs[i]))
    d2 <- d2 + geom_hex(binwidth = c(200000, 0.05)) +
               labs(x = paste(chrs[i], nrow(ld_sub)), 
                    y = paste("R^2 for", ld.window, "moving window", window_size)) +
               ylim(c(0,1)) + 
               scale_fill_gradientn(colors = two.colors(100, 
                                                        start = "#56B1F7", 
                                                        end = "red", 
                                                        middle = "black"))
    if (i %in% misassembled) {
      brk.tmp <- breakpoints[which(breakpoints$Chr == i), ]
      #ld_sub$Low <- brk.tmp$Low
      #ld_sub$High <- brk.tmp$High
      for (j in 1:nrow(brk.tmp)) {
        d2 <- d2 + annotate("rect",
                            xmin = brk.tmp$Low[j], 
                            xmax = brk.tmp$High[j], 
                            ymin = -Inf, 
                            ymax = Inf,
                            fill = "green", 
                            alpha = 0.3)
      }
    }
    d2 <- d2 + geom_vline(data = snp.pos, color = "lightgrey")
    print(d2)

    hist(ld_sub$R.2, breaks = seq(0,1, by = 0.01))  
    abline(v = mean(ld_sub$R.2, na.rm=TRUE), col = "blue")

    dev.off()
  }
}

manhattan.ld.plot(path = "data/large_data/ldAnalysisData/excluding_LM/",
                  thinned.snps = "data/thinned_snps/thinnedMatrixAndMetaData50000Window_exclude_LM.rds",
                  ld.window = "49500-50000",
                  subset = "excluding_LM")