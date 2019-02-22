### KE Lotterhos
### 20181024
### plot moving average

setwd("/Users/katie/Desktop/OysterGenomeProject/popstructureOutliers/LDanalysis")

library(data.table)
library(tidyverse)
library(devtools)
library(fields)
devtools::install_github('johannesbjork/LaCroixColoR', force=TRUE)
mycol <- lacroix_palette("Pamplemousse", n = 50, type = "continuous")

ld <- fread("~/Google Drive/Eastern Oyster Genome/Population Structure/geno_ld_window_4500-5000.txt", header=TRUE)

str(ld)
chrs <- levels(as.factor(ld$CHR))
chrs
names(ld)[5] <- "R.2"

for (i in 1:length(chrs)){
  ld_sub <- ld %>% filter(CHR==chrs[i])
  print(c(chrs[i], (max(ld_sub$POS1) - min(ld_sub$POS1))/10^6))
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
#######

#system.time(mov_ave(ld_sub$POS1[1:10000], ld_sub$R.2[1:10000],  window_size))


for (i in 1:length(chrs)){

  ld_sub <- ld %>% filter(CHR==chrs[i])
  ld_sub <- ld_sub[order(ld_sub$POS1),]
  cat(chrs[i], "Number pairs in LD calculation:", nrow(ld_sub), "\n")
  
  ### parameters for moving average calculation
  window_size <- 10000
  
  by = nrow(ld_sub)/100000 # will give 100000 points for moving ave
    # 100,000 points takes about 7 minutes
  which_ave <- round(seq(1, nrow(ld_sub), by),0)
  every_bp <- ld_sub$POS1[which_ave[1:(length(which_ave)-1)]]- 
    ld_sub$POS1[which_ave[2:length(which_ave)]]
  cat("bp gap for moving ave",abs(mean(every_bp)), "\n")
  
  ld_sub$movave <- NA
  
  ld_sub$movave[which_ave] <- mov_ave(ld_sub$POS1[which_ave], ld_sub$R.2[which_ave],  window_size)

  d <- ggplot(ld_sub, aes(POS1, R.2))
  d2 <- ggplot(ld_sub, aes(POS1, movave))
  
  pdf(paste0(chrs[i],"LDplot_4500-5000",".pdf"), width=20, height=5)
  layout(1)
  par(mar=c(3,4,3,1))
  #plot(ld_sub$POS1, ld_sub$R.2)
  print(d + geom_hex(binwidth = c(200000, 0.05)) + 
    labs(x=paste(chrs[i], nrow(ld_sub)), 
         y = "R^2 for 4500-5000") +
    ylim(c(0,1))  + scale_fill_gradientn(colors=two.colors(100, start="#56B1F7", end="red",
                                                  middle="black"))
      #scale_fill_gradient(low = "#56B1F7", high = "black")

  )  

  # first element of bin width is bp
  
  #plot(ld_sub$POS1, ld_sub$movave, pch=19, col=rgb(0,0,1,0.1),
  #     main = paste(i, chrs[i]))
  print(d2 + geom_hex(binwidth = c(200000, 0.05)) +
  labs(x=paste(chrs[i], nrow(ld_sub)), 
       y = paste("R^2 for 4500-5000 moving window",
                              window_size)) +
    ylim(c(0,1)) + scale_fill_gradientn(colors=two.colors(100, start="#56B1F7", end="red",
                                                          middle="black"))
  )
  
  hist(ld_sub$R.2, breaks=seq(0,1, by=0.01))  
  abline(v=mean(ld_sub$R.2, na.rm=TRUE), col="blue")
  
  dev.off()
}




