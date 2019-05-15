### KE Lotterhos
### 20181024
### plot moving average

setwd("/Users/katie/Desktop/OysterGenomeProject/popstructureOutliers/LDanalysis")

ld <- read.csv("geno_ld_window_5000-5000.csv")

str(ld)
ld$R.2 <- as.numeric(as.character(ld$R.2))
chrs <- levels(ld$CHR)

hist(ld$R.2)


for (i in 1:length(chrs)){

  ld_sub <- ld %>% filter(CHR==chrs[i])
  ld_sub <- ld_sub[order(ld_sub$POS1),]
  cat(chrs[i], "Number loci in LD calculation:", nrow(ld_sub), "\n")
  window_size <- 10000
  ld_sub$movave <- mov_ave(ld_sub$POS1, ld_sub$R.2,  window_size)

  d <- ggplot(ld_sub, aes(POS1, R.2))
  d2 <- ggplot(ld_sub, aes(POS1, movave))
  
  pdf(paste0(chrs[i],"_LDplot_5000-5000_",nrow(ld_sub),".pdf"), width=20, height=5)
  layout(1)
  par(mar=c(3,4,3,1))
  #plot(ld_sub$POS1, ld_sub$R.2)
  print(d + geom_hex(binwidth = c(200000, 0.05)) + 
    labs(x=paste(chrs[i], nrow(ld_sub)), 
         y = "R^2 for 5000-5000") +
    ylim(c(0,1))
  )
  # first element of bin width is bp
  
  #plot(ld_sub$POS1, ld_sub$movave, pch=19, col=rgb(0,0,1,0.1),
  #     main = paste(i, chrs[i]))
  print(d2 + geom_hex(binwidth = c(200000, 0.05)) +
  labs(x=paste(chrs[i], nrow(ld_sub)), 
       y = paste("R^2 for 5000-5000 moving window",
                              window_size)) +
    ylim(c(0,1))
  )
  
  dev.off()
}


###


mov_ave <- function(x_pos, y_response, window_size_bp){
  
  getmeanwindow <- function(j){
    lower <- (x_pos[j]-window_size_bp)
    upper <- (x_pos[j]+window_size_bp)
    return(mean(y_response[which( x_pos < upper & x_pos > lower)]))
  }
  
  return(sapply(1:length(x_pos), getmeanwindow))
}


