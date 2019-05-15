###### Fst functions for OutFLANK

####### Fst for haploids
#' 
#' Calculates FST with correction for local sample sizes, for haploid biallelic data. Based on Weir (1996 - Genetic Data Analysis II)
#' 
#' @title FST calculation for biallelic haploid data
#'
#' @param AllCounts This is an array with a row for each population, and two values per row: Number of alleles in the sample of one type,  number of alleles of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item   p_ave: the average allele frequency
#'  \item   FST:  Fst (with sample size correction)
#'  \item   T1: The numerator of the Fst calculation
#'  \item   T2: The denominator of the Fst calculation
#'  }
#'@export
#'  
WC_FST_FiniteSample_Haploids_2AllelesB_MCW <- function(AllCounts){
  #Input a matrix of the counts of each allele (columns) in each population (rows)
  #returns vector instead of list of Fst values, according to Weir
  
  n_pops<-dim(AllCounts)[1]
  r<-n_pops
  counts1 <- AllCounts[,1]
  sample_sizes <- rowSums(AllCounts)
  n_ave <- mean(as.numeric(sample_sizes))
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)  
  p_freqs = counts1/sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)

  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  
  T1 <- s2 - 1/(n_ave-1)*(p_ave*(1-p_ave) -(s2*(r-1)/r))
  T2 <- (n_c - 1)*(p_ave*(1-p_ave))/(n_ave-1) + (1 + (r-1)*(n_ave-n_c)/(n_ave-1))*s2/r
  
  FST <- T1/T2 
  
  return(c(He,p_ave, FST, T1, T2))
  
}

#####The following function does not correct for finite local sample size, and
#####therefore creates a biased estimate of Fst that does not go negative.

#' 
#' Calculates FST without correction for local sample sizes, for haploid biallelic data. This is necessary for using OutFLANK, which depends on these uncorrected values for reliable function. (Otherwise, sampling corrections can sometimes cause negative estimates of FST.)
#' 
#' @title FSTNoCorr calculation for biallelic haploid data
#'
#' @param AllCounts This is an array with a row for each population, and two values per row: Number of alleles in the sample of one type,  number of alleles of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   HeNoCorr:  the expected heterozygosity of the locus
#'  \item   p_aveNoCorr: the average allele frequency
#'  \item   FSTNoCorr:  Fst (without sample size correction)
#'  \item 	T1NoCorr: The numerator of the FSTNoCOrr calculation 
#'  \item   T2NoCorr: The denominator of the FSTNoCOrr calculation
#'  }
#'@export
#'  
WC_FST_FiniteSample_Haploids_2AllelesB_NoSamplingCorrection<-function(AllCounts){
  #Input a matrix of the counts of each allele (columns) in each population (rows)
  
  n_pops<-dim(AllCounts)[1]
  r<-n_pops
  counts1 <- AllCounts[,1]
  sample_sizes <- rowSums(AllCounts)
  n_ave <- mean(as.numeric(sample_sizes))
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)  
  p_freqs = counts1/sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)

  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  
  T1NoCorr <- s2 
  T2NoCorr <- s2/r+(p_ave*(1 - p_ave))
  
  FSTNoCorr <- T1NoCorr/T2NoCorr 
  
  return(c(HeNoCorr=He,p_aveNoCorr=p_ave, FSTNoCorr=FSTNoCorr, T1NoCorr=T1NoCorr, T2NoCorr=T2NoCorr))
  
}

fstBarCalculatorNoCorr=function(DataList){
  #Calculates mean FstNoCorr from the dataframe, using sum(T1NoCorr) / sum(T2NoCorr) as the estimate of mean Fst.
  #Uses only data for which qvalues > qthreshold (i.e. $OutlierFlag==FALSE)
  #Does not internally screen for low MAF or low He values (but that can be added by only sending the
  #  high MAF rows to this function)
  sum(DataList$T1NoCorr[which(!DataList$OutlierFlag)])/sum(DataList$T2NoCorr[which(!DataList$OutlierFlag)])
}

fstBarCalculator=function(DataList){
  #Calculates mean Fst from the dataframe, using sum(T1) / sum(T2) as the estimate of mean Fst.
  #Uses only data for which qvalues > qthreshold (i.e. $OutlierFlag==FALSE)
  #Does not internally screen for low MAF or low He values (but that can be added by only sending the
  #  high MAF rows to this function)
  sum(DataList$T1[which(!DataList$OutlierFlag)])/sum(DataList$T2[which(!DataList$OutlierFlag)])
}

#####From FODR-- Diploid Fst from Weir and Cockerham 1984##########
##########################################
## WC FST for infinite sample of diploid allele freqs; without a correction for local sample size
###########################################  


#' 
#' Calculates FST without correction for local sample sizes, for diploid biallelic data. This is necessary for using OutFLANK, which depends on these uncorrected values for reliable function. (Otherwise, sampling corrections can sometimes cause negative estimates of FST.)
#' 
#' @title FSTNoCorr calculation for biallelic diploid data
#'
#' @param Sample_Mat This is an array with a row for each population, and three values per row: Number of Homozygotes of one type, number of heterozygotes, number of homozygotes of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item 	FSTNoCorr:  Fst (without sample size correction)
#'  \item 	T1NoCorr: The numerator of the uncorrected sample size correction (similar to Weir and Cockerham 1984)
#'  \item   T2NoCorr: The denominator of the uncorrected sample size correction
#'  }
#'@export
#'  
WC_FST_FiniteSample_Diploids_2Alleles_NoCorr<-function(Sample_Mat){
  
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  if(s2==0){return(1); break}  
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  
  a <- n_ave/n_c*(s2)
  
  b <- (p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  return(list(He=He,FSTNoCorr=FST, T1NoCorr=a, T2NoCorr=(a+b+c)))
}

#############FSt for diploids with local sample size corrections###############

##########################################
## WC FST for infinite sample of diploid allele freqs
###########################################  
#' 
#' Calculates FST with correction for local sample sizes, for diploid biallelic data. 
#' 
#' @title FST calculation for biallelic diploid data
#'
#' @param Sample_Mat This is an array with a row for each population, and three values per row: Number of Homozygotes of one type, number of heterozygotes, number of homozygotes of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item   FST:  Fst (with sample size correction)
#'  \item 	T1: The numerator of the Fst calculation (a from Weir and Cockerham 1984)
#'  \item   T2NoCorr: The denominator of the Fst calculation (a+b+c from Weir and Cockerham 1984)
#'  }
#'@export
#' 
WC_FST_FiniteSample_Diploids_2Alleles<-function(Sample_Mat){
  
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)

  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  if(s2==0){return(1); break}	
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))
  
  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c)))
}

####Diploid Fst functions for OutFLANK

####### Fst functions for OutFLANK
#' 
#' Calculates FST both with and without correction for local sample sizes, for diploid biallelic data. Based on Weir and Cockerham (1984)
#' 
#' @title FST calculation for biallelic diploid data
#'
#' @param Sample_Mat This is an array with a row for each population. There should be three columns, with the numbers of individuals from that population which are homozygotes for one allele, heterozygotes, and homozygotes for the other allele. 
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  The expected heterozygosity of the locus
#'  \item   FST:  Fst (with sample size correction)
#'  \item   T1: The numerator of the Fst calculation
#'  \item   T2: The denominator of the Fst calculation
#'  \item   FSTNoCorr:  Fst (without sample size correction)
#'  \item   T1NoCorr: The numerator of the Fst calculation without sample size correction
#'  \item   T2NoCorr: The denominator of the Fst calculation without sample size correction
#'  \item   meanAlleleFreq: The mean allele frequency over all populations for this locus
#'  }
#'@export
#'  
WC_FST_Diploids_2Alleles<-function(Sample_Mat){
  ##Calculate both Fst and Fst NoCorr at the same time, from WC84
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  if(s2==0){return(0); break}  
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))
  
  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  aNoCorr <- n_ave/n_c*(s2)
  
  bNoCorr <- (p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave)/(4*n_ave)*h_ave)
  
  cNoCorr <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  FSTNoCorr = aNoCorr/(aNoCorr+bNoCorr+cNoCorr)
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c),FSTNoCorr=FSTNoCorr, T1NoCorr=aNoCorr, T2NoCorr=(aNoCorr+bNoCorr+cNoCorr),meanAlleleFreq = p_ave))
}

####### Create input file for OutFLANK
#' 
#' Creates OutFLANK input file from individual genotype info. 
#' 
#' @title Create OutFLANK input file
#'
#' @param SNPmat This is an array of genotypes with a row for each individual. There should be a column for each SNP, with the number of copies of the focal allele (0, 1, or 2) for that individual. If that individual is missing data for that SNP, there should be a 9, instead. 
#' @param locusNames A list of names for each SNP locus. There should be the same number of locus names as there are columns in SNPmat.
#' @param popNames A list of population names to give location for each individual. Typically multiple individuals will have the same popName. The list popNames should have the same length as the number of rows in SNPmat.
#' 
#' @return Returns a data frame in the form needed for the main OutFLANK function.
#' 
#'@export
#'  
MakeDiploidFSTMat.v2 = function(SNPmat,locusNames,popNames){
  # SNPmat is a matrix with individuals in rows and snps in columns
  # 0, 1, or 2 represent the number of copies of the focal allele, and 9 is for missing data
  # locusNames is a character vector of names of each SNP
  # popNames is a character vector with the population identifier for each individual 
  
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  
  ### Check that SNPmat has appropriate values (0, 1, 2, or 9, only)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  ls <- paste(snplevs, collapse="")
  if(ls!="012" & ls!="0129"){print("Error: Your snp matrix does not have 0,1, and 2"); break}
  
  ### Checking that locusNames and popNames have the same lengths as the columns and rows of SNPmat
  if(dim(SNPmat)[1]!=length(popname) ){
    print("Error: your population names do not match your SNP matrix")
    break}
  
  if(dim(SNPmat)[2]!=length(locusname)){
    print("Error:  your locus names do not match your SNP matrix")
    break}
  
  writeLines("Calculating FSTs, may take a few minutes...")
  
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow=nloci, ncol=8)
  for (i in 1:nloci){
    FSTmat[i,]=unlist(getFSTs_diploids(popname,SNPmat[,i]))
    if (i%%10000==0){print(paste(i, "done of", nloci))}
  }
  outTemp=as.data.frame(FSTmat)
  outTemp = cbind(locusname,outTemp)
  
  colnames(outTemp)= c("LocusName","He", "FST", "T1", "T2", "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return (outTemp)
  
}


### Calculates FST etc. from a single locus from a column of individual data
getFSTs_diploids = function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

####### Likelihood functions needed by OutFLANK


EffectiveNumberSamplesMLE=function(FstVect, Fstbar, NumberOfSamples, SmallestFstInTrimmedList, LargestFstInTrimmedList){
  #This function should find the maximum likelihood value 
  #of the effective number of samples, for a given list of
  #Fst values.
  
  #The FstVect should already have been purged of NaN values and of loci with 
  #too low heterozygosity or MAF. 
  
  sortedFst=FstVect[order(FstVect)]  
  
  #The Minimum Fst considered in the trimmed data is the larger of the amount
  #specified by the user or the mean FSt over 100. This is to prevent extremely
  #small Fsts from causing estimation errors (Especially when R interprets a
  #small Fst as FSt=0.)
  LowTrimPoint=max(Fstbar/100,SmallestFstInTrimmedList)
  
  trimmedFstVect =FstVect[which((FstVect>=LowTrimPoint)&(FstVect<=LargestFstInTrimmedList))]
  
  trimmedFstArray=as.array(trimmedFstVect)
  
  localNLLAllData=function(dfInferred){
    localNLLOneLocus=function(Fst){
      negLLdfFstTrim(Fst,dfInferred,Fstbar,LowTrimPoint,LargestFstInTrimmedList)
    }
    sum(localNLLOneLocus(trimmedFstVect))
  }
  
  optim(NumberOfSamples, localNLLAllData, lower=2, method="L-BFGS-B")$par
}

IncompleteGammaFunction=function(a, z) {
  #equivalence to Mathematica Gamma[a,z] according to 
  #   http://r.789695.n4.nabble.com/Incomplete-Gamma-function-td833545.html
  pgamma(z,a,lower=FALSE)*gamma(a)
}

negLLdfFstTrim=function(Fst, dfInferred, Fstbar, LowTrimPoint, HighTrimPoint){
  #Fst is the Fst from a locus, and dfInferred is the candidate value for the
  #degrees of freedom for the chi-squared distribution of neutral Fst, and
  #Fstbar is the mean Fst of all neutral loci (sequentially inferred from
  #non-outlier loci) LowTrimPoint and HighTrimPoint are the values of the lowest
  #and highest Fst values allowed to be included in the Fst list.
  #
  #Finds contribution to the negative log likelihood of a given locus' Fst for a
  #given dfInferred #CHECKED AGAINST MATHEMATICA DERIVATION##
  
  df=dfInferred
  
  1/(2*Fstbar)*(df * Fst +df * Fstbar * log(2) - df * Fstbar *log(df)-(df-2)*Fstbar * log(Fst)+df * Fstbar * log(Fstbar) + 2*Fstbar * log(-IncompleteGammaFunction(df/2,df*HighTrimPoint/(2*Fstbar))+IncompleteGammaFunction(df/2,df*LowTrimPoint/(2*Fstbar))))
}

#OutFLANK:  An Fst outlier approach by Mike Whitlock and Katie Lotterhos, University of British Columbia.
#Development supported by AdapTree, Genome Canada, Genome BC, and an NSERC Discovery Grant to MCW.

########################How to use OutFLANK##################

#This method looks for Fst outliers from a list of Fst's for different loci. It
#assumes that each locus has been genotyped in all populations with approximately equal coverage. 

#OutFLANK estimates the distribution of Fst based on a trimmed sample of Fst's. It
#assumes that the majority of loci in the center of the distribution are
#neutral and infers the shape of the distribution of neutral Fst using a trimmed set of
#loci. Loci with the highest and lowest Fst's are trimmed from the data set
#before this inference, and the distribution of Fst df/(mean Fst)  is assumed to
#follow a chi-square distribution. Based on this inferred distribution, each
#locus is given a q-value based on its quantile in the inferred null
#distribution.

#The main procedure is called OutFLANK -- see comments in that
#function immediately below for input and output formats. The other functions
#here are necessary and must be uploaded, but are not necessarily needed by the
#user directly.

#Steps:
# 1. Make sure you have the biocLite package on your computer.  Code for getting it is commented in the next section.
#
# 2. Load all the functions in this script.
#
# 3. Create a file that has a row for each locus in your data set, with the following columns:
#     $LocusName: a character string that uniquely names each locus. 
#     $FST: Fst calculated for this locus. (Kept here to report the unbiased Fst of the results) 
#     $T1: The numerator of the estimator for Fst (necessary, with $T2, to calculate mean Fst) 
#     $T2: The denominator of the estimator of Fst 
#     $FSTNoCorr: Fst calculated for this locus without sample size correction. (Used to find outliers) 
#     $T1NoCorr: The numerator of the estimator for Fst without sample size correction (necessary, 
#                with $T2NoCorr, to calculate mean Fst) 
#     $T2NoCorr: The denominator of the estimator of Fst without sample size correction 
#     $He: The heterozygosity of the locus (used to screen out low heterozygosity loci that have 
#                a different distribution) 

#FstNoCorr, T1NoCorr, and T2NoCorr can be calculated from function given below: 
#WC_FST_FiniteSample_Haploids_2AllelesB_MCW  for the haploid case or
#WC_FST_FiniteSample_Diploids_2Alleles_NoCorr for diploids.



#The procedure will return a list with the following elements:
#     $FSTbar:  the mean Fst of the data from loci with high enough heterozygosity
#     $dfInferred:     the effective number of populations in the data (equals df + 1)
#     $numberLowFstOutliers:  the number of loci flagged by the OutFLANK procedure as having significantly low Fst (i.e. with a q-value less than 0.05)
#     $numberHighFstOutliers:  the number of loci flagged by the OutFLANK procedure as having significantly high Fst (i.e. with a q-value less than 0.05)
#     $results:   a data frame with information about each locus


# This results dataframe includes all of the input data, plus the following columns:
#     $indexOrder: integer index giving the original order of rows in the input file
#     $GoodH: TRUE if the heterozygosity is above the threshold set; FALSE otherwise
#     $qvalues: q-value for locus against null hypothesis of neutrality
#     $pvalues: p-value for locus against null hypothesis of neutrality
#     $pvaluesRightTail: p-value for locus against null hypothesis of neutrality, based only on the right tail
#     $OutlierFlag: TRUE if locus is an outlier; FALSE otherwise 




#############LOAD NECESSARY PACKAGES############# 

#Download the biocLite package at first use. On subsequent uses, run library(qvalue) before 
#using functions in the rest of this file.

#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
#source("FST functions.R")

#############FUNCTIONS###################################
 #'
 #'Takes Fst data for a list of loci to find outliers, using a trimmed likelihood approach.
 #'
 #'This function should take in a dataframe ("FstDataFrame") that 
#'has columns for $LocusName,$Fst,$T1,$T2,$FstNoCorr, $T1NoCorr, $T2NoCorr,$H. It should return a dataframe 
#'with those same columns but also new columns for $LowOutlierFlag, $HighOutlierFlag, and $q.
#'
#'This function requires Fst's calculated without sample size correction. These
#'can be calculated, for example, with WC_FST_FiniteSample_Haploids_2AllelesB_NoSamplingCorrection in this package.
#'
#'This use of the biased FSTs is necessary for the trimming outlier approach 
#'with small samples, because the debiasing sometimes creates negative Fsts 
#'which do not fit into the chi-square distribution.

#'This will use FST's calculated without sample size correction for outlier tests.
#'Such FSTs will be biased upwards, but as long as the sample size is similar for
#'all loci, the resulting measures ought to be give similar results.

#'This use of the biased FSTs is necessary for the trimming outlier approach with
#'small samples, because the debiasing sometimes creates negative Fsts which do
#'not fit into the chi-square distribution.
#'
#'@title Fst outliers with trimming
#'
#'@param FstDataFrame A data frame that includes a row for each locus, with columns as follows: 
#'\itemize{
#'                   \item $LocusName: a character string that uniquely names each locus. 
#'                    \item $FST: Fst calculated for this locus. (Kept here to report the unbiased Fst of the results) 
#'                    \item $T1: The numerator of the estimator for Fst (necessary, with $T2, to calculate mean Fst) 
#'                    \item $T2: The denominator of the estimator of Fst 
#'                    \item $FSTNoCorr: Fst calculated for this locus without sample
#'                    size correction. (Used to find outliers) 
#'                    \item $T1NoCorr: The numerator of the estimator for Fst without sample size correction (necessary, with $T2, to 
#'                    calculate mean Fst) 
#'                    \item $T2NoCorr: The denominator of the estimator of Fst 
#'                    without sample size correction 
#'                    \item $He: The heterozygosity of the locus (used to screen out low heterozygosity loci that have a different distribution) 
#'                    }
#'                    
#'@param LeftTrimFraction The proportion of loci that are trimmed from the lower end of the range of Fst before the likelihood function is applied.
#' 
#'@param RightTrimFraction The proportion of loci that are trimmed from the upper end of the range of Fst before the likelihood funciton is applied.
#' 
#'@param Hmin The minimum heterozygosity required before including calculations from a locus.
#' 
#'@param NumberOfSamples The number of spatial locations included in the data set.
#' 
#'@param qthreshold The desired false discovery rate threshold for calculating q-values.
#' 
#'@return
#' 
#' The function returns a list with seven elements:
#' \itemize{
#'  \item   FSTbar: the mean FST inferred from loci not marked as outliers 
#'  \item 	FSTNoCorrbar: the mean FST (not corrected for sample size---gives an upwardly biased estimate of FST)
#'  \item 	dfInferred: the inferred number of degrees of freedom for the chi-square distribution of neutral FST
#'   \item  numberLowFstOutliers: Number of loci flagged as having a significantly low FST (not reliable)
#'   \item  numberHighFstOutliers: Number of loci identified as having significantly high FST
#'   \item  results:  a data frame with a row for each locus. This data frame includes all the original columns in the 
#'                    data set, and six new ones: 
#'                    \itemize{
#'              \item $indexOrder (the original order of the input data set),
#'              \item $GoodH (Boolean variable which is TRUE if the expected heterozygosity is greater than the Hmin set by input),
#'              \item $OutlierFlag (TRUE if the method identifies the locus as an outlier, FALSE otherwise), and 
#'              \item $q (the q-value for the test of neutrality for the locus)
#'              \item $pvalues (the p-value for the test of neutrality for the locus)
#'              \item $pvaluesRightTail the one-sided (right tail) p-value for a locus
#'              }
#'  }
#'@export
#'  
OutFLANK.v2=function(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples, qthreshold=0.05){
  
  #
  #
  #Setting up necessary columns in dataframe
  Fstdata= outputDFStarterNoCorr(FstDataFrame,Hmin)

  
  # making working dataframe with real Fst (no NAs), storing NAs to add back later
  # This also removes loci with He values lower than Hmin from the working data frame
  nonkeepers = which((is.na(Fstdata$FSTNoCorr))|(Fstdata$He<Hmin))
  if(length(nonkeepers)>0) 
      workingDataFrame = Fstdata[-nonkeepers,]
  else
      workingDataFrame = Fstdata
 
  storedDataFrameNA = Fstdata[nonkeepers,]
  
  
  #Finding upper and lower bounds for trimming (eliminating NAs, but not negative FSTs)
  sortedDataFrame=workingDataFrame[order(workingDataFrame$FSTNoCorr),]
  
  NLociTotal=length(sortedDataFrame$FSTNoCorr)
  SmallestKeeper=ceiling(NLociTotal*LeftTrimFraction)
  LargestKeeper=floor(NLociTotal*(1-RightTrimFraction))
  LowTrimPoint=sortedDataFrame$FSTNoCorr[[SmallestKeeper]]
  HighTrimPoint=sortedDataFrame$FSTNoCorr[[LargestKeeper]]
    
  
  if(LowTrimPoint<0) {writeLines("ERROR: The smallest FST in the trimmed set must be > 0. Please use a larger LeftTrimFraction."); return()}
  if(HighTrimPoint>=1) {writeLines("ERROR: The largest FST in the trimmed set must be < 1. Please use a larger RightTrimFraction."); return()}
  
  #finding dfInferred and Fstbar iteratively  
  putativeNeutralListTemp=ifelse(workingDataFrame$FSTNoCorr>0,TRUE,FALSE)
  
  oldOutlierFlag=rep(FALSE,NLociTotal)
  
  
  #Note: All negative FST loci are marked as putative outliers, which will need
  #to be tested with the coalescent model later. In the meantime, they are
  #removed so as to not confuse the likelihood function.
  
  keepGoing=TRUE
  count = 0
  #writeLines(paste(mean(workingDataFrame$FSTNoCorr[putativeNeutralListTemp])))
  
  while(keepGoing){
    count=count+1
    if(count>19) {
      keepGoing=FALSE  
      writeLines("Exceeded iteration maximum.") ###Try with increased maximum value for count two lines above.
    }
    
    FstbarNoCorrTemp=fstBarCalculatorNoCorr(workingDataFrame[putativeNeutralListTemp,])  

    dfInferredTemp=EffectiveNumberSamplesMLE(workingDataFrame$FSTNoCorr[putativeNeutralListTemp],FstbarNoCorrTemp,NumberOfSamples,LowTrimPoint,HighTrimPoint)
    workingDataFrame=pOutlierFinderChiSqNoCorr.v2(workingDataFrame,FstbarNoCorrTemp,dfInferredTemp,qthreshold, Hmin)

    #### mark all negative FSTs as outliers if lowest nonneg FST is outlier
    #### (because negative Fst estimates can't be evaluated through the
    #### chi-square approach on their own)
    if(any(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<LowTrimPoint])) workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<0]=TRUE
    
    ####Any loci previously marked as $OutlierFlag=TRUE remain so, even if the new iteration doesn't flag them as outliers
    #     workingDataFrame$OutlierFlag=!as.logical((!workingDataFrame$OutlierFlag)*(!oldOutlierFlag))
    
    #Resetting neutral list, and checking whether the outlier list has stabilized
    putativeNeutralListTemp=ifelse((!workingDataFrame$OutlierFlag),TRUE,FALSE)
    if(sum(putativeNeutralListTemp)==0) {writeLines("No loci in neutral list..."); return("FAIL")}
    
    if(identical(oldOutlierFlag,workingDataFrame$OutlierFlag)) keepGoing=FALSE
    
    ######if all in trimmed get IDed as outlier - return to user with warning
    if(all(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<LowTrimPoint])){
      writeLines("All loci with Fst below the lower (lefthand) trim point were marked as outliers. Re-run with larger LeftTrimFraction or smaller qthreshold.")
      return(0)
    }
    
    if(all(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr>HighTrimPoint])){
      writeLines("All loci with Fst above the upper (righthand) trim point were marked as outliers. Re-run with smaller RightTrimFraction or smaller qthreshold.")
      return(0)
    }
    
    oldOutlierFlag=workingDataFrame$OutlierFlag
    
    #writeLines(paste(as.character(count),"   ",as.character(sum(putativeNeutralListTemp))))
  }
  
  if(count>19) writeLines("Loop iteration limit exceeded.")
  
  numberLowFstOutliers=sum(workingDataFrame$OutlierFlag[(workingDataFrame$FSTNoCorr<LowTrimPoint)])
  numberHighFstOutliers=sum(workingDataFrame$OutlierFlag[(workingDataFrame$FSTNoCorr>HighTrimPoint)])
  
  FSTbar=fstBarCalculator(workingDataFrame[putativeNeutralListTemp,])  
  
  
  #merge NA list back to working list, and sort by original order	
  resultsDataFrame=rbind(workingDataFrame,storedDataFrameNA)
  resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrder),]
  #return new dataframe
  list(FSTbar=FSTbar,FSTNoCorrbar=FstbarNoCorrTemp,dfInferred=dfInferredTemp,numberLowFstOutliers=numberLowFstOutliers,numberHighFstOutliers=numberHighFstOutliers,results=resultsDataFrame)
}


outputDFStarterNoCorr=function(FstDataFrame,Hmin=0.1) {
  #This will take a given dataframe with $LocusName, $FST,$He, $T1,  $T2, etc. and 
  #    initialize $indexOrder,$GoodH,$OutlierFlag (to 0), and $q (to 1).
  
  #Hmin is the smallest allowable He for which a locus should be included in 
  #the initial calculations. By default this requires that a locus have 
  #heterozygosity equal to 10% or more.
  
  len=length(FstDataFrame$FSTNoCorr)
  indexOrder=seq(1,len)
  GoodH=ifelse(FstDataFrame$He<Hmin,"lowH","goodH")
  OutlierFlag=ifelse(is.na(FstDataFrame$FSTNoCorr),NA,FALSE)
  qvalues=rep(NA,len)
  pvalues=rep(NA,len)
  pvaluesRightTail=rep(NA,len)
  cbind(FstDataFrame, indexOrder, GoodH, qvalues,pvalues,pvaluesRightTail,OutlierFlag )
  
}


#' 
#' Calculates q-values for test of neutrality for a list of loci, using input of an inferred degrees of freedom for the chi-square and mean Neutral FST
#' 
#'@title q values for test of neutrality
#'
#'@param DataList A data frame with a row for each locus, that includes at least a column for $FSTNoCorr. It also helps if there is a column with an identifier for the locus. This dataframe should have empty columns called $qvalues and $OutlierFlag as well.
#' 
#'@param Fstbar Mean Fst (without sample size correction) as inferred from neutral loci or OutFLank 
#'  
#'@param dfInferred The inferred degrees of freedom of the chi-square distribution describing neutral Fst values.
#'  
#'@param qthreshold The threshold False Discovery Rate for calling a locus an outlier ( default = 0.05)
#'@param Hmin The threshold heterozygosity (H) below which loci will be removed  
#'@return Returns a data frame with the original data, and two new columns appended:
#' \itemize{
#' \item $qvalues the q-value for a locus
#' \item $OutlierFlag TRUE if q is less than the qthreshold; FALSE otherwise
#' \item $pvalues the p-value for a locus
#' \item $pvaluesRightTail the one-sided (right tail) p-value for a locus
#' }
#' 
#'@export
#'

pOutlierFinderChiSqNoCorr.v2=function(DataList, Fstbar, dfInferred, qthreshold=0.05, Hmin=0.1){
  #Finds outliers based on chi-squared distribution
  #Takes given values of dfInferred and Fstbar, and returns a list of p-values and q-values for all loci based on chi-square.
  #Assumes that the DataList input has a column called $FSTNoCorr and that empty columns exist for $qvalues and $OutlierFlag 
  
  #
  #
  #Divide DataList into 3 lists:  DataListGood has $FST>0; DataListNeg has cases where $FST <=0; and
  #   DataListNA has cases where $FST is NA.
  #DataListNeg is necessary to keep separate here because these cases do not have meaningful results with the chi-square approach;
  #   however, they do carry information.
  
  keepers = which((DataList$FSTNoCorr > 0) & (DataList$He >= Hmin))
  DataListGood = DataList[keepers,]
  DataListOthers = DataList[-keepers,]
  numOthers = length(DataListOthers[,1])

  #Putting NAs in the results columns for all loci that don'tmeet Hmin or positive Fst criteria
  DataListOthers$pvalues = rep(NA,numOthers)
  DataListOthers$pvaluesRightTail = rep(NA,numOthers)
  DataListOthers$qvalues = rep(NA,numOthers)
  DataListOthers$OutlierFlag = rep(NA,numOthers)
  
  #Calculating p values and q-values for loci with high enough He and postive Fst
  pList = pTwoSidedFromChiSq(DataListGood$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
  pListRightTail = 1-pchisq(DataListGood$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
  # Edited for the outlierAnalysis. pi0 changed from NULL to 1.  
  qtemp=qvalue(pListRightTail, fdr.level=qthreshold, pi0 = 1)
  #Note:  Using the bootstrap method here seems OK, but if this causes problems remove the pi0.method="bootstrap" in the previous line to revert to the default.
  
  DataListGood$pvalues = pList
  DataListGood$pvaluesRightTail = pListRightTail
  DataListGood$qvalues = qtemp$qvalues
  DataListGood$OutlierFlag = qtemp$significant
  
  #Combining the good and bad loci back and sorting
  resultsDataFrame = rbind(DataListGood,DataListOthers) 
  #resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrder),]
}


#' 
#' Calculates q-values for test of neutrality for a list of loci, using input of an inferred degrees of freedom for the chi-square and mean Neutral FST, and returns the results in the same row order as the input
#' 
#'@title q values for test of neutrality
#'
#'@param DataList A data frame with a row for each locus, that includes at least a column for $FSTNoCorr. It also helps if there is a column with an identifier for the locus. 
#' 
#'@param Fstbar Mean Fst (without sample size correction) as inferred from neutral loci or OutFLank 
#'  
#'@param dfInferred The inferred degrees of freedom of the chi-square distribution describing neutral Fst values.
#'  
#'@param qthreshold The threshold False Discovery Rate for calling a locus an outlier ( default = 0.05)
#'@param Hmin The threshold heterozygosity (H) below which loci will be removed  
#'@return Returns a data frame with the original data, and two new columns appended:
#' \itemize{
#' \item $qvalues the q-value for a locus
#' \item $OutlierFlag TRUE if q is less than the qthreshold; FALSE otherwise
#' \item $pvalues the p-value for a locus
#' \item $pvaluesRightTail the one-sided (right tail) p-value for a locus
#' }
#' 
#'@export
#'
pOutlierFinderInOrder=function(DataList, Fstbar, dfInferred, qthreshold=0.05, Hmin=0.1){
  
  #Assign a temporary index to each row
  
  len = length(DataList$FSTNoCorr)
  indexOrderTEMP = seq(1,len)
  
  DataListTEMP = cbind(DataList, indexOrderTEMP)
  
  #Calculate p and q values using pOutlierFinderChiSqNoCorr
  
  resultsDataFrame = pOutlierFinderChiSqNoCorr.v2(DataListTEMP, Fstbar, dfInferred, qthreshold, Hmin)
  
  #Sort to index and delete temporary index
  
  resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrderTEMP),]
  
  within(resultsDataFrame, rm(indexOrderTEMP))
  
}



#' 
#' Calculates P-values for test of neutrality for a list of loci, using input of an inferred degrees of freedom for the chi-square and mean Neutral FST
#' 
#'@title P-values for test of neutrality
#'
#'@param DataList A data frame with a row for each locus, that includes at least a column for $FSTNoCorr and $He. 
#'@param Fstbar Mean Fst (without sample size correction) as inferred from neutral loci or OutFLank 
#'  
#'@param dfInferred The inferred degrees of freedom of the chi-square distribution describing neutral Fst values.
#'@param Hmin Minimum heterozygosity (H) to exclude low H alleles
#'    
#'@return Returns a data frame with the original data, and two new columns appended:
#' \itemize{
#' \item $pvalues the p-value for a locus, with extremely large values of FST near 0
#' \item $pvaluesRightTail the one-sided (right tail) p-value for a locus
#' }
#' 
#' #@export
#'
#pChiSqNoCorr=function(DataList, Fstbar, dfInferred, Hmin=0.1){
  #Finds outliers based on chi-squared distribution
  #Takes given values of dfInferred and Fstbar, and returns a list of p-values and q-values for all loci based on chi-square.
  #Assumes that the DataList input has a column called $FSTNoCorr and that empty columns exist for $qvalues and $OutlierFlag 

  #Divide DataList into 3 lists:  DataListGood has $FST>0; DataListNeg has cases where $FST <=0; and
  #   DataListNA has cases where $FST is NA.
  #DataListNeg is necessary to keep separate here because these cases do not have meaningful results with the chi-square approach;
  #   however, they do carry information.

#  pList=1-pchisq(DataList$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
#  pList[DataList$He < Hmin] = NA
  # add negative FST 
#  return(data.frame(DataList, Pval=pList))
#}  




pTwoSidedFromChiSq=function(x,df){
  #Takes a value x, finds the two-sided p-value for comparison to a chi-square distribution with df degrees of freedom.
  pOneSided=pchisq(x,df)
  ifelse(pOneSided>.5,(1-pOneSided)*2,pOneSided*2)
}

