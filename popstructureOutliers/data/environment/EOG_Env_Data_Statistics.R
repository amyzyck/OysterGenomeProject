# Script to Organize Environmental Data

############################## TO DO LIST ###############################
#To do: in file named: "EOG_Reseq_Sample Info_with_env_data"

# 1) Put sample IDs for those sequenced and remove the ones that were not sequenced (I started doing it, but only did a couple of sites).  
#   Delete rows for samples not sequenced (only delete in that file named above).

# 2) Latitude and Longitude = change to an homogeneous format (choose one, usually better if just one number).

# 3) Add the following columns to the file (do not delete the ones there)

#Columns in file (some to add):
#Individual_ID (individuals IDs same names and order as vcf)
#latitude
#longitude
#MAT(mean annual temperature)
#max_temp
#min_temp
#MAS (mean annual salinity)
#max_salinity
#min_salinity
#notes
#dd_0 (degree days below 0C)
#dd_15 (degree days above 15C)
#dd_30 (degree days above 30C)
#Dermo disease pressure since (approximate year - Gulf of Mexico, Chesapeake and Delaware since the 60s and high, Maine since the 90s and low)
#MSX disease pressure since (approximate year - no disease for Gulf of Mexico, Chesapeake and Delaware since the 60s, constant, Maine since the 90s, sporadic)

# 4) Process the environmental data in the folder “environmental data” to get the data for the columns above.  I have monthly data for Louisiana and Maine sites, still working on the other sites.  I will put them in the folder as I get them.

########################################################################################

# Load required packages
library(readxl)
library(tidyverse)
library(DescTools)
library(Rmisc)
library(emdbook)
library(bbmle)
library(lme4)
library(lattice)
library(ggplot2)
library(stringr)
library(data.table)

# Load in data

list_of_sequenced_samples <- read.csv(file="list_of_sequenced_samples.csv", header=TRUE)

full_sample_info <- read.csv(file="EOG_Reseq_Sample Info_9.17.2018_SAMPLE_info.csv", header=TRUE)

#### ITEM 1 ####

# Remove sample columns from full_sample_info that were not sequenced using Natural Join
sequenced_sample_info <- merge(list_of_sequenced_samples, full_sample_info, by="Sample.Name", all=FALSE)


#### ITEM 2 ####

# Change Latitude and longitude to homogenous format
sequenced_sample_info$Lat
unique(sequenced_sample_info$Lat)

# going to convert latitude and longitude in Excel!
write.csv(file="sequenced_sample_info.csv", sequenced_sample_info)

# Saved updated latitude and longitude in "sequenced_sampl_ingo_Updated_LAT_LONG.xlsx"
  
### Maine Sites emperature readings are in celsius not farenheit, not going to conver the values
# Remove erroneous reader temperature values of 99C for both spreadsheets 
library(readxl)
Hog_island_combined <- read_xlsx("./Environmental_Data/HI_Hog_Island_Maine/Temp_Sal_Hog Island_calculations.xlsx", sheet="Hog Island Combined")
Sherman_Marsh_Sheepscot <- read_xlsx("./Environmental_Data/HI_Hog_Island_Maine/Temp_Sal_Hog Island_calculations.xlsx", sheet="Sheepscot_Combined")

Hog_island_combined_removed <- Hog_island_combined %>% filter(Temp != "99")
  Sherman_Marsh_Sheepscot_removed <- Sherman_Marsh_Sheepscot %>% filter(TEMP_C != "99.9")

# Average temperature in C, min temperature in C, max temperature
summary(Hog_island_combined_removed$Temp)
summary(Hog_island_combined_removed$Sal)
summary(Sherman_Marsh_Sheepscot_removed$TEMP_C)
summary(Sherman_Marsh_Sheepscot_removed$SALINITY_PCT)

# Sherman Marsh average Temp and salinity in 8/17 when oysters were collected
Sherman_Marsh_Sheepscot_removed[grepl("2017-",Sherman_Marsh_Sheepscot_removed$EFFORT_START_DATE),] 
  # Latest measurement is 7/25/2017, filling in the temp and salinity measurements from this date
  
### Average Max, min, temp, Average salinity for Calcasieu lake and Caillou Lake from provided spreadsheets rather than online

CL_Calcasieu_Lake <- read_excel("./Environmental_Data/CL_Calcasieu_Lake_Env_OBOYS2/Louisiana_Calcasieu Lake LOU monthly env data_edited.xlsx", sheet="Calcasieu_lake")
CL_Calcasieu_Lake_temp <- CL_Calcasieu_Lake %>% filter(parameter_cd == "temp_C") 
summary(CL_Calcasieu_Lake_temp$mean_va)
CL_Calcasieu_Lake_sal <- CL_Calcasieu_Lake %>% filter(parameter_cd =="salinity_ppt")
summary(CL_Calcasieu_Lake_sal$mean_va)

SL_Sister_Lake_Caillou_Lake <- read_xlsx("./Environmental_Data/SL Sister Lake_Caillou Lake/Louisiana_Caillou Lake monthly env data.xlsx", sheet="Caillou_Lake")
SL_Sister_Lake_Caillou_Lake_temp <- SL_Sister_Lake_Caillou_Lake %>% filter(parameter_cd == "temp_C") 
summary(SL_Sister_Lake_Caillou_Lake_temp$mean_va)

SL_Sister_Lake_Caillou_Lake_sal <- SL_Sister_Lake_Caillou_Lake %>% filter(parameter_cd =="salinity_ppt")
summary(SL_Sister_Lake_Caillou_Lake_sal$mean_va)

### Compile and calculate Cape Shore Average Max, min, temp and average salinity 

# 2013 2014 data with temp in F
Cape_Shore_2013_2014 <- read_excel("./Environmental_Data/CS Cape Shore DelBay High Sal_NEH_RU inbred/Cape_Shore_2013_May_2014_Jun.xlsx", skip = 2)
Cape_Shore_2013_2014 <- Cape_Shore_2013_2014[,c(2:3)] #select needed rows
Cape_Shore_2013_2014 <- Cape_Shore_2013_2014 %>% drop_na() # drop empty rows
colnames(Cape_Shore_2013_2014) <- c("Date_Time", "temp_F")
# convert farenheit to celsius
Cape_Shore_2013_2014$temp_C <- UnitConv(Cape_Shore_2013_2014$temp_F, "F","C")

# 2016 Apr_Oct data also has temp in F
Cape_Shore_2016_Apr_Oct <- read_excel("./Environmental_Data/CS Cape Shore DelBay High Sal_NEH_RU inbred/CapeShore_2016_Apr_Oct.xlsx",skip=2)
Cape_Shore_2016_Apr_Oct <- Cape_Shore_2016_Apr_Oct[,c(2:3)]
Cape_Shore_2016_Apr_Oct <- Cape_Shore_2016_Apr_Oct %>% drop_na()
colnames(Cape_Shore_2016_Apr_Oct) <- c("Date_Time", "temp_F")
# convert farenheit to celsius
Cape_Shore_2016_Apr_Oct$temp_C <- UnitConv(Cape_Shore_2016_Apr_Oct$temp_F, "F","C")

# 2017 Jul_Sep
Cape_Shore_2017_Jul_Sep <- read_excel("./Environmental_Data/CS Cape Shore DelBay High Sal_NEH_RU inbred/Cape_Shore_2017_Jul_Sep.xlsx", skip=2)
Cape_Shore_2017_Jul_Sep <- Cape_Shore_2017_Jul_Sep[,c(2,4)] 
Cape_Shore_2017_Jul_Sep <-Cape_Shore_2017_Jul_Sep %>% drop_na()
colnames(Cape_Shore_2017_Jul_Sep) <- c("Date_Time", "temp_F")
Cape_Shore_2017_Jul_Sep$temp_C <- UnitConv(Cape_Shore_2017_Jul_Sep$temp_F, "F","C")

# 2017 Sep_Dec
Cape_Shore_2017_Sep_Dec <- read_excel("./Environmental_Data/CS Cape Shore DelBay High Sal_NEH_RU inbred/Cape_Shore_2017_Sep-Dec.xlsx", skip=2)
Cape_Shore_2017_Sep_Dec <- Cape_Shore_2017_Sep_Dec[,c(2,4)]
Cape_Shore_2017_Sep_Dec <- Cape_Shore_2017_Sep_Dec %>% drop_na()
colnames(Cape_Shore_2017_Sep_Dec) <- c("Date_Time", "temp_C")
  
# Combine Cape Shore spreadsheets
Cape_Shore_combined <- rbind(Cape_Shore_2013_2014[,c(1,3)], Cape_Shore_2016_Apr_Oct[,c(1,3)], Cape_Shore_2017_Jul_Sep[,c(1,3)], Cape_Shore_2017_Sep_Dec)
summary(Cape_Shore_combined)

# look up average Cape Shore temperature for Oct. 2017 
Cape_Shore_combined_Oct_2017 <- Cape_Shore_combined[grepl("2017-10-", Cape_Shore_combined$Date_Time),]
summary(Cape_Shore_combined_Oct_2017$temp_C)

#### ITEM 3 Loading in Remaining spreadsheets and check calculations ##### 
# load in remaining spreadsheets for each site or group of sites, 16 populations total

CB_Low_Salinity_LOLA_CLP <- read_csv("./Environmental_Data/Chesapeake Bay Low Salinity/Chesapeake low salinity 605949/CBMIPWQ_R_formatted.csv", col_names =TRUE)
CB_Low_Salinity_LOLA_CLP_subset <- CB_Low_Salinity_LOLA_CLP[,c("DateTimeStamp","Temp_C","Sal_PSU")]
summary(CB_Low_Salinity_LOLA_CLP_subset$Temp_C)
summary(CB_Low_Salinity_LOLA_CLP_subset$Sal_PSU)

# filling in LOLA temperature from 10/2017
library(utils) # use glob2rx package to convert wildcard to regex
CB_Low_Salinity_LOLA_CLP_subset_split <- str_split_fixed(CB_Low_Salinity_LOLA_CLP_subset$DateTimeStamp, " ", 3)
CB_Low_Salinity_LOLA_CLP_subset_split <- str_split_fixed(CB_Low_Salinity_LOLA_CLP_subset_split[,1], "/",3)
CB_Low_Salinity_LOLA_CLP_subset_split_combined <- cbind(CB_Low_Salinity_LOLA_CLP_subset_split, CB_Low_Salinity_LOLA_CLP_subset)
colnames(CB_Low_Salinity_LOLA_CLP_subset_split_combined) <- c("Month","Day","Year","DateTimeStamp","Temp_C","Sal_PSU")
CB_Low_Salinity_LOLA_CLP_subset_split_combined_Oct <- CB_Low_Salinity_LOLA_CLP_subset_split_combined %>% filter(Year == "17") 
CB_Low_Salinity_LOLA_CLP_subset_split_combined_Oct <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_Oct %>% filter(Month == "10")
summary(CB_Low_Salinity_LOLA_CLP_subset_split_combined_Oct$Temp_C)
summary(CB_Low_Salinity_LOLA_CLP_subset_split_combined_Oct$Sal_PSU)

CB_High_Salinity_DEBY_York_Like <- read_csv("./Environmental_Data/ChesapeakeBay High Sal/Chesapeake Bay york like 629420/CBMOCWQ_R_formatted.csv", col_names=TRUE)
CB_High_Salinity_DEBY_York_Like_subset <- CB_High_Salinity_DEBY_York_Like[,c("DateTimeStamp","Temp_C","Sal_PSU")]
summary(CB_High_Salinity_DEBY_York_Like_subset$Temp_C)
summary(CB_High_Salinity_DEBY_York_Like_subset$Sal_PSU)

CB_High_Salinity_HC_VA <- read_csv("./Environmental_Data/ChesapeakeBay High Sal/Chesapeake high salinity 827793/CBVGIWQ_R_formatted.csv", col_names = TRUE)
CB_High_Salinity_HC_VA_subset <- CB_High_Salinity_HC_VA[,c("DateTimeStamp","Temp_C","Sal_PSU")]
summary(CB_High_Salinity_HC_VA_subset$Temp_C)
summary(CB_High_Salinity_HC_VA_subset$Sal_PSU)

HC_Hope_Creek_2012 <- read_excel("./Environmental_Data/HC Hope Creek DelBay Low Sal/URI_Gomez-Chiarra_HopeCreek_2012through2015.xlsx", sheet ="2012")
colnames(HC_Hope_Creek_2012)[5] <- "Sal_ppt"
HC_Hope_Creek_2013 <- read_excel("./Environmental_Data/HC Hope Creek DelBay Low Sal/URI_Gomez-Chiarra_HopeCreek_2012through2015.xlsx", sheet ="2013")
HC_Hope_Creek_2014 <- read_excel("./Environmental_Data/HC Hope Creek DelBay Low Sal/URI_Gomez-Chiarra_HopeCreek_2012through2015.xlsx", sheet ="2014")
colnames(HC_Hope_Creek_2014)[5] <- "Sal_ppt"
colnames(HC_Hope_Creek_2014)[3] <- "TempC"
HC_Hope_Creek_2015 <- read_excel("./Environmental_Data/HC Hope Creek DelBay Low Sal/URI_Gomez-Chiarra_HopeCreek_2012through2015.xlsx", sheet ="2015")
HC_Hope_Creek_Combined <- rbind(HC_Hope_Creek_2012[,c("Date","TempC","Sal_ppt")], HC_Hope_Creek_2013[,c("Date","TempC","Sal_ppt")], HC_Hope_Creek_2014[,c("Date","TempC","Sal_ppt")], HC_Hope_Creek_2015[,c("Date","TempC","Sal_ppt")])
summary(HC_Hope_Creek_Combined$TempC)
summary(HC_Hope_Creek_Combined$Sal_ppt)

LM_Laguna_Madre <- read_csv("./Environmental_Data/LM Laguna Madre/Texas_Aransas_569727/MARABWQ_R_formatted.csv", col_names = TRUE)
LM_Laguna_Madre_subset <- LM_Laguna_Madre[,c("DateTimeStamp","Temp_C","Sal_PSU")]
summary(LM_Laguna_Madre_subset$Temp_C)
summary(LM_Laguna_Madre_subset $Sal_PSU)

#### ITEM 4 Calculate Degree days below 0,15,30 ####
# Split the dates into individual columns for each spreadsheet (day, month and year)
splitenv <- str_split_fixed(env$Daily.Averages, "/", 3)
env <- cbind(env, splitenv)
colnames(env)[9] <- "Month"
colnames(env)[10] <- "Day"
colnames(env)[11] <- "Year"

### CL and SL degree days calculations
# CL Calcasieu Lake and SL Sister Lake temp data only given as monthly averages 
# Going to make a key where each month = a certain number of days
CL_Calcasieu_Lake_temp
SL_Sister_Lake_Caillou_Lake_temp

# Add key for number of days in each month to table
month_day_number_key <- read.csv("./Environmental_Data/Month_day_number_key.csv", header=TRUE)

CL_Calcasieu_Lake_temp_day_number <- merge(CL_Calcasieu_Lake_temp, month_day_number_key, by ="month_nu")
SL_Sister_Lake_Caillou_Lake_temp_day_number <- merge(SL_Sister_Lake_Caillou_Lake_temp, month_day_number_key, by ="month_nu")
write.csv(file="./Environmental_Data/CL_Calcasieu_Lake_temp_day_number.csv",CL_Calcasieu_Lake_temp_day_number)
write.csv(file="./Environment_Data/SL_Sister_Lake_Caillou_Lake_temp_day_number.csv", SL_Sister_Lake_Caillou_Lake_temp_day_number)
unique(CL_Calcasieu_Lake_temp_day_number$year_nu)

#1999 
CL_Calcasieu_Lake_temp_day_number_1999 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "1999")
CL_Calcasieu_Lake_temp_day_number_1999
which(CL_Calcasieu_Lake_temp_day_number_1999$mean_va <=0)
CL_Calcasieu_Lake_temp_day_number_1999_dd_0 <- 0

which(CL_Calcasieu_Lake_temp_day_number_1999$mean_va <=15) # row 10
CL_Calcasieu_Lake_temp_day_number_1999_dd_15 <- CL_Calcasieu_Lake_temp_day_number_1999[10,"Number_Days"] 

which(CL_Calcasieu_Lake_temp_day_number_1999$mean_va <=30) # 1  2  3  4  7  8  9 10
CL_Calcasieu_Lake_temp_day_number_1999_dd_30<- sum(CL_Calcasieu_Lake_temp_day_number_1999[c(1,2,3,4,7,8,9,10),"Number_Days"])

#2000 
CL_Calcasieu_Lake_temp_day_number_2000 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2000")
which(CL_Calcasieu_Lake_temp_day_number_2000$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2000$mean_va <=15) # row 1
CL_Calcasieu_Lake_temp_day_number_2000_dd_15 <- CL_Calcasieu_Lake_temp_day_number_2000[1,"Number_Days"] 

which(CL_Calcasieu_Lake_temp_day_number_2000$mean_va <=30) # 1 2 3 4 5 6
CL_Calcasieu_Lake_temp_day_number_2000_dd_30<- sum(CL_Calcasieu_Lake_temp_day_number_2000[c(1, 2, 3, 4, 5, 6),"Number_Days"])

#2001 
CL_Calcasieu_Lake_temp_day_number_2001 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2001")
which(CL_Calcasieu_Lake_temp_day_number_2001$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2001$mean_va <=15)
CL_Calcasieu_Lake_temp_day_number_2001_dd_15 <- 0

which(CL_Calcasieu_Lake_temp_day_number_2001$mean_va <=30) # 1 2 3 
CL_Calcasieu_Lake_temp_day_number_2001_dd_30<- sum(CL_Calcasieu_Lake_temp_day_number_2001[c(1, 2, 3),"Number_Days"])

#2002
CL_Calcasieu_Lake_temp_day_number_2002 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2002")
which(CL_Calcasieu_Lake_temp_day_number_2002$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2002$mean_va <=15) #1,8
CL_Calcasieu_Lake_temp_day_number_2002_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2002[c(1, 8),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2002$mean_va <=30) # 1 2 3 4 5 7 8 
CL_Calcasieu_Lake_temp_day_number_2002_dd_30<- sum(CL_Calcasieu_Lake_temp_day_number_2002[c(1, 2, 3, 4 ,5, 7, 8),"Number_Days"])

# 2003
CL_Calcasieu_Lake_temp_day_number_2003 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2003")
which(CL_Calcasieu_Lake_temp_day_number_2003$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2003$mean_va <=15) #1,2, 8
CL_Calcasieu_Lake_temp_day_number_2003_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2003[c(1,2, 8),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2003$mean_va <=30) # 1 2 3 4 5 6 7 8
CL_Calcasieu_Lake_temp_day_number_2003_dd_30<- sum(CL_Calcasieu_Lake_temp_day_number_2003[c(1, 2, 3, 4, 5, 6, 7, 8),"Number_Days"])

# 2004
CL_Calcasieu_Lake_temp_day_number_2004 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2004")
which(CL_Calcasieu_Lake_temp_day_number_2004$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2004$mean_va <=15) #1,2, 10
CL_Calcasieu_Lake_temp_day_number_2004_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2004[c(1,2, 10),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2004$mean_va <=30) # 1 2 3 4 5 6 7 8
CL_Calcasieu_Lake_temp_day_number_2004_dd_30<- sum(CL_Calcasieu_Lake_temp_day_number_2004[c(1,2,3,4,5,6,7,8,9, 10),"Number_Days"])

#2005 
CL_Calcasieu_Lake_temp_day_number_2005 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2005")
which(CL_Calcasieu_Lake_temp_day_number_2005$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2005$mean_va <=15) #1,2
CL_Calcasieu_Lake_temp_day_number_2005_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2005[c(1,2),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2005$mean_va <=30) # 1 2 3 4 
CL_Calcasieu_Lake_temp_day_number_2005_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2005[c(1,2,3,4),"Number_Days"])

# 2006
CL_Calcasieu_Lake_temp_day_number_2006 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2006")
which(CL_Calcasieu_Lake_temp_day_number_2006$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2006$mean_va <=15) #6
CL_Calcasieu_Lake_temp_day_number_2006_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2006[c(6),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2006$mean_va <=30) # 1 2 3 4 5 6 
CL_Calcasieu_Lake_temp_day_number_2006_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2006[c(1, 2, 3, 4, 5, 6),"Number_Days"])

#2007
CL_Calcasieu_Lake_temp_day_number_2007 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2007")
which(CL_Calcasieu_Lake_temp_day_number_2007$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2007$mean_va <=15) #1,2
CL_Calcasieu_Lake_temp_day_number_2007_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2007[c(1,2),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2007$mean_va <=30) # 1 2 3 4 5 6 7 8
CL_Calcasieu_Lake_temp_day_number_2007_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2007[c(1, 2, 3, 4, 5, 6, 7, 8),"Number_Days"])

#2008
CL_Calcasieu_Lake_temp_day_number_2008 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2008")
which(CL_Calcasieu_Lake_temp_day_number_2008$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2008$mean_va <=15) #10
CL_Calcasieu_Lake_temp_day_number_2008_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2008[c(10),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2008$mean_va <=30) #  1  2  3  4  5  6  7  8  9 10
CL_Calcasieu_Lake_temp_day_number_2008_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2008[c(1, 2,3,4,5,6,7,8,9,10),"Number_Days"])

#2009
CL_Calcasieu_Lake_temp_day_number_2009 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2009")
which(CL_Calcasieu_Lake_temp_day_number_2009$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2009$mean_va <=15) #0
CL_Calcasieu_Lake_temp_day_number_2009_dd_15 <- 0

which(CL_Calcasieu_Lake_temp_day_number_2009$mean_va <=30) #  1 2 3
CL_Calcasieu_Lake_temp_day_number_2009_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2009[c(1, 2,3),"Number_Days"])

#2010
CL_Calcasieu_Lake_temp_day_number_2010 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2010")
which(CL_Calcasieu_Lake_temp_day_number_2010$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2010$mean_va <=15) #9
CL_Calcasieu_Lake_temp_day_number_2010_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2010[c(9),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2010$mean_va <=30) #  1 2 6 7 8 9
CL_Calcasieu_Lake_temp_day_number_2010_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2010[c(1, 2, 6, 7, 8, 9),"Number_Days"])

#2011
CL_Calcasieu_Lake_temp_day_number_2011 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2011")
which(CL_Calcasieu_Lake_temp_day_number_2011$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2011$mean_va <=15) #1,2,9
CL_Calcasieu_Lake_temp_day_number_2011_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2011[c(1,2,9),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2011$mean_va <=30) #  1 2 3 4 5 6 8 9
CL_Calcasieu_Lake_temp_day_number_2011_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2011[c(1, 2, 3, 4, 5, 6, 8, 9),"Number_Days"])

#2012
CL_Calcasieu_Lake_temp_day_number_2012 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2012")
which(CL_Calcasieu_Lake_temp_day_number_2012$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2012$mean_va <=15) 
CL_Calcasieu_Lake_temp_day_number_2012_dd_15 <- 0

which(CL_Calcasieu_Lake_temp_day_number_2012$mean_va <=30) #  1  2  3  4  5  6  7  9 10
CL_Calcasieu_Lake_temp_day_number_2012_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2012[c(1, 2,3,4,5,6,7,9,10),"Number_Days"])

#2013
CL_Calcasieu_Lake_temp_day_number_2013 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2013")
which(CL_Calcasieu_Lake_temp_day_number_2013$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2013$mean_va <=15) #1
CL_Calcasieu_Lake_temp_day_number_2013_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2013[c(1),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2013$mean_va <=30) #  1 2 3 4 5 7 9
CL_Calcasieu_Lake_temp_day_number_2013_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2013[c(1, 2, 3, 4, 5, 7, 9),"Number_Days"])

#2014
CL_Calcasieu_Lake_temp_day_number_2014 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2014")
which(CL_Calcasieu_Lake_temp_day_number_2014$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2014$mean_va <=15) #1,2
CL_Calcasieu_Lake_temp_day_number_2014_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2014[c(1,2),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2014$mean_va <=30) # 1  2  3  4  7  8  9 10
CL_Calcasieu_Lake_temp_day_number_2014_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2014[c(1,  2,  3,  4,  7,  8,  9, 10),"Number_Days"])

#2015
CL_Calcasieu_Lake_temp_day_number_2015 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2015")
which(CL_Calcasieu_Lake_temp_day_number_2015$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2015$mean_va <=15) #1,2
CL_Calcasieu_Lake_temp_day_number_2015_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2015[c(1,2),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2015$mean_va <=30) # 1  2  3  4  5  6  8  9 10 11
CL_Calcasieu_Lake_temp_day_number_2015_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2015[c(1,  2,  3,  4,  5,  6,  8,  9, 10, 11),"Number_Days"])

#2016
CL_Calcasieu_Lake_temp_day_number_2016 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2016")
which(CL_Calcasieu_Lake_temp_day_number_2016$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2016$mean_va <=15) #1,9
CL_Calcasieu_Lake_temp_day_number_2016_dd_15 <- sum(CL_Calcasieu_Lake_temp_day_number_2016[c(1,9),"Number_Days"])

which(CL_Calcasieu_Lake_temp_day_number_2016$mean_va <=30) #  1 2 3 4 5 7 8 9
CL_Calcasieu_Lake_temp_day_number_2016_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2016[c( 1, 2, 3, 4, 5, 7, 8, 9),"Number_Days"])

#2017
CL_Calcasieu_Lake_temp_day_number_2017 <- CL_Calcasieu_Lake_temp_day_number %>%
  filter(year_nu == "2017")
which(CL_Calcasieu_Lake_temp_day_number_2017$mean_va <=0)

which(CL_Calcasieu_Lake_temp_day_number_2017$mean_va <=15) #0
CL_Calcasieu_Lake_temp_day_number_2017_dd_15 <- 0

which(CL_Calcasieu_Lake_temp_day_number_2017$mean_va <=30) #  1 2 3 4 5 6 7 8
CL_Calcasieu_Lake_temp_day_number_2017_dd_30 <- sum(CL_Calcasieu_Lake_temp_day_number_2017[c(1, 2, 3, 4, 5, 6, 7, 8),"Number_Days"])

#2018 (only 1 reading, not going to be used!)

# Find average degree days over each year data collected
CL_dd_0_average <- sum(CL_Calcasieu_Lake_temp_day_number_1999_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2000_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2001_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2002_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2003_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2004_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2005_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2006_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2007_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2008_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2009_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2010_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2011_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2012_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2013_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2014_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2015_dd_15 , 
                        CL_Calcasieu_Lake_temp_day_number_2016_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2017_dd_15)/19

# no sites had a degree day below freezing

CL_dd_15_average <- sum(CL_Calcasieu_Lake_temp_day_number_1999_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2000_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2001_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2002_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2003_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2004_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2005_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2006_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2007_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2008_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2009_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2010_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2011_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2012_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2013_dd_15, 
                        CL_Calcasieu_Lake_temp_day_number_2014_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2015_dd_15 , 
                        CL_Calcasieu_Lake_temp_day_number_2016_dd_15,
                        CL_Calcasieu_Lake_temp_day_number_2017_dd_15)/19

CL_dd_30_average <- sum(CL_Calcasieu_Lake_temp_day_number_1999_dd_30, CL_Calcasieu_Lake_temp_day_number_2000_dd_30, 
CL_Calcasieu_Lake_temp_day_number_2001_dd_30,
CL_Calcasieu_Lake_temp_day_number_2002_dd_30,  
CL_Calcasieu_Lake_temp_day_number_2003_dd_30, CL_Calcasieu_Lake_temp_day_number_2004_dd_30,
CL_Calcasieu_Lake_temp_day_number_2005_dd_30, CL_Calcasieu_Lake_temp_day_number_2006_dd_30,
CL_Calcasieu_Lake_temp_day_number_2007_dd_30, CL_Calcasieu_Lake_temp_day_number_2008_dd_30,
CL_Calcasieu_Lake_temp_day_number_2009_dd_30, CL_Calcasieu_Lake_temp_day_number_2010_dd_30,
CL_Calcasieu_Lake_temp_day_number_2011_dd_30,CL_Calcasieu_Lake_temp_day_number_2012_dd_30,
CL_Calcasieu_Lake_temp_day_number_2013_dd_30, CL_Calcasieu_Lake_temp_day_number_2014_dd_30,
CL_Calcasieu_Lake_temp_day_number_2015_dd_30, CL_Calcasieu_Lake_temp_day_number_2016_dd_30, 
CL_Calcasieu_Lake_temp_day_number_2017_dd_30)/19

### SL Caillou Lake
SL_Sister_Lake_Caillou_Lake_temp_day_number
sort(SL_Sister_Lake_Caillou_Lake_temp_day_number$year_nu)

#1997
SL_Sister_Lake_Caillou_Lake_temp_day_number_1997 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "1997")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997$mean_va <=15)  #5
SL_Sister_Lake_Caillou_Lake_temp_day_number_1997_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997[c(5),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997$mean_va <=30) #1 2 3 4 5
SL_Sister_Lake_Caillou_Lake_temp_day_number_1997_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997[c(1, 2, 3, 4, 5),"Number_Days"])

#1998
SL_Sister_Lake_Caillou_Lake_temp_day_number_1998 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "1998")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1998$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1998$mean_va <=15)  #1
SL_Sister_Lake_Caillou_Lake_temp_day_number_1998_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1998[c(1),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1998$mean_va <=30) #1  2  3  4  5  6  8  9 10
SL_Sister_Lake_Caillou_Lake_temp_day_number_1998_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1998[c(1,  2,  3,  4,  5,  6,  8,  9, 10),"Number_Days"])

#1999
SL_Sister_Lake_Caillou_Lake_temp_day_number_1999 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "1999")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1999$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1999$mean_va <=15)  #1,12
SL_Sister_Lake_Caillou_Lake_temp_day_number_1999_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1999[c(1,12),"Number_Days"])
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_1999$mean_va <=30) #1  2  3  4  5  6  7  9 10 11 12
SL_Sister_Lake_Caillou_Lake_temp_day_number_1999_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1999[c(1,  2,  3,  4,  5,  6,  7,  9, 10, 11, 12),"Number_Days"])

#2000
SL_Sister_Lake_Caillou_Lake_temp_day_number_2000 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2000")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2000$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2000$mean_va <=15)  #1
SL_Sister_Lake_Caillou_Lake_temp_day_number_2000_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2000[c(1),"Number_Days"])
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2000$mean_va <=30) #1 2 3 4 5 6
SL_Sister_Lake_Caillou_Lake_temp_day_number_2000_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2000[c(1, 2, 3, 4, 5, 6),"Number_Days"])

#2001
SL_Sister_Lake_Caillou_Lake_temp_day_number_2001 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2001")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2001$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2001$mean_va <=15)  
SL_Sister_Lake_Caillou_Lake_temp_day_number_2001_dd_15 <- 0
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2001$mean_va <=30) #1 2 3
SL_Sister_Lake_Caillou_Lake_temp_day_number_2001_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2001[c(1, 2, 3),"Number_Days"])

#2002
SL_Sister_Lake_Caillou_Lake_temp_day_number_2002 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2002")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2002$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2002$mean_va <=15)  #1,2
SL_Sister_Lake_Caillou_Lake_temp_day_number_2002_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2002[c(1,2),"Number_Days"])
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2002$mean_va <=30) #1 2 3 4 5 7 8
SL_Sister_Lake_Caillou_Lake_temp_day_number_2002_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2002[c(1, 2, 3, 4, 5, 7, 8),"Number_Days"])

#2003
SL_Sister_Lake_Caillou_Lake_temp_day_number_2003 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2003")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2003$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2003$mean_va <=15)  #1,11
SL_Sister_Lake_Caillou_Lake_temp_day_number_2003_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2003[c(1,11),"Number_Days"])
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2003$mean_va <=30) #1  2  3  4  5  6  7  9 10 11
SL_Sister_Lake_Caillou_Lake_temp_day_number_2003_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2003[c(1 , 2,  3,  4,  5,  6,  7,  9, 10, 11),"Number_Days"])

#2004
SL_Sister_Lake_Caillou_Lake_temp_day_number_2004 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2004")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2004$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2004$mean_va <=15)  #1,2,12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2004_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2004[c(1,2,12),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2004$mean_va <=30) #1  2  3  4  5  6  8  9 10 11 12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2004_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2004[c(1,  2,  3,  4,  5,  6,  8,  9, 10, 11, 12),"Number_Days"])

#2005
SL_Sister_Lake_Caillou_Lake_temp_day_number_2005 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2005")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2005$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2005$mean_va <=15)  #11
SL_Sister_Lake_Caillou_Lake_temp_day_number_2005_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2005[c(11),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2005$mean_va <=30) #1  2  3  4  5  6  9 10 11
SL_Sister_Lake_Caillou_Lake_temp_day_number_2005_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2005[c(1,  2,  3,  4,  5,  6,  9, 10, 11),"Number_Days"])

#2006
SL_Sister_Lake_Caillou_Lake_temp_day_number_2006 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2006")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2006$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2006$mean_va <=15)  #7
SL_Sister_Lake_Caillou_Lake_temp_day_number_2006_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2006[c(7),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2006$mean_va <=30) #1 2 5 6 7
SL_Sister_Lake_Caillou_Lake_temp_day_number_2006_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2006[c(1, 2, 5, 6, 7),"Number_Days"])

#2007
SL_Sister_Lake_Caillou_Lake_temp_day_number_2007 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2007")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2007$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2007$mean_va <=15)  #1,2
SL_Sister_Lake_Caillou_Lake_temp_day_number_2007_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2007[c(1,2),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2007$mean_va <=30) #1 2 3 4 5 6 7 9
SL_Sister_Lake_Caillou_Lake_temp_day_number_2007_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2007[c(1, 2, 3, 4, 5, 6, 7, 9),"Number_Days"])

#2008
SL_Sister_Lake_Caillou_Lake_temp_day_number_2008 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2008")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2008$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2008$mean_va <=15)  #10
SL_Sister_Lake_Caillou_Lake_temp_day_number_2008_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2008[c(10),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2008$mean_va <=30) #1  2  3  4  7  8  9 10
SL_Sister_Lake_Caillou_Lake_temp_day_number_2008_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2008[c(1,2,3,4,7,8, 9,10),"Number_Days"])

#2009
SL_Sister_Lake_Caillou_Lake_temp_day_number_2009 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2009")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2009$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2009$mean_va <=15)  #11
SL_Sister_Lake_Caillou_Lake_temp_day_number_2009_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2009[c(11),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2009$mean_va <=30) #1  2  3  4  8  9 10 11
SL_Sister_Lake_Caillou_Lake_temp_day_number_2009_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2009[c(1, 2,3,4,8,9,10,11),"Number_Days"])

#2010
SL_Sister_Lake_Caillou_Lake_temp_day_number_2010 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2010")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2010$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2010$mean_va <=15)  #1,2,12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2010_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2010[c(1,2,12),"Number_Days"])
w
hich(SL_Sister_Lake_Caillou_Lake_temp_day_number_2010$mean_va <=30) #1  2  3  4  5  9 10 11 12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2010_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2010[c(1, 2,3,4,5,9,10,11,12),"Number_Days"])

#2011
SL_Sister_Lake_Caillou_Lake_temp_day_number_2011 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2011")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2011$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2011$mean_va <=15)  #1,2,12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2011_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2011[c(1,2,12),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2011$mean_va <=30) #1  2  3  4  5  9 10 11 12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2011_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2011[c(1, 2,3,4,5,9,10,11,12),"Number_Days"])

#2012
SL_Sister_Lake_Caillou_Lake_temp_day_number_2012 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2012")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2012$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2012$mean_va <=15)
SL_Sister_Lake_Caillou_Lake_temp_day_number_2012_dd_15 <- 0
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2012$mean_va <=30) #1  2  3  4  5  6  8  9 10 11 12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2012_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2012[c(1,2,3,4,5,6,8,9,10,11,12),"Number_Days"])

#2013
SL_Sister_Lake_Caillou_Lake_temp_day_number_2013 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2013")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2013$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2013$mean_va <=15) #1,10
SL_Sister_Lake_Caillou_Lake_temp_day_number_2013_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2013[c(1,10),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2013$mean_va <=30) #1  2  3  4  5  7  8  9 10
SL_Sister_Lake_Caillou_Lake_temp_day_number_2013_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2013[c(1, 2,3,4,5,7,8,9,10),"Number_Days"])

#2014
SL_Sister_Lake_Caillou_Lake_temp_day_number_2014 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2014")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2014$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2014$mean_va <=15) #1,2
SL_Sister_Lake_Caillou_Lake_temp_day_number_2014_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2014[c(1,2),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2014$mean_va <=30) #1  2  3  4  5  6  7  9 10 11 12
SL_Sister_Lake_Caillou_Lake_temp_day_number_2014_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2014[c(1,2,3,4,5,6,7,9,10,11,12),"Number_Days"])

#2015
SL_Sister_Lake_Caillou_Lake_temp_day_number_2015 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2015")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2015$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2015$mean_va <=15) #1,2
SL_Sister_Lake_Caillou_Lake_temp_day_number_2015_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2015[c(1,2),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2015$mean_va <=30) #1  2  3  4  5  8  9 10
SL_Sister_Lake_Caillou_Lake_temp_day_number_2015_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2015[c(1,2,3,4,5,8,9, 10),"Number_Days"])

#2016
SL_Sister_Lake_Caillou_Lake_temp_day_number_2016 <- SL_Sister_Lake_Caillou_Lake_temp_day_number %>% 
  filter(year_nu == "2016")
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2016$mean_va <= 0)
which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2016$mean_va <=15) #1
SL_Sister_Lake_Caillou_Lake_temp_day_number_2016_dd_15 <-  sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2016[c(1),"Number_Days"])

which(SL_Sister_Lake_Caillou_Lake_temp_day_number_2016$mean_va <=30) #1  2  3  4  6  7  8  9 10
SL_Sister_Lake_Caillou_Lake_temp_day_number_2016_dd_30 <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_2016[c(1,  2,  3,  4,  6,  7,  8,  9, 10),"Number_Days"])

#2017 and 2018 have very limited data
#SL_dd_0 no months are below zero

SL_dd_15_average <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_1998_dd_15, 
                SL_Sister_Lake_Caillou_Lake_temp_day_number_1999_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2000_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2001_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2002_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2003_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2004_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2005_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2006_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2007_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2008_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2009_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2010_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2011_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2012_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2013_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2014_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2015_dd_15,
                SL_Sister_Lake_Caillou_Lake_temp_day_number_2016_dd_15)/20

SL_dd_30_average <- sum(SL_Sister_Lake_Caillou_Lake_temp_day_number_1997_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_1998_dd_30, 
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_1999_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2000_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2001_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2002_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2003_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2004_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2005_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2006_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2007_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2008_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2009_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2010_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2011_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2012_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2013_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2014_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2015_dd_30,
                        SL_Sister_Lake_Caillou_Lake_temp_day_number_2016_dd_30)/20


### Hog Island and Sheepscot data given as one reading per month
Hog_island_combined_removed_datesplit <- str_split_fixed(Hog_island_combined_removed$Date," ", 2)
Hog_island_combined_removed_datesplit <- str_split_fixed(Hog_island_combined_removed_datesplit[,1],"-", 3)
Hog_island_combined_removed_datesplit_combined <- cbind(Hog_island_combined_removed_datesplit, Hog_island_combined_removed)
colnames(Hog_island_combined_removed_datesplit_combined) <- c("Year","Month","Day","Station","Date","Tide","Temp","Sal","LATITUDE","LONGITUDE")
head(Hog_island_combined_removed_datesplit_combined)

head(Sherman_Marsh_Sheepscot_removed)
Sherman_Marsh_Sheepscot_removed_datesplit <- str_split_fixed(Sherman_Marsh_Sheepscot_removed$EFFORT_START_DATE, "-",3)
Sherman_Marsh_Sheepscot_removed_combined<- cbind(Sherman_Marsh_Sheepscot_removed_datesplit,Sherman_Marsh_Sheepscot_removed)
colnames(Sherman_Marsh_Sheepscot_removed_combined) <- c("Year","Month","Day","LOCATION_ID","EFFORT_START_DATE","TIDE_STAGE","TEMP_C",
                                                        "SALINITY_PCT","LATITUDE_DECIMAL","LONGITUDE_DECIMAL")

# Add day number per month column to both Hog Island and Sheepscot
class(month_day_number_key$month_nu)
month_day_number_key_hog_island <- month_day_number_key
colnames(month_day_number_key_hog_island) <- c("Month", "Number_Days")
Hog_island_combined_removed_datesplit_combined$Month <- as.integer(Hog_island_combined_removed_datesplit_combined$Month)
Hog_island_combined_removed_datesplit_combined_days <- merge(Hog_island_combined_removed_datesplit_combined, month_day_number_key_hog_island, by ="Month")

class(Sherman_Marsh_Sheepscot_removed_combined$Month) #class is factor need to change to integer
Sherman_Marsh_Sheepscot_removed_combined$Month <- as.integer(Sherman_Marsh_Sheepscot_removed_combined$Month)
Sherman_Marsh_Sheepscot_removed_combined_days <- merge(Sherman_Marsh_Sheepscot_removed_combined, month_day_number_key_hog_island, by ="Month")
Sherman_Marsh_Sheepscot_removed_combined_days <- na.omit(Sherman_Marsh_Sheepscot_removed_combined_days)

# Hog Island degree days per year
levels(Hog_island_combined_removed_datesplit_combined_days$Year)

#1990 has only two days sampled
Hog_island_combined_removed_datesplit_combined_days_1991 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1991")
Hog_island_combined_removed_datesplit_combined_days_1992 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1992")
Hog_island_combined_removed_datesplit_combined_days_1993 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1993")
Hog_island_combined_removed_datesplit_combined_days_1994 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1994")
Hog_island_combined_removed_datesplit_combined_days_1995 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1995")
Hog_island_combined_removed_datesplit_combined_days_1996<- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1996")
Hog_island_combined_removed_datesplit_combined_days_1997 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1997")
Hog_island_combined_removed_datesplit_combined_days_1998 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1998")
Hog_island_combined_removed_datesplit_combined_days_1999 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "1999")
Hog_island_combined_removed_datesplit_combined_days_2000 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2000")
Hog_island_combined_removed_datesplit_combined_days_2001 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2001")
Hog_island_combined_removed_datesplit_combined_days_2002 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2002")
Hog_island_combined_removed_datesplit_combined_days_2003 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2003")
Hog_island_combined_removed_datesplit_combined_days_2004 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2004")
Hog_island_combined_removed_datesplit_combined_days_2005 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2005")
Hog_island_combined_removed_datesplit_combined_days_2006 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2006")
Hog_island_combined_removed_datesplit_combined_days_2007 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2007")
Hog_island_combined_removed_datesplit_combined_days_2008 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2008")
Hog_island_combined_removed_datesplit_combined_days_2009 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2009")
Hog_island_combined_removed_datesplit_combined_days_2010 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2010")
Hog_island_combined_removed_datesplit_combined_days_2011 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2011")
Hog_island_combined_removed_datesplit_combined_days_2012 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2012")
Hog_island_combined_removed_datesplit_combined_days_2013 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2013")
Hog_island_combined_removed_datesplit_combined_days_2014 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2014")
Hog_island_combined_removed_datesplit_combined_days_2015 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2015")
Hog_island_combined_removed_datesplit_combined_days_2016 <- Hog_island_combined_removed_datesplit_combined_days %>% 
  filter(Year == "2016")

# HI dd_0

which(Hog_island_combined_removed_datesplit_combined_days_1991$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_1991_dd_0 <- 0
  
which(Hog_island_combined_removed_datesplit_combined_days_1992$Temp <=0) # 1  2  5 26
Hog_island_combined_removed_datesplit_combined_days_1992_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_1992[c(1,  2,  5, 26),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1993$Temp <=0) #1  2 24
Hog_island_combined_removed_datesplit_combined_days_1993_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_1993[c(1,  2, 24),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1994$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_1994_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_1995$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_1995_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_1996$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_1996_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_1997$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_1997_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_1998$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_1998_dd_0 <-0

which(Hog_island_combined_removed_datesplit_combined_days_1999$Temp <=0) # 1,2
Hog_island_combined_removed_datesplit_combined_days_1999_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_1999[c(1,  2),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2000$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2000_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_2001$Temp <=0) # 1,2 
Hog_island_combined_removed_datesplit_combined_days_2001_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2001[c(1,  2),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2002$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2002_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_2003$Temp <=0) # 1,2,3,4
Hog_island_combined_removed_datesplit_combined_days_2003_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2003[c(1,2,3,4),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2004$Temp <=0) # 4  5 33 34 35 36
Hog_island_combined_removed_datesplit_combined_days_2004_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2004[c(4,  5, 33, 34, 35, 36),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2005$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2005_dd_0 <-0

which(Hog_island_combined_removed_datesplit_combined_days_2006$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2006_dd_0 <-0

which(Hog_island_combined_removed_datesplit_combined_days_2007$Temp <=0)# 21
Hog_island_combined_removed_datesplit_combined_days_2007_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2007[c(21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2008$Temp <=0) # 1,4,6
Hog_island_combined_removed_datesplit_combined_days_2008_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2008[c(1,4,6),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2009$Temp <=0) # 1,2,3
Hog_island_combined_removed_datesplit_combined_days_2009_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2009[c(1,2,3),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2010$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2010_dd_0 <-0

which(Hog_island_combined_removed_datesplit_combined_days_2011$Temp <=0) #1
Hog_island_combined_removed_datesplit_combined_days_2011_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2011[c(1),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2012$Temp <=0) # 3
Hog_island_combined_removed_datesplit_combined_days_2012_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2012[c(3),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2013$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2013_dd_0 <- 0

which(Hog_island_combined_removed_datesplit_combined_days_2014$Temp <=0) # 1
Hog_island_combined_removed_datesplit_combined_days_2014_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2014[c(1),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2015$Temp <=0) # 1,2
Hog_island_combined_removed_datesplit_combined_days_2015_dd_0 <- sum(Hog_island_combined_removed_datesplit_combined_days_2015[c(1,2),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2016$Temp <=0)
Hog_island_combined_removed_datesplit_combined_days_2016_dd_0 <-0

HI_dd_0_average <- sum(Hog_island_combined_removed_datesplit_combined_days_1991_dd_0,
  Hog_island_combined_removed_datesplit_combined_days_1992_dd_0, 
                       Hog_island_combined_removed_datesplit_combined_days_1993_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_1994_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_1995_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_1996_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_1997_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_1998_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_1999_dd_0, 
  Hog_island_combined_removed_datesplit_combined_days_2000_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2001_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2002_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2003_dd_0, 
                       Hog_island_combined_removed_datesplit_combined_days_2004_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2005_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2006_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2007_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2008_dd_0, 
                       Hog_island_combined_removed_datesplit_combined_days_2009_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2010_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2011_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2012_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2013_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2014_dd_0, 
                       Hog_island_combined_removed_datesplit_combined_days_2015_dd_0,
                       Hog_island_combined_removed_datesplit_combined_days_2016_dd_0)/26

# HI dd_15
which(Hog_island_combined_removed_datesplit_combined_days_1991$Temp <=15) # 1  2  3 10 11
Hog_island_combined_removed_datesplit_combined_days_1991_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1991[c(1,  2,  3, 10, 11),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1992$Temp <=15) # 1  2  3  4  5  6  7  8  9 20 21 22 23 24 25 26
Hog_island_combined_removed_datesplit_combined_days_1992_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1992[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 20, 21, 22, 23, 24, 25, 26),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1993$Temp <=15) #1  2  3  4  5 14 15 17 18 19 20 21 22 23 24 25
Hog_island_combined_removed_datesplit_combined_days_1993_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1993[c(1,  2,  3,  4,  5, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1994$Temp <=15)
Hog_island_combined_removed_datesplit_combined_days_1994_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1994[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1995$Temp <=15) #  1  2  3  4 17 18 19
Hog_island_combined_removed_datesplit_combined_days_1995_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1995[c( 1,  2,  3,  4, 17, 18, 19),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1996$Temp <=15) # 1  2  3  4  5  6 19 20 21
Hog_island_combined_removed_datesplit_combined_days_1996_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1996[c(1,  2,  3,  4,  5,  6, 19, 20, 21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1997$Temp <=15) # 1  2  3  4  5  6 18 20
Hog_island_combined_removed_datesplit_combined_days_1997_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1997[c(1,  2,  3,  4,  5,  6, 18, 20),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1998$Temp <=15) # 1  2  3  4  5  6 16 17 18 19 20 21
Hog_island_combined_removed_datesplit_combined_days_1998_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1998[c(1,  2,  3,  4,  5,  6, 16, 17, 18, 19, 20, 21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1999$Temp <=15) # 1  2  3  4  5  6  7  8 21 22 23 24 25 26 27 28 29
Hog_island_combined_removed_datesplit_combined_days_1999_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_1999[c(1,  2,  3,  4,  5,  6,  7,  8, 21, 22, 23, 24, 25, 26, 27, 28, 29),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2000$Temp <=15) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 24 25 26 27 28 29
Hog_island_combined_removed_datesplit_combined_days_2000_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2000[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 24, 25, 26, 27, 28, 29),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2001$Temp <=15) # 1,2 
Hog_island_combined_removed_datesplit_combined_days_2001_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2001[c(1,  2),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2002$Temp <=15)# 1  2  3  4  5  6  7  8  9 10 27 28 29 30 31
Hog_island_combined_removed_datesplit_combined_days_2002_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2002[c(1,  2,  3,  4,  5,  6,  7, 8,  9, 10, 27, 28, 29, 30, 31),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2003$Temp <=15) # 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 20, 21, 23, 24, 25, 26, 27, 28, 29
Hog_island_combined_removed_datesplit_combined_days_2003_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2003[c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 20, 21, 23, 24, 25, 26, 27, 28, 29),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2004$Temp <=15) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 24 25 26 27 28 29 30 31 32 33 34 35 36
Hog_island_combined_removed_datesplit_combined_days_2004_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2004[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2005$Temp <=15) # 1  2  3  4  6 12 13 14 15 16 17 18 19 20
Hog_island_combined_removed_datesplit_combined_days_2005_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2005[c(1,  2,  3,  4,  6, 12, 13, 14, 15, 16, 17, 18, 19, 20),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2006$Temp <=15) # 1  2  3  4  5  6 12 13 14 15 16 17 18 19
Hog_island_combined_removed_datesplit_combined_days_2006_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2006[c(1,  2,  3,  4,  5,  6, 12, 13, 14, 15, 16, 17, 18, 19),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2007$Temp <=15) # 1  2  3  4  5  6  7  8  9 10 11 12 13 18 19 20 21
Hog_island_combined_removed_datesplit_combined_days_2007_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2007[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 18, 19, 20, 21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2008$Temp <=15) #  1  2  3  4  5  6  7  8  9 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2008_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2008[c( 1,  2,  3,  4 , 5,  6 , 7,  8 , 9 ,16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2009$Temp <=15) # 1  2  3  4  5  6  7  8  9 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2009_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2009[c(1 , 2,  3,  4,  5,  6,  7,  8,  9, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2010$Temp <=15) # 1  2  3  4  5 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2010_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2010[c(1,  2,  3,  4,  5, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2011$Temp <=15)  # 1  2  3  4  5  6 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2011_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2011[c(1,  2,  3,  4,  5,  6, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2012$Temp <=15) #1  2  3  4  5  6  7  8  9 16 17
Hog_island_combined_removed_datesplit_combined_days_2012_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2012[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 16, 17),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2013$Temp <=15) # 1  2  3  4 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2013_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2013[c(1,  2,  3,  4, 16, 17, 18 ),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2014$Temp <=15) # 1  2  3  4  5  6  7  8 18
Hog_island_combined_removed_datesplit_combined_days_2014_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2014[c(1,  2,  3,  4,  5,  6,  7,  8, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2015$Temp <=15) # 1  2  3  4  5  6 13 14 16 17 18 19
Hog_island_combined_removed_datesplit_combined_days_2015_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2015[c(1,  2,  3,  4,  5,  6, 13, 14, 16, 17, 18, 19),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2016$Temp <=15) # 1  2  3  4  5  6  7  8  9 10 11 12
Hog_island_combined_removed_datesplit_combined_days_2016_dd_15 <- sum(Hog_island_combined_removed_datesplit_combined_days_2016[c(1 , 2,  3,  4,  5,  6,  7,  8  ,9, 10, 11, 12),"Number_Days"])

HI_dd_15_average <- sum(Hog_island_combined_removed_datesplit_combined_days_1991_dd_15, Hog_island_combined_removed_datesplit_combined_days_1992_dd_15, 
Hog_island_combined_removed_datesplit_combined_days_1993_dd_15, Hog_island_combined_removed_datesplit_combined_days_1994_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_1995_dd_15 ,Hog_island_combined_removed_datesplit_combined_days_1996_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_1997_dd_15 ,Hog_island_combined_removed_datesplit_combined_days_1998_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_1999_dd_15 ,Hog_island_combined_removed_datesplit_combined_days_2000_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2001_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2002_dd_15,
Hog_island_combined_removed_datesplit_combined_days_2003_dd_15,
Hog_island_combined_removed_datesplit_combined_days_2004_dd_15,
Hog_island_combined_removed_datesplit_combined_days_2005_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2006_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2007_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2008_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2009_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2010_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2011_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2012_dd_15,
Hog_island_combined_removed_datesplit_combined_days_2013_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2014_dd_15,
Hog_island_combined_removed_datesplit_combined_days_2015_dd_15 ,
Hog_island_combined_removed_datesplit_combined_days_2016_dd_15)/26

## HI_dd_30 
which(Hog_island_combined_removed_datesplit_combined_days_1991$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11
Hog_island_combined_removed_datesplit_combined_days_1991_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1991[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1992$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
Hog_island_combined_removed_datesplit_combined_days_1992_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1992[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1993$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
Hog_island_combined_removed_datesplit_combined_days_1993_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1993[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1994$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Hog_island_combined_removed_datesplit_combined_days_1994_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1994[c(1, 2,  3 , 4,  5 , 6,  7,  8 , 9 ,10, 11, 12, 13, 14, 15, 16 ,17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31, 32),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1995$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
Hog_island_combined_removed_datesplit_combined_days_1995_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1995[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1996$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
Hog_island_combined_removed_datesplit_combined_days_1996_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1996[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1997$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
Hog_island_combined_removed_datesplit_combined_days_1997_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1997[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1998$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
Hog_island_combined_removed_datesplit_combined_days_1998_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1998[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_1999$Temp <=30) #  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
Hog_island_combined_removed_datesplit_combined_days_1999_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_1999[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2000$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
Hog_island_combined_removed_datesplit_combined_days_2000_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2000[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2001$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
Hog_island_combined_removed_datesplit_combined_days_2001_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2001[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2002$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 
Hog_island_combined_removed_datesplit_combined_days_2002_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2002[c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2003$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
Hog_island_combined_removed_datesplit_combined_days_2003_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2003[c( 1,  2,  3,  4 , 5,  6,  7 , 8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2004$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
Hog_island_combined_removed_datesplit_combined_days_2004_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2004[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                                                                                                                  27, 28, 29, 30, 31, 32, 33, 34, 35, 36),"Number_Days"])
which(Hog_island_combined_removed_datesplit_combined_days_2005$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
Hog_island_combined_removed_datesplit_combined_days_2005_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2005[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2006$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
Hog_island_combined_removed_datesplit_combined_days_2006_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2006[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2007$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
Hog_island_combined_removed_datesplit_combined_days_2007_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2007[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2008$Temp <=30) #  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2008_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2008[c( 1,  2,  3,  4 , 5,  6 , 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2009$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2009_dd_30<- sum(Hog_island_combined_removed_datesplit_combined_days_2009[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2010$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2010_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2010[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2011$Temp <=30)  # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2011_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2011[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2012$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2012_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2012[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2013$Temp <=30) #  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2013_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2013[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2014$Temp <=30) #1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2014_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2014[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2015$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
Hog_island_combined_removed_datesplit_combined_days_2015_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2015[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),"Number_Days"])

which(Hog_island_combined_removed_datesplit_combined_days_2016$Temp <=30) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
Hog_island_combined_removed_datesplit_combined_days_2016_dd_30 <- sum(Hog_island_combined_removed_datesplit_combined_days_2016[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),"Number_Days"])

HI_dd30_average <- sum(Hog_island_combined_removed_datesplit_combined_days_1991_dd_30,Hog_island_combined_removed_datesplit_combined_days_1992_dd_30 ,
                       Hog_island_combined_removed_datesplit_combined_days_1993_dd_30, Hog_island_combined_removed_datesplit_combined_days_1994_dd_30,
                       Hog_island_combined_removed_datesplit_combined_days_1995_dd_30, Hog_island_combined_removed_datesplit_combined_days_1996_dd_30 ,
                       Hog_island_combined_removed_datesplit_combined_days_1997_dd_30, Hog_island_combined_removed_datesplit_combined_days_1998_dd_30,
                       Hog_island_combined_removed_datesplit_combined_days_1999_dd_30, Hog_island_combined_removed_datesplit_combined_days_2000_dd_30 ,
                       Hog_island_combined_removed_datesplit_combined_days_2001_dd_30, Hog_island_combined_removed_datesplit_combined_days_2002_dd_30,
                       Hog_island_combined_removed_datesplit_combined_days_2003_dd_30, Hog_island_combined_removed_datesplit_combined_days_2004_dd_30, 
                       Hog_island_combined_removed_datesplit_combined_days_2005_dd_30,  Hog_island_combined_removed_datesplit_combined_days_2006_dd_30,
                       Hog_island_combined_removed_datesplit_combined_days_2007_dd_30, Hog_island_combined_removed_datesplit_combined_days_2008_dd_30,
                       Hog_island_combined_removed_datesplit_combined_days_2009_dd_30 , Hog_island_combined_removed_datesplit_combined_days_2010_dd_30 ,
                       Hog_island_combined_removed_datesplit_combined_days_2011_dd_30 , Hog_island_combined_removed_datesplit_combined_days_2012_dd_30,
                       Hog_island_combined_removed_datesplit_combined_days_2013_dd_30 , Hog_island_combined_removed_datesplit_combined_days_2014_dd_30 ,
                       Hog_island_combined_removed_datesplit_combined_days_2015_dd_30, Hog_island_combined_removed_datesplit_combined_days_2016_dd_30)/26
unique(Hog_island_combined_removed_datesplit_combined_days$Year)  
        
## SM Calculation 
Sherman_Marsh_Sheepscot_removed_combined_days
Sherman_Marsh_Sheepscot_removed_combined_days_summary <- aggregate(TEMP_C ~ Month + Year, Sherman_Marsh_Sheepscot_removed_combined_days, mean)
Sherman_Marsh_Sheepscot_removed_combined_days_summary <- merge(Sherman_Marsh_Sheepscot_removed_combined_days_summary,month_day_number_key_hog_island, by="Month") 

levels(Sherman_Marsh_Sheepscot_removed_combined_days$Year)
colnames(Sherman_Marsh_Sheepscot_removed_combined_days)[11] <- "Location_name"
Sherman_Marsh_Sheepscot_removed_combined_days_1991 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1991")
Sherman_Marsh_Sheepscot_removed_combined_days_1992 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1992")
Sherman_Marsh_Sheepscot_removed_combined_days_1993 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1993")
Sherman_Marsh_Sheepscot_removed_combined_days_1994 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1994")
Sherman_Marsh_Sheepscot_removed_combined_days_1995 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1995")
Sherman_Marsh_Sheepscot_removed_combined_days_1996 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1996")
Sherman_Marsh_Sheepscot_removed_combined_days_1997 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1997")
Sherman_Marsh_Sheepscot_removed_combined_days_1998 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1998")
Sherman_Marsh_Sheepscot_removed_combined_days_1999 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="1999")
Sherman_Marsh_Sheepscot_removed_combined_days_2000 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2000")
Sherman_Marsh_Sheepscot_removed_combined_days_2001 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2001")
Sherman_Marsh_Sheepscot_removed_combined_days_2002 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2002")
Sherman_Marsh_Sheepscot_removed_combined_days_2003 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2003")
Sherman_Marsh_Sheepscot_removed_combined_days_2004 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2004")
Sherman_Marsh_Sheepscot_removed_combined_days_2005 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2005")
Sherman_Marsh_Sheepscot_removed_combined_days_2006 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2006")
Sherman_Marsh_Sheepscot_removed_combined_days_2007 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2007")
Sherman_Marsh_Sheepscot_removed_combined_days_2008 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2008")
Sherman_Marsh_Sheepscot_removed_combined_days_2009 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2009")
Sherman_Marsh_Sheepscot_removed_combined_days_2010 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2010")
Sherman_Marsh_Sheepscot_removed_combined_days_2011 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2011")
Sherman_Marsh_Sheepscot_removed_combined_days_2012 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2012")
Sherman_Marsh_Sheepscot_removed_combined_days_2013 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2013")
Sherman_Marsh_Sheepscot_removed_combined_days_2014 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2014")
Sherman_Marsh_Sheepscot_removed_combined_days_2015 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2015")
Sherman_Marsh_Sheepscot_removed_combined_days_2016 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2016")
Sherman_Marsh_Sheepscot_removed_combined_days_2017 <- Sherman_Marsh_Sheepscot_removed_combined_days_summary %>% filter(Year =="2017")


# SM_dd0
which(Sherman_Marsh_Sheepscot_removed_combined_days_1991$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1991_dd_0 <-0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1992$TEMP_C <=0) #1
Sherman_Marsh_Sheepscot_removed_combined_days_1992_dd_0<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1992[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1993$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1993_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1994$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1994_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1995$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1995_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1996$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1996_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1997$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1997_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1998$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1998_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_1999$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_1999_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2000$TEMP_C <=0) #1
Sherman_Marsh_Sheepscot_removed_combined_days_2000_dd_0<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2000[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2001$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2001_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2002$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2002_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2003$TEMP_C <=0) #2,3
Sherman_Marsh_Sheepscot_removed_combined_days_2003_dd_0<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2003[c(2,3),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2004$TEMP_C <=0) 
Sherman_Marsh_Sheepscot_removed_combined_days_2004_dd_0<- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2005$TEMP_C <=0) # 1
Sherman_Marsh_Sheepscot_removed_combined_days_2005_dd_0<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2005[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2006$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2006_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2007$TEMP_C <=0) # 1 2
Sherman_Marsh_Sheepscot_removed_combined_days_2007_dd_0 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2007[c(1,  2),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2008$TEMP_C <=0) # 
Sherman_Marsh_Sheepscot_removed_combined_days_2008_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2009$TEMP_C <=0) #  
Sherman_Marsh_Sheepscot_removed_combined_days_2009_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2010$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2010_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2011$TEMP_C <=0) #1 
Sherman_Marsh_Sheepscot_removed_combined_days_2011_dd_0 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2011[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2012$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2012_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2013$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2013_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2014$TEMP_C <=0) # 1
Sherman_Marsh_Sheepscot_removed_combined_days_2014_dd_0 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2014[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2015$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2015_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2016$TEMP_C <=0) # 0
Sherman_Marsh_Sheepscot_removed_combined_days_2016_dd_0 <- 0

which(Sherman_Marsh_Sheepscot_removed_combined_days_2017$TEMP_C <=0)
Sherman_Marsh_Sheepscot_removed_combined_days_2017_dd_0 <- 0

SM_dd_0_average <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1991_dd_0,
                       Sherman_Marsh_Sheepscot_removed_combined_days_1992_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_1993_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_1994_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_1995_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_1996_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_1997_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_1998_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_1999_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2000_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_2001_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2002_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_2003_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_2004_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_2005_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_2006_dd_0,
Sherman_Marsh_Sheepscot_removed_combined_days_2007_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2008_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2009_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2010_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2011_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2012_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2013_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2014_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2015_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2016_dd_0, 
Sherman_Marsh_Sheepscot_removed_combined_days_2017_dd_0)/27

# 28 levels
unique(Sherman_Marsh_Sheepscot_removed_combined_days$Year)

# SM_dd_15
which(Sherman_Marsh_Sheepscot_removed_combined_days_1991$TEMP_C <=15) #  1
Sherman_Marsh_Sheepscot_removed_combined_days_1991_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1991[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1992$TEMP_C <=15) #1 2 3
Sherman_Marsh_Sheepscot_removed_combined_days_1992_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1992[c(1,2,3),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1993$TEMP_C <=15) #5
Sherman_Marsh_Sheepscot_removed_combined_days_1993_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1993[c(5),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1994$TEMP_C <=15) #1 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_1994_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1994[c(1, 4,5),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1995$TEMP_C <=15) # 1  2  3  4  5  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_1995_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1995[c(1,  2,  3,  4,  5,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1996$TEMP_C <=15) # 1 2 3 4 5 6 7
Sherman_Marsh_Sheepscot_removed_combined_days_1996_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1996[c(1, 2, 3, 4, 5, 6, 7),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1997$TEMP_C <=15) # 1 2 3 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_1997_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1997[c(1, 2, 3, 4, 5),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1998$TEMP_C <=15) # 1 2 3 4 5 6
Sherman_Marsh_Sheepscot_removed_combined_days_1998_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1998[c(1, 2, 3, 4, 5, 6),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1999$TEMP_C <=15) # 1 2 5 7
Sherman_Marsh_Sheepscot_removed_combined_days_1999_dd_15 <-  sum(Sherman_Marsh_Sheepscot_removed_combined_days_1999[c(1, 2, 5, 7),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2000$TEMP_C <=15) #1  2  3  4  5  6  7 10 11
Sherman_Marsh_Sheepscot_removed_combined_days_2000_dd_15<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2000[c(1,  2,  3,  4,  5,  6,  7, 10, 11),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2001$TEMP_C <=15) #1  2  3  4  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2001_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2001[c(1,  2,  3,  4,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2002$TEMP_C <=15) #1 2 3 5 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2002_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2002[c(1, 2, 3, 5, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2003$TEMP_C <=15) #1 2 3 4 5 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2003_dd_15<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2003[c(1, 2, 3, 4, 5, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2004$TEMP_C <=15) #1  2  3  4  5  6  7  9 10 11
Sherman_Marsh_Sheepscot_removed_combined_days_2004_dd_15<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2004[c(1,  2,  3,  4,  5,  6,  7,  9, 10, 11),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2005$TEMP_C <=15) # 1  2  3  4  5  6  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2005_dd_15<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2005[c(1,  2,  3,  4,  5,  6,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2006$TEMP_C <=15) #1 2 3 4 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2006_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2006[c(1, 2, 3, 4, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2007$TEMP_C <=15) # 1  2  3  4  5  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2007_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2007[c(1,  2,  3,  4,  5,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2008$TEMP_C <=15) # 1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2008_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2008[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2009$TEMP_C <=15) # 1 2 3 6
Sherman_Marsh_Sheepscot_removed_combined_days_2009_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2009[c(1, 2, 3, 6),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2010$TEMP_C <=15) # 1 2 3 4 7 8
Sherman_Marsh_Sheepscot_removed_combined_days_2010_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2010[c(1, 2, 3, 4, 7, 8),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2011$TEMP_C <=15) #1 2 3 4 8 
Sherman_Marsh_Sheepscot_removed_combined_days_2011_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2011[c(1, 2, 3, 4, 8),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2012$TEMP_C <=15) # 1 2 7 9
Sherman_Marsh_Sheepscot_removed_combined_days_2012_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2012[c(1, 2, 7, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2013$TEMP_C <=15) # 1 2 3 4 5 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2013_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2013[c(1, 2, 3, 4, 5, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2014$TEMP_C <=15) # 1 2 3 7
Sherman_Marsh_Sheepscot_removed_combined_days_2014_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2014[c(1, 2, 3, 7),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2015$TEMP_C <=15) # 1 2 3 4 5 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2015_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2015[c(1, 2, 3, 4, 5, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2016$TEMP_C <=15) # 1  2  3  4  5  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2016_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2016[c(1,  2,  3,  4,  5,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2017$TEMP_C <=15) # 1 2 3 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_2017_dd_15 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2017[c(1, 2, 3, 4, 5),"Number_Days"])

SM_dd_15_average <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1991_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_1992_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_1993_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_1994_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1995_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1996_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1997_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1998_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1999_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2000_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2001_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2002_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2003_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2004_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2005_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2006_dd_15,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2007_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2008_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2009_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2010_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2011_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2012_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2013_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2014_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2015_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2016_dd_15, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2017_dd_15)/27

# SM_dd_30
which(Sherman_Marsh_Sheepscot_removed_combined_days_1991$TEMP_C <=30) # 1 
Sherman_Marsh_Sheepscot_removed_combined_days_1991_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1991[c(1),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1992$TEMP_C <=30) #1 2 3
Sherman_Marsh_Sheepscot_removed_combined_days_1992_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1992[c(1,2,3),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1993$TEMP_C <=30) #1 2 3 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_1993_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1993[c(1, 2, 3, 4, 5),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1994$TEMP_C <=30) #1 2 3 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_1994_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1994[c(1, 2, 3, 4, 5),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1995$TEMP_C <=30) # 1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_1995_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1995[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1996$TEMP_C <=30) # 1 2 3 4 5 6 7
Sherman_Marsh_Sheepscot_removed_combined_days_1996_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1996[c(1, 2, 3, 4, 5, 6, 7),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1997$TEMP_C <=30) # 1 2 3 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_1997_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1997[c(1, 2, 3, 4, 5),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1998$TEMP_C <=30) # 1 2 3 4 5 6
Sherman_Marsh_Sheepscot_removed_combined_days_1998_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1998[c(1, 2, 3, 4, 5, 6),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_1999$TEMP_C <=30) # 1 2 3 4 5 6 7
Sherman_Marsh_Sheepscot_removed_combined_days_1999_dd_30 <-  sum(Sherman_Marsh_Sheepscot_removed_combined_days_1999[c(1, 2, 3, 4, 5, 6, 7),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2000$TEMP_C <=30) #1  2  3  4  5  6  7  8  9 10 11
Sherman_Marsh_Sheepscot_removed_combined_days_2000_dd_30<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2000[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2001$TEMP_C <=30) #1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2001_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2001[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2002$TEMP_C <=30) #1 2 3 4 5 6 7 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2002_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2002[c(1, 2, 3, 4, 5, 6, 7, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2003$TEMP_C <=30) #1 2 3 4 5 6 7 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2003_dd_30<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2003[c(1, 2, 3, 4, 5, 6, 7, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2004$TEMP_C <=30) #1  2  3  4  5  6  7  8  9 10 11
Sherman_Marsh_Sheepscot_removed_combined_days_2004_dd_30<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2004[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2005$TEMP_C <=30) # 1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2005_dd_30<- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2005[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2006$TEMP_C <=30) # 1 2 3 4 5 6 7 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2006_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2006[c(1, 2, 3, 4, 5, 6, 7, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2007$TEMP_C <=30) # 1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2007_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2007[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2008$TEMP_C <=30) # 1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2008_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2008[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2009$TEMP_C <=30) # 1 2 3 4 5 6 
Sherman_Marsh_Sheepscot_removed_combined_days_2009_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2009[c(1, 2, 3, 4, 5, 6),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2010$TEMP_C <=30) # 1 2 3 4 5 6 7 8
Sherman_Marsh_Sheepscot_removed_combined_days_2010_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2010[c(1, 2, 3, 4, 5, 6, 7, 8),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2011$TEMP_C <=30) # 1 2 3 4 5 6 7 8 
Sherman_Marsh_Sheepscot_removed_combined_days_2011_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2011[c(1, 2, 3, 4, 5, 6, 7, 8),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2012$TEMP_C <=30) # 1 2 3 4 5 6 7 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2012_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2012[c(1, 2, 3, 4, 5, 6, 7 ,8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2013$TEMP_C <=30) # 1 2 3 4 5 6 7 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2013_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2013[c(1, 2, 3, 4, 5, 6, 7, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2014$TEMP_C <=30) # 1 2 3 4 5 6 7
Sherman_Marsh_Sheepscot_removed_combined_days_2014_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2014[c(1, 2, 3, 4, 5, 6, 7),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2015$TEMP_C <=30) # 1 2 3 4 5 6 7 8 9
Sherman_Marsh_Sheepscot_removed_combined_days_2015_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2015[c(1, 2, 3, 4, 5, 6, 7, 8, 9),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2016$TEMP_C <=30) # 1  2  3  4  5  6  7  8  9 10
Sherman_Marsh_Sheepscot_removed_combined_days_2016_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2016[c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10),"Number_Days"])

which(Sherman_Marsh_Sheepscot_removed_combined_days_2017$TEMP_C <=30) # 1 2 3 4 5
Sherman_Marsh_Sheepscot_removed_combined_days_2017_dd_30 <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_2017[c(1, 2, 3, 4, 5),"Number_Days"])

SM_dd_30_average <- sum(Sherman_Marsh_Sheepscot_removed_combined_days_1991_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_1992_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_1993_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_1994_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1995_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1996_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1997_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1998_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_1999_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2000_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2001_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2002_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2003_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2004_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2005_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2006_dd_30,
                        Sherman_Marsh_Sheepscot_removed_combined_days_2007_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2008_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2009_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2010_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2011_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2012_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2013_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2014_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2015_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2016_dd_30, 
                        Sherman_Marsh_Sheepscot_removed_combined_days_2017_dd_30)/27


# Cape Shore has typically four readings per hour per day
Cape_Shore_combined

# Split cape shore date into separate columns
Cape_Shore_combined_datesplit <- str_split_fixed(Cape_Shore_combined$Date_Time," ", 2)
Cape_Shore_combined_datesplit <- str_split_fixed(Cape_Shore_combined_datesplit[,1],"-", 3)
Cape_Shore_combined_datesplit_combined <- cbind(Cape_Shore_combined, Cape_Shore_combined_datesplit)
colnames(Cape_Shore_combined_datesplit_combined) <- c("Date_Time","temp_C","Year","Month","Day")
head(Cape_Shore_combined_datesplit_combined)

# Calculate Daily average
Cape_Shore_combined_datesplit_combined_daily_average <- aggregate(temp_C ~ Year + Month + Day, Cape_Shore_combined_datesplit_combined, mean)

# Tally observations each year with temp below threshold
unique(Cape_Shore_combined_datesplit_combined_daily_average$Year) # 2013 2014 2016 2017

Cape_Shore_combined_datesplit_combined_daily_average_dd_0_2013 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2013" &
                                                                  temp_C <= "0") %>% summarize(degreedays_0 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_0_2014 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2014" &
                                                          temp_C <= "0") %>% summarize(degreedays_0 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_0_2016 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2016" &
                                                                          temp_C <= "0") %>% summarize(degreedays_0 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_0_2017 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2017" &
                                                                              temp_C <= "0") %>% summarize(degreedays_0 = n())

# dd_15
Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2013 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2013" &
                                                                                                                                    temp_C <= "15") %>% summarize(degreedays_15 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2014 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2014" &
                                                                                                                                    temp_C <= "15") %>% summarize(degreedays_15 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2016 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2016" &
                                                                                                                                    temp_C <= "15") %>% summarize(degreedays_15 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2017 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2017" &
                                                                                                                                    temp_C <= "15") %>% summarize(degreedays_15 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_15_mean <- sum(Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2013,
                                                                       Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2014,
                                                                       Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2016,
                                                                       Cape_Shore_combined_datesplit_combined_daily_average_dd_15_2017)/4
# dd_30
Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2013 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2013" &
                                                                                                                                     temp_C <= "30") %>% summarize(degreedays_30 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2014 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2014" &
                                                                                                                                     temp_C <= "30") %>% summarize(degreedays_30 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2016 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2016" &
                                                                                                                                     temp_C <= "30") %>% summarize(degreedays_30 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2017 <- Cape_Shore_combined_datesplit_combined_daily_average %>% filter(Year=="2017" &
                                                                                                                                     temp_C <= "30") %>% summarize(degreedays_30 = n())
Cape_Shore_combined_datesplit_combined_daily_average_dd_30_mean <- sum(Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2013,
                                                                       Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2014,
                                                                       Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2016,
                                                                       Cape_Shore_combined_datesplit_combined_daily_average_dd_30_2017)/4

# Hope Creek has many readings each day, calculate Daily average
# Split cape shore date into separate columns
HC_Hope_Creek_Combined_daysplit <- str_split_fixed(HC_Hope_Creek_Combined$Date," ", 2)
HC_Hope_Creek_Combined_datesplit <- str_split_fixed(HC_Hope_Creek_Combined_daysplit[,1],"-", 3)
HC_Hope_Creek_Combined_datesplit_combined <- cbind(HC_Hope_Creek_Combined, HC_Hope_Creek_Combined_datesplit)
colnames(HC_Hope_Creek_Combined_datesplit_combined) <- c("Date","TempC","Sal_ppt","Year","Month","Day")
head(HC_Hope_Creek_Combined_datesplit_combined)

HC_Hope_Creek_Combined_datesplit_combined_daily_average <- aggregate(TempC ~ Year + Month + Day, HC_Hope_Creek_Combined_datesplit_combined, mean)
unique(HC_Hope_Creek_Combined_datesplit_combined_daily_average$Year) # 2012 2013 2014 2015

# dd_0
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_0_2012 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2012" &
                                                                                                                                          TempC <= "0") %>% summarize(degreedays_0 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_0_2013 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2013" &
                                                                                                                                          TempC <= "0") %>% summarize(degreedays_0 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_0_2014 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2014" &
                                                                                                                                          TempC <= "0") %>% summarize(degreedays_0 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_0_2015 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2015" &
                                                                                                                                          TempC <= "0") %>% summarize(degreedays_0 = n())
# dd 15
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2012 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2012" &
                                                                                                                                          TempC <= "15") %>% summarize(degreedays_15 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2013 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2013" &
                                                                                                                                          TempC <= "15") %>% summarize(degreedays_15 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2014 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2014" &
                                                                                                                                          TempC <= "15") %>% summarize(degreedays_15 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2015 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2015" &
                                                                                                                                          TempC <= "15") %>% summarize(degreedays_15 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_average <- sum(HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2012,
                                                                             HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2013,
                                                                             HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2014,
                                                                             HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_15_2015)/4
# dd 30
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2012 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2012" &
                                                                                                                                           TempC <= "30") %>% summarize(degreedays_30 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2013 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2013" &
                                                                                                                                           TempC <= "30") %>% summarize(degreedays_30 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2014 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2014" &
                                                                                                                                           TempC <= "30") %>% summarize(degreedays_30 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2015 <- HC_Hope_Creek_Combined_datesplit_combined_daily_average %>% filter(Year=="2015" &
                                                                                                                                           TempC <= "30") %>% summarize(degreedays_30 = n())
HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_average <- sum(HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2012,
                                                                             HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2013,
                                                                             HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2014,
                                                                             HC_Hope_Creek_Combined_datesplit_combined_daily_average_dd_30_2015)/4


# CB Low Salinity and High Salinity have daily readings every fifteen minutes
CB_Low_Salinity_LOLA_CLP_subset_split_combined
unique(CB_Low_Salinity_LOLA_CLP_subset_split_combined$Year) # 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18
CB_Low_Salinity_LOLA_CLP_subset_split_combined$Sal_PSU

# Calculating daily averages 
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average <- aggregate(Temp_C ~ Year + Month + Day, CB_Low_Salinity_LOLA_CLP_subset_split_combined,mean)

# dd_0
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_03 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="03" &
                                                                       Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_04 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="04" &
                                                                        Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_05 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="05" &
                                                                           Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_06 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="06" &
                                                                              Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_07<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="07" &
                                                                      Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_08<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="08" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_09<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="09" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_10<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="10" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_11<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="11" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_12<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="12" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_13<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="13" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_14<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="14" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_15<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="15" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_16<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="16" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_17<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="17" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_0_18<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="18" &
                                                                                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
#dd_15
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_03 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="03" &
                                                                                                                                                  Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_04 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="04" &
                                                                                                                                                  Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_05 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="05" &
                                                                                                                                                  Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_06 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="06" &
                                                                                                                                                  Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_07<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="07" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_08<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="08" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_09<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="09" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_10<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="10" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_11<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="11" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_12<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="12" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_13<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="13" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_14<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="14" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_15<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="15" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_16<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="16" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_17<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="17" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_18<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="18" &
                                                                                                                                                 Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_average <- sum(CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_03,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_04,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_05,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_06,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_07,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_08,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_09,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_10,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_11,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_12,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_13,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_14,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_15,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_16,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_17,
                                                                                  CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_15_18)/16
# dd_30
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_03 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="03" &
                                                                                                                                                   Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_04 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="04" &
                                                                                                                                                   Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_05 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="05" &
                                                                                                                                                   Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_06 <- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="06" &
                                                                                                                                                   Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_07<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="07" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_08<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="08" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_09<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="09" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_10<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="10" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_11<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="11" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_12<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="12" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_13<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="13" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_14<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="14" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_15<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="15" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_16<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="16" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_17<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="17" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_18<- CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average %>% filter(Year=="18" &
                                                                                                                                                  Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30 <- sum(CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_03,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_04,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_05,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_06,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_07,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_08,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_09,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_10,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_11,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_12,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_13,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_14,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_15,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_16,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_17,
                                                                          CB_Low_Salinity_LOLA_CLP_subset_split_combined_daily_average_dd_30_18)/16 
  

# CB High salinity daily readings every fifteen minutes
CB_High_Salinity_DEBY_York_Like_subset

# split the DEBY dates up
CB_High_Salinity_DEBY_York_Like_subset_daysplit <- str_split_fixed(CB_High_Salinity_DEBY_York_Like_subset$DateTimeStamp," ", 2)
CB_High_Salinity_DEBY_York_Like_subset_datesplit <- str_split_fixed(CB_High_Salinity_DEBY_York_Like_subset_daysplit[,1],"/", 3)
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined <- cbind(CB_High_Salinity_DEBY_York_Like_subset, CB_High_Salinity_DEBY_York_Like_subset_datesplit)
colnames(CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined) <- c("Date","Temp_C","Sal_PSU","Month","Day","Year")
head(CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined)

CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average <- aggregate(Temp_C ~ Year + Month + Day, CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined, mean)
unique(CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average$Year) # 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18

# dd_0
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_04 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="04" &
                                                                                 Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_05 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="05" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_06 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="06" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_07 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="07" &
                                                                                    Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_08 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="08" &
                                                                                     Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_09 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="09" &
                                                                                  Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_10 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="10" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_11 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="11" &
                                                                                    Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_12 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="12" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_13 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="13" &
                                                                                    Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_14 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="14" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_15 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="15" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_16 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="16" &
                                                                                   Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_17 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="17" &
                                                                                    Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_18 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="18" & Temp_C <= "0") %>% summarize(degreedays_0 = n())

CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_mean <- sum(CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_04,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_05,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_06,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_07,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_08,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_09,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_10,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_11,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_12,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_13,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_14,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_15,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_16,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_17,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_0_18)/14
# dd_15
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_04 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="04" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_05 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="05" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_06 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="06" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_07 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="07" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_08 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="08" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_09 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="09" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_10 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="10" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_11 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="11" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_12 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="12" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_13 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="13" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_14 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="14" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_15 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="15" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_16 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="16" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_17 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="17" &
                                                                                                                                                                        Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_18 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="18" & Temp_C <= "15") %>% summarize(degreedays_15 = n())

CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_mean <- sum(CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_04,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_05,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_06,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_07,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_08,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_09,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_10,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_11,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_12,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_13,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_14,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_15,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_16,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_17,
                                                                                         CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_15_18)/14
# dd_30
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_04 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="04" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_05 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="05" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_06 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="06" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_07 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="07" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_08 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="08" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_09 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="09" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_10 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="10" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_11 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="11" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_12 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="12" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_13 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="13" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_14 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="14" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_15 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="15" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_16 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="16" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_17 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="17" &
                                                                                                                                                                         Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_18 <- CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average %>% filter(Year=="18" & Temp_C <= "30") %>% summarize(degreedays_30 = n())

CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_mean <- sum(CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_04,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_05,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_06,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_07,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_08,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_09,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_10,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_11,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_12,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_13,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_14,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_15,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_16,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_17,
                                                                                          CB_High_Salinity_DEBY_York_Like_subset_datesplit_combined_daily_average_dd_30_18)/14

# CB high salinity multiple readings per day
head(CB_High_Salinity_HC_VA_subset)

# split the HC-VA dates up
CB_High_Salinity_HC_VA_subset_daysplit <- str_split_fixed(CB_High_Salinity_HC_VA_subset$DateTimeStamp," ", 2)
CB_High_Salinity_HC_VA_subset_datesplit <- str_split_fixed(CB_High_Salinity_HC_VA_subset_daysplit[,1],"/", 3)
CB_High_Salinity_HC_VA_subset_datesplit_combined <- cbind(CB_High_Salinity_HC_VA_subset, CB_High_Salinity_HC_VA_subset_datesplit)
colnames(CB_High_Salinity_HC_VA_subset_datesplit_combined) <- c("Date","Temp_C","Sal_PSU","Month","Day","Year")
head(CB_High_Salinity_HC_VA_subset_datesplit_combined)

CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average <- aggregate(Temp_C ~ Year + Month + Day, CB_High_Salinity_HC_VA_subset_datesplit_combined, mean)
unique(CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average$Year) #03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18

# dd_0
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_03 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="03" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_04 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="04" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_05 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="05" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_06 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="06" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_07 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="07" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_08 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="08" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_09 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="09" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_10 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="10" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_11 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="11" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_12 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="12" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_13 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="13" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_14 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="14" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_15 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="15" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_16 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="16" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_17 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="17" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_18 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="18" & Temp_C <= "0") %>% summarize(degreedays_0 = n())

CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_mean <- sum(CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_03,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_04,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_05,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_06,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_07,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_08,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_09,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_10,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_11,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_12,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_13,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_14,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_15,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_16,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_17,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_0_18)/16

# dd_15
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_03 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="03" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_04 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="04" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_05 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="05" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_06 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="06" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_07 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="07" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_08 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="08" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_09 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="09" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_10 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="10" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_11 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="11" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_12 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="12" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_13 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="13" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_14 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="14" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_15 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="15" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_16 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="16" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_17 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="17" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_18 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="18" & Temp_C <= "15") %>% summarize(degreedays_15 = n())

CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_mean <- sum(CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_03,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_04,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_05,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_06,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_07,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_08,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_09,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_10,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_11,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_12,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_13,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_14,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_15,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_16,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_17,
                                                                                CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_15_18)/16

#dd_30
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_03 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="03" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_04 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="04" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_05 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="05" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_06 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="06" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_07 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="07" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_08 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="08" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_09 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="09" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_10 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="10" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_11 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="11" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_12 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="12" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_13 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="13" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_14 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="14" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_15 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="15" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_16 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="16" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_17 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="17" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_18 <- CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average %>% filter(Year=="18" & Temp_C <= "30") %>% summarize(degreedays_30 = n())

CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_mean <- sum(CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_03,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_04,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_05,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_06,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_07,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_08,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_09,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_10,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_11,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_12,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_13,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_14,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_15,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_16,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_17,
                                                                                 CB_High_Salinity_HC_VA_subset_datesplit_combined_daily_average_dd_30_18)/16


# Laguna madre has multiple readings per day as well
LM_Laguna_Madre_subset

# split up the days 
LM_Laguna_Madre_subset_daysplit <- str_split_fixed(LM_Laguna_Madre_subset$DateTimeStamp," ", 2)
LM_Laguna_Madre_subset_datesplit <- str_split_fixed(LM_Laguna_Madre_subset_daysplit[,1],"/", 3)
LM_Laguna_Madre_subset_datesplit_combined <- cbind(LM_Laguna_Madre_subset, LM_Laguna_Madre_subset_datesplit)
colnames(LM_Laguna_Madre_subset_datesplit_combined) <- c("Date","Temp_C","Sal_PSU","Month","Day","Year")
head(LM_Laguna_Madre_subset_datesplit_combined)

# remove NA's
LM_Laguna_Madre_subset_datesplit_combined_na_omited <- na.omit(LM_Laguna_Madre_subset_datesplit_combined)

# aggregate data 
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average <- aggregate(Temp_C ~ Year + Month + Day, LM_Laguna_Madre_subset_datesplit_combined_na_omited, mean)
unique(LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average$Year) #07 08 09 10 11 12 13 14 15 16 17

# dd_0
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_07 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="07" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_08 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="08" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_09 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="09" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_10 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="10" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_11 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="11" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_12 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="12" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_13 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="13" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_14 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="14" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_15 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="15" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_16 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="16" & Temp_C <= "0") %>% summarize(degreedays_0 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_17 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="17" & Temp_C <= "0") %>% summarize(degreedays_0 = n())

LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_average <- sum(LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_07,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_08,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_09,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_10,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_11,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_12,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_13,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_14,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_15,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_16,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_0_17)/11
# dd_15
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_07 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="07" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_08 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="08" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_09 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="09" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_10 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="10" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_11 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="11" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_12 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="12" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_13 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="13" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_14 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="14" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_15 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="15" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_16 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="16" & Temp_C <= "15") %>% summarize(degreedays_15 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_17 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="17" & Temp_C <= "15") %>% summarize(degreedays_15 = n())

LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_average <- sum(LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_07,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_08,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_09,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_10,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_11,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_12,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_13,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_14,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_15,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_16,
                                                                                      LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_15_17)/11
# dd_30
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_07 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="07" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_08 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="08" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_09 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="09" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_10 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="10" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_11 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="11" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_12 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="12" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_13 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="13" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_14 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="14" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_15 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="15" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_16 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="16" & Temp_C <= "30") %>% summarize(degreedays_30 = n())
LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_17 <- LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average %>% filter(Year=="17" & Temp_C <= "30") %>% summarize(degreedays_30 = n())

LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_average <- sum(LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_07,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_08,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_09,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_10,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_11,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_12,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_13,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_14,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_15,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_16,
                                                                                       LM_Laguna_Madre_subset_datesplit_combined_na_omited_daily_average_dd_30_17)/11


##### PLOTTING #####

# Plotting longitude and latitude of samples
#install.packages(c("maps","mapdata"))
library(maps)
library(mapdata)
library(ggplot2)
library(Hmisc)
#install.packages('ggrepel')
library(ggrepel)
library(tidyr)
library(reshape2)
library(ggpubr)

# helpful tutorial
#https://sarahleejane.github.io/learning/r/2014/09/21/plotting-data-points-on-maps-with-r.html

# Points to plot for each population and the population name
pop_coordinates <- read.csv("./Environmental_Data/population_coordinates.csv")
pop_coordinates$Lat <- as.character(pop_coordinates$Lat)
pop_coordinates$Lat <- as.numeric(pop_coordinates$Lat)
pop_coordinates$Long <- as.character(pop_coordinates$Long)
pop_coordinates$Long <- as.numeric(pop_coordinates$Long)
pop_coordinates <- na.omit(pop_coordinates)

#plot east coast of the US
states <- map_data("state")
east_coast <- subset(states, region %in% c("maine","connecticut","new hampshire","new york",
                                          "delaware","virginia","georgia","florida","north carolina",
                                          "south carolina","texas","louisiana","alabama","rhode island","maryland","massachusetts",
                                         "pennsylvania","mississippi"))
east_coast_map <- ggplot(data=east_coast) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "light blue", color = "black") + 
  coord_fixed(1.3) + ggtitle("Oyster Genome Resequencing Sample Populations") 
#Add coordinates for samples   
east_coast_map + geom_point(data=pop_coordinates, aes(x=Long,y=Lat, color="red"))  
  
# Add labels and Repel ID's to make it a finished map
#Finished Map
finished_map <- east_coast_map + 
  geom_point(data=pop_coordinates, aes(x=Long,y=Lat, color="red"))  +
geom_label_repel(data=pop_coordinates, 
          aes(x=Long, y=Lat, label=Pop.ID.), box.padding = 0.2, point.padding = 0.2, 
          segment.colour = 'grey40') + xlab("Longitude") + ylab("Latitude") 

## Graph for Temperature, Salinity, dd_0, dd_15, dd_30 

# Load in edited spreadsheet to retrieve this information

combined_spreadsheet <- read.csv("./Environmental_Data/sequenced_sample_info_Jan.23_2019_V1.csv", header=TRUE)
head(combined_spreadsheet)
combined_spreadsheet_subset <- combined_spreadsheet[,c(3,5,6,7,10,11,20:26)]
head(combined_spreadsheet_subset)
combined_spreadsheet_subset<- unique(combined_spreadsheet_subset)
combined_spreadsheet_subset_maxtemp_sorted <- combined_spreadsheet_subset[order(combined_spreadsheet_subset$Max_temperature_Celsius),]

combined_spreadsheet_subset_maxtemp_sorted$Pop.ID. <- factor(combined_spreadsheet_subset_maxtemp_sorted$Pop.ID., levels = combined_spreadsheet_subset_maxtemp_sorted$Pop.ID.)

# reshape data with tidyr
combined_spreadsheet_subset_maxtemp_sorted_long <- melt(combined_spreadsheet_subset_maxtemp_sorted, id="Pop.ID.")

# Make temperature data frame
combined_spreadsheet_subset_maxtemp_sorted_long_temp <- combined_spreadsheet_subset_maxtemp_sorted_long[c(81:128),]
combined_spreadsheet_subset_maxtemp_sorted_long_temp <- na.omit(combined_spreadsheet_subset_maxtemp_sorted_long_temp) # remove empty UMFS
combined_spreadsheet_subset_maxtemp_sorted_long_temp$value <- as.numeric(combined_spreadsheet_subset_maxtemp_sorted_long_temp$value)


# Plotting Temperature
temp_plot <- ggplot(combined_spreadsheet_subset_maxtemp_sorted_long_temp, aes(x=Pop.ID., y=value, fill=variable)) + geom_col(position="dodge") +
  xlab("Population ID") + ylab("Temperature (C)") + ggtitle("Maximum, Mininum, Annual Temperature Across Populations") + coord_flip() + 
  scale_x_discrete(limits=c("LM","SL","CL", "OBOYS2","DEBY","HC-VA","LOLA","CLP","CS","HG","NEH","NG","HC","HI","SM")) 
temp_plot_finished <- temp_plot + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Temperature", labels=c("Mean Annual Temperature",
                                                      "Maximum Temperature","Minimum Temperature"))
# Salinity data frame 
combined_spreadsheet_subset_maxtemp_sorted_long_salinity <- combined_spreadsheet_subset_maxtemp_sorted_long[c(129:144),]
combined_spreadsheet_subset_maxtemp_sorted_long_salinity <- na.omit(combined_spreadsheet_subset_maxtemp_sorted_long_salinity)
combined_spreadsheet_subset_maxtemp_sorted_long_salinity$value <- as.numeric(combined_spreadsheet_subset_maxtemp_sorted_long_salinity$value)

# Plotting salinity
salinity_plot <- ggplot(combined_spreadsheet_subset_maxtemp_sorted_long_salinity, aes(x=Pop.ID., y=value, fill=variable)) + geom_col(position="dodge") +
  xlab("Population ID") + ylab("Salinity (ppt)") + ggtitle("Mean Annual Salinity Across Populations") + coord_flip() + 
  scale_x_discrete(limits=c("LM","SL","CL", "OBOYS2","DEBY","HC-VA","LOLA","CLP","CS","HG","NEH","NG","HC","HI","SM")) 

salinity_plot_finished <- salinity_plot + scale_fill_manual(values="plum4", name="Salinity", labels="Mean Annual Salinity")

# Disease pressure data frame
combined_spreadsheet_subset_maxtemp_sorted_long_disease_pressure <- combined_spreadsheet_subset_maxtemp_sorted_long[c(145:192),]
combined_spreadsheet_subset_maxtemp_sorted_long_disease_pressure <- na.omit(combined_spreadsheet_subset_maxtemp_sorted_long_disease_pressure) #omit empty UMFS
combined_spreadsheet_subset_maxtemp_sorted_long_disease_pressure$value <- as.numeric(combined_spreadsheet_subset_maxtemp_sorted_long_disease_pressure$value)

# Plotting disease pressure
disease_plot <- ggplot(combined_spreadsheet_subset_maxtemp_sorted_long_disease_pressure, aes(x=Pop.ID., y=value, fill=variable)) + geom_col(position="dodge") +
  xlab("Population ID") + ylab("Number of Days") + ggtitle("Mean Degrees Days Below Temperature Across Populations") + coord_flip() + 
  scale_x_discrete(limits=c("LM","SL","CL", "OBOYS2","DEBY","HC-VA","LOLA","CLP","CS","HG","NEH","NG","HC","HI","SM")) 


disease_plot_finished <- disease_plot + scale_fill_manual(values=c("#cb6a49","#a46cb7","#7aa457"), name="Number of Degree Days", labels=c("Degree Days Below 0 C",
                                                "Degree Days Below 15 C", "Degree Days Below 30 C"))

# combine plots side by side
# using library grid extra
library(gridExtra)
grid.arrange(finished_map, 
          arrangeGrob(temp_plot_finished,  disease_plot_finished,salinity_plot_finished, ncol=3), nrow=2)




