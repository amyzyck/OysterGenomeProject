library(ape)
library(tidyr)

#set working directory to a location that has this gff file and the RNAseq results files
setwd()
#Load in the gff file
GFF_table_Cvirginica <- read.gff("GCF_002022765.2_C_virginica-3.0_genomic.gff")
#subset and keep mRNA sequences
GFF_table_Cvirginica_mRNA <- subset(GFF_table_Cvirginica, type == "mRNA" ) 
head(GFF_table_Cvirginica_mRNA)
dim(GFF_table_Cvirginica_mRNA)
#[1] 60201     9
#There are 60,201 mRNA transcripts in the genome
#let's break up the attributes column in order to get the Parent ID's we want to use for this analysis
att <- tidyr::separate(data.frame(text = GFF_table_Cvirginica_mRNA$attributes), text, into = c("ID","Parent","Dbxref","Genbank","Name","gbkey","gene","model_evidence","transcript_id"), sep = ";", fill = "warn", extra = "warn")
att <- att[,1:6]
att[,1] <- sub("^ID=","",att[,1])
att[,2] <- sub("^Parent=","",att[,2])
att[,4] <- sub("^Name=","",att[,4])
att[,5] <- sub("^gbkey=","",att[,5])
att[,6] <- sub("^gene=","",att[,6])
colnames(att)[4] <- c("XM")
head(att)
dim(att)
#[1] 60201     6

#We also want this information for each gene
GFF_table_Cvirginica_gene <- subset(GFF_table_Cvirginica, type == "gene" ) 
head(GFF_table_Cvirginica_gene)
att2 <- tidyr::separate(data.frame(text = GFF_table_Cvirginica_gene$attributes), text, into = c("ID","Dbxref","Name","gbkey","gene_biotype"), sep = ";", fill = "warn", extra = "warn")
dim(att2)
#[1] 38838     5
#There are 38,838 genes in the genome
att2[,1] <- sub("^ID=","",att2[,1])
att2[,2] <- sub("^Dbxref=","",att2[,2])
att2[,3] <- sub("^Name=","",att2[,3])
att2[,4] <- sub("^gbkey=","",att2[,4])
att2[,5] <- sub("^gene=","",att2[,5])
colnames(att2)[1] <- c("Parent")
#let's cbind the att2 with GFF_table_Cvirginica_gene. We have kept all rows in the same order so there will be no mixup.
#This can be confirmed by comparing the Parent ID column with the attribut colun in the final GFF_table_Cvirginica_gene2 table.
GFF_table_Cvirginica_gene2 <- cbind(GFF_table_Cvirginica_gene,att2)
head(GFF_table_Cvirginica_gene2)
dim(GFF_table_Cvirginica_gene2)
#[1] 38838     14


#OK, now we are ready to read the files in.
resistant <- read.csv("Proestou_resistant_day7_with_position.csv")
head(resistant)
resistant_w_gene <- merge(resistant,att, by = "XM")
head(resistant_w_gene)
dim(resistant_w_gene)
#[1] 28726    15
#OK we have 28,726 mRNA sequences, lets sort by both Parent ID and increasing FDR value.  
resistant_w_gene2 <- resistant_w_gene[order(resistant_w_gene$Parent, resistant_w_gene$FDR ), ]
#Now we can remove all duplicated Parent names keeping the first value that occures, which after sorting by FDR will be the most significant Parent ID.
resistant_w_gene3 <- resistant_w_gene2[ !duplicated(resistant_w_gene2$Parent), ]   
dim(resistant_w_gene3)
# 20246    15
head(resistant_w_gene3)
#So we ahve 20,246 genes with scores.
#Let's make a smaller file to share, but we will keep the XM and XP sequences with these for now as well.
resistant_gene <- resistant_w_gene3[,c(12,1,3,4,5,6,7)]
head(resistant_w_gene)
#Now let's add start and stop for each gene
resistant_gene_final <- merge(resistant_gene,GFF_table_Cvirginica_gene2,by="Parent")
dim(resistant_gene_final)
#[1] 20246    20
head(resistant_gene_final)
write.csv(resistant_gene_final,"Proestou_resistant_day7_gene.csv",row.names = F,quote=F)

susceptible <- read.csv("Proestou_susceptible_day7_with_position.csv")
dim(susceptible)
#[1] 30499    10
susceptible_w_gene <- merge(susceptible,att, by = "XM")
head(susceptible_w_gene)
dim(susceptible_w_gene)
#[1] 30499    15
#OK we have 30,499 mRNA sequences, lets sort by both Parent ID and increasing FDR value.  
susceptible_w_gene2 <- susceptible_w_gene[order(susceptible_w_gene$Parent, susceptible_w_gene$FDR ), ]
#Now we can remove all duplicated Parent names keeping the first value that occures, which after sorting by FDR will be the most significant Parent ID.
susceptible_w_gene3 <- susceptible_w_gene2[ !duplicated(susceptible_w_gene2$Parent), ]   
dim(susceptible_w_gene3)
# 21216    15
head(susceptible_w_gene3)
#So we ahve 21,216 genes with scores.
#Let's make a smaller file to share, but we will keep the XM and XP sequences with these for now as well.
susceptible_gene <- susceptible_w_gene3[,c(12,1,3,4,5,6,7)]
head(susceptible_gene)
#Now let's add start and stop for each gene
susceptible_gene_final <- merge(susceptible_gene,GFF_table_Cvirginica_gene2,by="Parent")
dim(susceptible_gene_final)
#[1] 21216    20
head(susceptible_gene_final)
write.csv(susceptible_gene_final,"Proestou_susceptible_day7_gene.csv",row.names = F,quote=F)

tolerant <- read.csv("Proestou_tolerant_day7_with_position.csv")
dim(tolerant)
tolerant_w_gene <- merge(tolerant,att, by = "XM")
head(tolerant_w_gene)
dim(tolerant_w_gene)
#[1] 32114    15
#OK we have 32,114 mRNA sequences, lets sort by both Parent ID and increasing FDR value.  
tolerant_w_gene2 <- tolerant_w_gene[order(tolerant_w_gene$Parent, tolerant_w_gene$FDR ), ]
#Now we can remove all duplicated Parent names keeping the first value that occures, which after sorting by FDR will be the most significant Parent ID.
tolerant_w_gene3 <- tolerant_w_gene2[ !duplicated(tolerant_w_gene2$Parent), ]   
dim(tolerant_w_gene3)
# 22150    15
head(tolerant_w_gene3)
#So we ahve 22,150 genes with scores.
#Let's make a smaller file to share, but we will keep the XM and XP sequences with these for now as well.
tolerant_gene <- tolerant_w_gene3[,c(12,1,3,4,5,6,7)]
head(tolerant_w_gene)
#Now let's add start and stop for each gene
tolerant_gene_final <- merge(tolerant_gene,GFF_table_Cvirginica_gene2,by="Parent")
dim(tolerant_gene_final)
#[1] 22150    20
head(tolerant_gene_final)
write.csv(tolerant_gene_final,"Proestou_tolerant_day7_gene.csv",row.names = F,quote=F)


salinity <- read.csv("Salinity.Diff.Expression.csv")
colnames(salinity)[1] <- c("XM")
salinity_w_gene <- merge(salinity,att, by = "XM")
head(salinity_w_gene)
dim(salinity_w_gene)
#[1] 42122    18
#OK we have 28,726 mRNA sequences, lets sort by both Parent ID and increasing FDR value.  
salinity_w_gene2 <- salinity_w_gene[order(salinity_w_gene$Parent, salinity_w_gene$FDR ), ]
#Now we can remove all duplicated Parent names keeping the first value that occures, which after sorting by FDR will be the most significant Parent ID.
salinity_w_gene3 <- salinity_w_gene2[ !duplicated(salinity_w_gene2$Parent), ]   
dim(salinity_w_gene3)
# 26163    18
head(salinity_w_gene3)
#So we ahve 26,613 genes with scores.
#Let's make a smaller file to share, but we will keep the XM and XP sequences with these for now as well.
salinity_w_gene <- salinity_w_gene3[,c(15,1,2,3,4,5,9)]
head(salinity_w_gene)
salinity_w_gene$Significant <- ifelse(salinity_w_gene$FDR <= 0.05,"True","False")
table(salinity_w_gene$Significant)
#Now let's add start and stop for each gene
salinity_gene_final <- merge(salinity_w_gene,GFF_table_Cvirginica_gene2,by="Parent")
dim(salinity_gene_final)
#[1] 26163    20
head(salinity_gene_final)

write.csv(salinity_gene_final,"Salinity.Diff.Expression_gene.csv",row.names = F,quote=F)


temperature <- read.csv("Temperature.Diff.Expression.csv")
colnames(temperature)[1] <- c("XM")
temperature_w_gene <- merge(temperature,att, by = "XM")
head(temperature_w_gene)
dim(temperature_w_gene)
#[1] 44314    12
#OK we have 44,314 mRNA sequences, lets sort by both Parent ID and increasing FDR value.  
temperature_w_gene2 <- temperature_w_gene[order(temperature_w_gene$Parent, temperature_w_gene$FDR ), ]
#Now we can remove all duplicated Parent names keeping the first value that occures, which after sorting by FDR will be the most significant Parent ID.
temperature_w_gene3 <- temperature_w_gene2[ !duplicated(temperature_w_gene2$Parent), ]   
dim(temperature_w_gene3)
# 27687    12
head(temperature_w_gene3)
#So we ahve 27,687 genes with scores.
#Let's make a smaller file to share, but we will keep the XM and XP sequences with these for now as well.
temperature_w_gene <- temperature_w_gene3[,c(9,1,2,3,4,5)]
head(temperature_w_gene)
temperature_w_gene$Significant <- ifelse(temperature_w_gene$FDR <= 0.05,"True","False")
table(temperature_w_gene$Significant)
#Now let's add start and stop for each gene
temperature_gene_final <- merge(temperature_w_gene,GFF_table_Cvirginica_gene2,by="Parent")
dim(temperature_gene_final)
#[1] 27687    20
head(temperature_gene_final)
write.csv(temperature_gene_final,"temperature.Diff.Expression_gene.csv",row.names = F,quote=F)

