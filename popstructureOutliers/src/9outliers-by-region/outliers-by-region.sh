#!/bin/bash
# Download data from Jon's KITT server
wget http://kitt.uri.edu/sorted.ref3.0*
# Finding which outliers intersect with each region identified
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.CDS.bed > Outliers_CDS.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.CDS.sc.bed > Outliers_CDS_sc.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.exon.bed > Outliers_exon.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.exon.sc.bed > Outliers_exon_sc.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.UTR.bed > Outliers_UTR.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.intergenic.bed > Outliers_Intergenic.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.gene.bed > Outliers_gene.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.gene.sc.bed > Outliers_gene_sc.txt
