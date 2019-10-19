#!/bin/bash
# Download data from Jon's KITT server
wget http://kitt.uri.edu/sorted.ref3.0*
mkdir GenomicRegionFiles
mv -r sorted.ref3.0* GenomicRegionFiles/
# Finding which outliers intersect with each region identified
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.CDS.bed > Outliers_CDS.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.exon.bed > Outliers_exon.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.UTR.bed > Outliers_UTR.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.intergenic.bed > Outliers_Intergenic.txt
bedtools intersect -u -a FstOutliers_Unrelated.bed -b GenomicRegionFiles/sorted.ref3.0.gene.bed > Outliers_gene.txt
