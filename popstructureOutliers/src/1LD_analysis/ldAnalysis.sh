#!bin/bash
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out geno_ld_window_50-50_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 200 --ld-window-bp 250 --out geno_ld_window_200-250_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 450 --ld-window-bp 500 --out geno_ld_window_450-500_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 4500 --ld-window-bp 5000 --out geno_ld_window_4500-5000_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 9500 --ld-window-bp 10000 --out geno_ld_window_9500-10000_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 49500 --ld-window-bp 50000 --out geno_ld_window_49500-50000_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 99500 --ld-window-bp 100000 --out geno_ld_window_99500-100000_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --keep unrelated.txt --geno-r2 --ld-window-bp-min 495000 --ld-window-bp 500000 --out geno_ld_window_495000-500000_unrelated
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out geno_ld_window_50-50_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 200 --ld-window-bp 250 --out geno_ld_window_200-250_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 450 --ld-window-bp 500 --out geno_ld_window_450-500_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 4500 --ld-window-bp 5000 --out geno_ld_window_4500-5000_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 9500 --ld-window-bp 10000 --out geno_ld_window_9500-10000_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 49500 --ld-window-bp 50000 --out geno_ld_window_49500-50000_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 99500 --ld-window-bp 100000 --out geno_ld_window_99500-100000_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove exclude_LM.txt --geno-r2 --ld-window-bp-min 495000 --ld-window-bp 500000 --out geno_ld_window_495000-500000_excludingLM
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out geno_ld_window_50-50_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 200 --ld-window-bp 250 --out geno_ld_window_200-250_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 450 --ld-window-bp 500 --out geno_ld_window_450-500_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 4500 --ld-window-bp 5000 --out geno_ld_window_4500-5000_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 9500 --ld-window-bp 10000 --out geno_ld_window_9500-10000_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 49500 --ld-window-bp 50000 --out geno_ld_window_49500-50000_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 99500 --ld-window-bp 100000 --out geno_ld_window_99500-100000_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeSelection.txt --geno-r2 --ld-window-bp-min 495000 --ld-window-bp 500000 --out geno_ld_window_495000-500000_excludingSelection
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out geno_ld_window_50-50_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 200 --ld-window-bp 250 --out geno_ld_window_200-250_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 450 --ld-window-bp 500 --out geno_ld_window_450-500_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 4500 --ld-window-bp 5000 --out geno_ld_window_4500-5000_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 9500 --ld-window-bp 10000 --out geno_ld_window_9500-10000_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 49500 --ld-window-bp 50000 --out geno_ld_window_49500-50000_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 99500 --ld-window-bp 100000 --out geno_ld_window_99500-100000_excludingWild
vcftools --gzvcf Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz --remove excludeWild.txt --geno-r2 --ld-window-bp-min 495000 --ld-window-bp 500000 --out geno_ld_window_495000-500000_excludingWild