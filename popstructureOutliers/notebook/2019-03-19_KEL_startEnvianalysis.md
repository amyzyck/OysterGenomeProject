# Start Environment analsysis with Kevin

### Step 1 organize data
`data/modified_samplemetadata.csv` same as "original" in spreadsheet

`data/PopPlotting.csv` info for colors and pop IDs for plotting

`data/environment/spreadsheet.xlsx` sequenced_sample tab has environment information for each population. We to convert to CSV and check if the sample_name IDs MATCH the sample name IDs in `data/modified_samplemetadata.csv`. Also output coordinates.csv tab into own csv file.

### Step 2 learn lfmm_ridge

https://github.com/TestTheTests/TTT_RecombinationGenomeScans/blob/master/src/b_Proc_Sims.R

Starting on line 519

Spearmans correlation - line 660

### Step 3 run association test

- load the relevant vcf file
- subset to relevant wild populations minus LM (Bodie can walk you through that)
- `data/modified_samplemetadata.csv` see column `Wild.Sel` and "W" for wild
- filter SNPs again for sites that are all the same genotype across populations (Bodie already doing)
- convert vcf to the "G matrix" (Bodie already doing)
- run 
    - lfmm_ridge, Y is genotype, X would be our environmental variable like temp
      - look into extracting PCs from the ridge (Alan)
    - raw Spearman's correlation between allele frequency and environment
    - do for all environment variables and obtain slopes and p-value for each SNP
    
- check to see how PCs from the population structure analysis correlate with the environments
 
 



