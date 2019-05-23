# This will serve as a central repository for all IGV data tracks for the Oyster Genome Project.

# XML Files

01-base.xml

`https://raw.githubusercontent.com/hputnam/FROGER/master/igv-xml/01-base.xml`

02-cafe.xml

`https://raw.githubusercontent.com/hputnam/FROGER/master/igv-xml/02-cafe.xml`

03-coverage.xml (coverage for each population)

`https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/IGV%20Tracks/coverage.xml`

04-IGVTracks_OutlierDetectionEnviAssocCombined.xml

`https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/IGV%20Tracks/IGVTracks_OutlierDetectionEnviAssocCombined.xml`

05 - Base_CAFE_dup_disease_loci_missassembly.xml
`https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/IGV%20Tracks/Base_CAFE_dup_disease_loci_missassembly.xml`


# Population Level Specific Files that can be accessed via `http://kitt.uri.edu/filename`:

### VCF Files for each population (with allele frequencies)

#### Locality labels
```
CL
CLP
CS
DEBY
HC
HG
HI
LM
LOLA
NEH
NG
OBOYS2
SL
SM
UMFS
```

#### Type of Files
|Extension| Data |
|---------|-------|
|`*.Combined.SNP.TRSdp5g1FnDNAmaf052alleles.empty.vcf.gz` |All biallelic SNPs called in all populations with a minor allele frequency greater than 0.05|
|`*.Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData50000Window_exclude_LM.empty.vcf.gz`| Thinned SNPs using 50,000 bp Window|
|`*.Thinned.SNP.TRSdp5g1FnDNAmaf052alleles.thinnedMatrixAndMetaData5000Window_exclude_LM.empty.vcf.gz.`| Thinned SNPs using 5,000 bp Window|
