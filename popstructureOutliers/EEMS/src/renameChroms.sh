#!/usr/bin/bash

##############################################################################
#
# File   : renameChroms.bash
# History: 10/16/2018 - Created by Kevin Freeman (KF)
#
##############################################################################
#
# This script uses sed to swap out refseq chromosome names with the corresponding
# chromosome number for the Eastern Oyster. It takes a vcf file as input, but 
# should be able to be used on any other kind of text file. It uses a hash 
# table so it requires bash 4
#
##############################################################################
file=$1

declare -A chroms
chroms=( ["1"]="NC_035780.1" ["2"]="NC_035781.1" ["3"]="NC_035782.1" \
	["4"]="NC_035783.1" ["5"]="NC_035784.1" ["6"]="NC_035785.1" \
	["7"]="NC_035786.1" ["8"]="NC_035787.1" ["9"]="NC_035788.1" \
	["10"]="NC_035789.1" )

for i in "${!chroms[@]}"
do
	echo $i
	sed -i "s/${chroms[$i]}/$i/g" $file
done
