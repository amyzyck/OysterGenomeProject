#!/usr/bin/env bash

cat namelist | parallel -j 8 "java -Xms4g -jar /usr/local/bin/picard.jar MarkDuplicates I={}-RG.bam O={}-RGmd.bam M={}_dup_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> md.{}.log"

echo -e "Picard has finished  in" `pwd` | mailx -s "Analysis has finished" jpuritz@uri.edu

