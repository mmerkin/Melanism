#!/bin/bash

datapath=/pub64/mattm/short_reads/Ob_reads/nigra

for file in "$datapath"/*; do 
filetag=${file##*/}
echo $filetag
mkdir Ob_output_$filetag
/pub64/mattm/apps/pangenie/build/src/PanGenie -f Ob_pangenie -i <(zcat $file/*R1*.fastq.gz $file/*R2*.fastq.gz $file/*R0*.fastq.gz ) -s Ob_output_$filetag -j 16 -t 30 -o Ob_output_$filetag
done

