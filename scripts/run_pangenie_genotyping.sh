#! /bin/bash

datapath=/path/to/short/reads
pangenie=/path/to/pangenie
outfile=outfile_name # Must be the same as the indexing step
kmer_threads=16
genotype_threads=32

for file in "$datapath"/*; do 
filetag=${file##*/}
echo $filetag
$pangenie/PanGenie -f $outfile -i <(zcat $file/*R1*.fastq.gz $file/*R2*.fastq.gz $file/*R0*.fastq.gz ) -s $filetag -j $kmer_threads -t $genotype_threads
done
