#! /bin/bash  

## Set variables

datapath=/path/to/short/reads
# file extension eg bam, .sorted.rm.bam, fastq.gz
filetype=bam
output=SPECIES_depth.tsv




## Code   

first=true
for file in "$datapath"/*."$filetype"; do 
if $first; then
        echo -e "Sample_id\tAverage_depth"
        first=false
fi
filebase=$(basename "$file" ".$filetype")
filetag=${filebase##*/}
depth=$(samtools depth "$file"  | awk '{sum+=$3} END { print sum/NR}')
echo -e "${filetag}\t${depth}"
done > $output
