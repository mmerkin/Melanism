#! /bin/bash  

## Set variables

datapath=/pub64/mattm/read_mapping/Ob/dedup
# file extension eg bam, fastq.gz
filetype=dedup.bam
output=Ob_depth.tsv




## Code   

first=true
for file in "$datapath"/*."$filetype"; do 
if $first; then
        echo -e "Sample_id\tAverage_depth"
else
        filebase=$(basename "$file" ".$filetype")
        filetag=${filebase##*/}
        depth=$(samtools depth "$file"  | awk '{sum+=$3} END { print sum/NR}')
        echo -e "${filetag}\t${depth}"
fi
first=false   
done > $output
