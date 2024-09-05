#!/bin/bash
threads=24

datapath=/pub64/mattm/read_mapping/Ob/bams/nigra

for file in "$datapath"/*; do 
filebase=$(basename "$file" .sorted.bam)
echo "Moving to sample: ${filebase}"
echo "running filter"
samtools view -b -F 0xc "${datapath}/${filebase}.sorted.bam"  -o "filter/${filebase}.filtered.bam"
echo "running nsort"
samtools sort -@ $threads -n "filter/${filebase}.filtered.bam" -o "nsort/${filebase}.sorted.n.bam"
echo "running fixmate"
samtools fixmate -m "nsort/${filebase}.sorted.n.bam" "fixmate/${filebase}.fixmate.bam"
echo "running psort"
samtools sort -@ $threads "fixmate/${filebase}.fixmate.bam" -o "psort/${filebase}.sorted.p.bam"
echo "running dedup"
samtools markdup -r -@ $threads "psort/${filebase}.sorted.p.bam" "dedup/${filebase}.dedup.bam"
echo "running index"
samtools index "dedup/${filebase}.dedup.bam"
done
