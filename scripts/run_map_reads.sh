#! /bin/bash

# Requirements: bwa-mem2

# Reference genome
REF=/path/to/bwa/indexed/reference.fa
# Short reads, each read pair should be within a directory named after the sample eg datapath/sample/sample_R1.fastq.gz
datapath=/path/to/reads/directory
# file to place final bams
output=/path/to/output/file
# Number of threads to use
threads=32
# Change to false to keep temporary files
remove_temp=true





# Code

mkdir -p $output

for file in "$datapath"/*; do 
filetag=${file##*/}

echo "Mapping reads"
bwa-mem2 mem -t $threads $REF $file/*R1*.fastq.gz $file/*R2*.fastq.gz > "$output/${filetag}.raw.bam"

echo "Filtering reads"

samtools view -@ $threads -b -F 2828 -q 20 "$output/$filetag/$filetag.raw.bam" -o "$output/$filetag/$filetag.filtered.bam"

echo "Sorting by name"
samtools sort -@ $threads -n "$output/${filetag}.filtered.bam" -o "$output/${filetag}.sorted.n.bam"

# Fixmate -r filters the reads and -m adds a mate score tag for markdup to select the best reads to keep

echo "Adding fixmate tags"
samtools fixmate -rm -@ $threads "$output/${filetag}.sorted.n.bam" "$output/${filetag}.fixmate.bam"

echo "Sorting by position"
samtools sort -@ $threads "$output/${filetag}.fixmate.bam" -o "$output/${filetag}.sorted.p.bam"

echo "Marking duplicates"
samtools markdup -r -@ $threads "$output/${filetag}.sorted.p.bam" "$output/${filetag}.bam"
echo "running index"
samtools index "$output/${filetag}.bam"

if $remove_temp; then
  rm "$output/${filetag}.raw.bam" "$output/${filetag}.sorted.n.bam" "$output/${filetag}.fixmate.bam" "$output/${filetag}.sorted.p.bam"
fi

done

