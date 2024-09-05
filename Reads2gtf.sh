#!/bin/bash


#Required tools: minimap2 bedtools ucsc_tools
#Install with conda: conda create -n Reads2gtf -c bioconda minimap2 bedtools ucsc-bedtogenepred ucsc-genepredtogtf 

read -p "What is the reference called? " REF 
read -p "What is the read fasta file called? " READS
BASE="${READS%.*}"

minimap2 -ax splice --cs $REF $READS | samtools sort -O BAM - > ${BASE}.bam

bedtools bamtobed -bed12 -i ${BASE}.bam > ${BASE}.bed

bedToGenePred ${BASE}.bed ${BASE}.genepred

genePredToGtf "file" ${BASE}.genepred ${BASE}.gtf

rm ${BASE}.bam
rm ${BASE}.bed
rm ${BASE}.genepred
