#!/bin/bash

REF=Pp_scaffolds.fa
THREADS=32

datapath=/pub64/mattm/read_mapping/Pp/dedup


for file in "$datapath"/*; do 
filebase=$(basename "$file" .dedup.bam)
if [[ "$file" =~ \.bai$ ]]; then
        continue
fi
if [[ "$filebase" =~ m$ ]]; then
        continue
fi

mkdir $filebase

insurveyor.py --threads $THREADS $file $filebase $REF

done

echo "completed"
