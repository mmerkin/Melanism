#!/bin/bash

# Initialize an empty string
concatenated_string=""
datapath=/pub64/mattm/read_mapping/Ob/dedup
REF=Ob_genome.fa
OUTPUT=Ob_SVs

# Loop through a set of elements
for file in "$datapath"/*; do 
        if [[ "$file" =~ \.bai$ ]]; then
                continue
        fi
        bam="--bam ${file} "
        # Concatenate the current element to the string
        concatenated_string="${concatenated_string}${bam}"
done

# Print the final concatenated string
echo $concatenated_string


# Start manta

configManta.py $concatenated_string  --referenceFasta $REF --runDir $OUTPUT
