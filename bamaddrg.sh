#!/bin/bash

# Initialize an empty string
concatenated_string=""
datapath=/pub64/mattm/read_mapping/Ob/dedup
REF=ref/Ob_genome.fa
base=Ob_freebayes


# Loop through a set of elements
for file in "$datapath"/*; do 
        if [[ "$file" =~ \.bai$ ]]; then
                continue
        fi
        bam="-b ${file} "
        # Concatenate the current element to the string
        concatenated_string="${concatenated_string}${bam}"
done

# Print the final concatenated string
echo $concatenated_string


bamaddrg $concatenated_string > "${base}_bamaddrg.bam"
