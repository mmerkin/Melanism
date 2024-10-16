#! /bin/bash

#Prerequisites: vcftools installed, population sample files created in the format SPECIESID_POP_samples.txt 


# Set variables

# Full path to .vcf or .vcf.gz file
VCF=/path/to/variants.vcf.gz 
# Designate species ID for output names
SPECIESID=example
# Space or tab separated list of population names
POPLIST=(POP1 POP2 POP3 POP4)








# Run code

# Generate a vcf file for all pairwise combinations of populations after determining whether the vcf is gzipped

if [[ $VCF =~ \.gz$ ]]; then
	for i in "${POPLISTt[@]}"; do 
		for j in "${POPLIST[@]}"; do 
			if [ $i = $j ]; then
				continue
			fi
			vcftools --gzvcf ${VCF} --weir-fst-pop ${SPECIESID}_${i}_samples.txt --weir-fst-pop ${SPECIESID}_${j}_samples.txt --out ${SPECIESID}_${i}_vs_${j}
		done
	done
else
    for i in "${POPLIST[@]}"; do 
		for j in "${POPLIST[@]}"; do 
			if [ $i = $j ]; then
				continue
			fi
			vcftools --vcf ${VCF} --weir-fst-pop ${SPECIESID}_${i}_samples.txt --weir-fst-pop ${SPECIESID}_${j}_samples.txt --out ${SPECIESID}_${i}_vs_${j}
		done
	done
fi

# Generate the names column

for i in "${POPLIST[@]}"; do 
	for j in "${POPLIST[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		echo "${i}-${j}"
	done
done > ${SPECIESID}_pairwise_comparisons.txt

# Generate the mean FST column

for i in "${POPLIST[@]}"; do 
	for j in "${POPLIST[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		grep 'Weir and Cockerham mean Fst estimate:' ${SPECIESID}_${i}_vs_${j}.log | sed 's/^.*: //'
	done
done > ${SPECIESID}_mean_fsts.txt

# Generate the weighted FST column

for i in "${POPLIST[@]}"; do 
	for j in "${POPLIST[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		grep 'Weir and Cockerham weighted Fst estimate:' ${SPECIESID}_${i}_vs_${j}.log | sed 's/^.*: //'
	done
done > ${SPECIESID}_weighted_fsts.txt
