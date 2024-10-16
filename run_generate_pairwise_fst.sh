#! /bin/bash

# Set variables

# Full path to gzipped vcf
gzvcf_path=~/final_vcfs/Bb_fb_background.vcf.gz 
# Designate species for output
species=Bb
# Space or tab separated list of population names
Pop_list=(BTW CDY HTH KKH MDP OGR WRT WBP)

# Run code

# Generate a vcf file for all pairwise combinations of populations. This does include repeats, so filter those out later

for i in "${Pop_list[@]}"; do 
	for j in "${Pop_list[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		vcftools --gzvcf ${gzvcf_path} --weir-fst-pop ${species}_${i}_samples.txt --weir-fst-pop ${species}_${j}_samples.txt --out ${species}_${i}_vs_${j}
	done
done

# Generate the names column

for i in "${Pop_list[@]}"; do 
	for j in "${Pop_list[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		echo "${i}-${j}"
	done
done > ${species}_pairwise_comparisons.txt

# Generate the mean FST column

for i in "${Pop_list[@]}"; do 
	for j in "${Pop_list[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		grep 'Weir and Cockerham mean Fst estimate:' ${species}_${i}_vs_${j}.log | sed 's/^.*: //'
	done
done > ${species}_mean_fsts.txt

# Generate the weighted FST column

for i in "${Pop_list[@]}"; do 
	for j in "${Pop_list[@]}"; do 
		if [ $i = $j ]; then
			continue
		fi
		grep 'Weir and Cockerham weighted Fst estimate:' ${species}_${i}_vs_${j}.log | sed 's/^.*: //'
	done
done > ${species}_weighted_fsts.txt
