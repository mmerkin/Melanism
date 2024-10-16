# Melanism
Files used in my melanism project



# Isolation-by-distance

Requirements: vcf, coordinates of each sample site

1) Create files listing the members of each population in the format SPECIESID_POP_samples.txt
2) Make a list of all populations to be compared separated by either a tab or space
3) Edit the script "run_generate_pairwise_fst.sh" using a text editor, such as nano
4) Replace the three variables at the top with your data, ensuring you are consistent with the files made in step #1
5) Run the script with bash run_generate_pairwise_fst.sh
6) Download the output files:
SPECIESID_pairwise_comparisons.txt contains a list of all the pairwise comparisons
SPECIESID_mean_fsts.txt contains the mean whole genome fst values for each pairwise comparison
SPECIESID_weighted_fsts.txt contains the weighted whole genome fst values for each pairwise comparison
