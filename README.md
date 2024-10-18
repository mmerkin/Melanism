# Melanism
Files used in my melanism project

# Long read assembly

WIP


# De novo variant calling

Requirements: reference genome, short reads in fastq format

1) Install bwa-mem2
2) Index the refence genome
```bash
bwa-mem2 index <in.fasta>
```
3) Edit the script run_map_reads.sh by replacing the variables at the top
4) Run the script to map, filter and index the reads, producing bam and bam.bai files


# Isolation-by-distance (IBD)

### Calculate FST

Requirements: vcf, coordinates of each sample site

1) Ensure vcftools is installed
2) Create files listing the members of each population in the format SPECIESID_POP_samples.txt. Sample names can be viewed using
```bash
vcf-query -l VARIANTS.vcf > all_samples.txt
```
3) Make a list of all populations to be compared separated by either a tab or space
4) Edit the script "run_generate_pairwise_fst.sh" using a text editor, such as nano
5) Replace the three variables at the top with your data, ensuring you are consistent with the files made in step #1
6) Run the script with bash run_generate_pairwise_fst.sh
7) Download the output files:
SPECIESID_pairwise_comparisons.txt contains a list of all the pairwise comparisons
SPECIESID_mean_fsts.txt contains the mean whole genome fst values for each pairwise comparison
SPECIESID_weighted_fsts.txt contains the weighted whole genome fst values for each pairwise comparison

### Plot IBD
