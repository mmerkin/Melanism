# Melanism
Files used in my melanism project

# Long read assembly

WIP


# De novo variant calling

### Map and filter reads

Required files: reference genome, short reads in fastq format

Required installations: bwa-mem2, samtools

1) Index the refence genome
```bash
bwa-mem2 index <in.fasta>
```
2) Edit the script run_map_reads.sh by replacing the variables at the top
3) Run the script to produce filtered bam and bam.bai files

### Produce vcf file

WIP

# Isolation-by-distance (IBD)

### Calculate FST

Required files: vcf, coordinates of each sample site

Required installations: vcftools

1) Ensure vcftools is installed
2) Create files listing the members of each population in the format SPECIESID_POP_samples.txt as in the example. A list of all sample names can be created with vcftools using
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
