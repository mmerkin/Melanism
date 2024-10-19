# Melanism

The code and files used to compare melanic haplotypes across three species of geometrid moth

Conda environment yaml files can be loaded with this command:

```bash
conda env create -n NAME -f FILE.yml
```

# Long read assembly

Required files: HiFi reads in fastq.gz format, HiC data in fastq.gz format

Required installations: yak,

1) Count 

./yak count -b37 -t16 -o pat.yak <(cat /pub64/mattm/Ob_typica/Ob_pat_reads/16-Ob_05_Pmale_AATCCGGAAT-AACTGTAGGT_L004_R1_001.fastq.gz /pub64/mattm/Ob_typica/Ob_pat_reads/16-Ob_05_Pmale_AATCCGGAAT-AACTGTAGGT_L004_R2_001.fastq.gz) <(cat /pub64/mattm/Ob_typica/Ob_pat_reads/16-Ob_05_Pmale_AATCCGGAAT-AACTGTAGGT_L004_R1_001.fastq.gz /pub64/mattm/Ob_typica/Ob_pat_reads/16-Ob_05_Pmale_AATCCGGAAT-AACTGTAGGT_L004_R2_001.fastq.gz)

# Read quality control

WIP

# Pangenome-based variant calling

# Association analysis

# De novo variant calling

### Map and filter reads

Required files: reference genome, short reads in fastq.gz format

Required installations: bwa-mem2, samtools

1) Index the refence genome
```bash
bwa-mem2 index <in.fasta>
```
2) Edit the script run_map_reads.sh by replacing the variables at the top
3) Run the script to produce filtered bam and bam.bai files

### Produce vcf file

WIP

# vcf filtering

# Principle component analysis

# Isolation-by-distance (IBD)

### Calculate FST

Required files: vcf, coordinates of each sample site

Required installations: vcftools

1) Create files listing the members of each population in the format SPECIESID_POP_samples.txt as in the example. A list of all sample names can be created with vcftools using
```bash
vcf-query -l VARIANTS.vcf > all_samples.txt
```
2) Make a list of all populations to be compared separated by either a tab or space
3) Edit the script [run_generate_pairwise_fst.sh](run_generate_pairwise_fst.sh) using a text editor, such as nano
4) Replace the three variables at the top with your data, ensuring you are consistent with the files made in step #1
5) Run the script
6) Download the output files and save them into your R working directory:

SPECIESID_pairwise_comparisons.txt contains a list of all the pairwise comparisons

SPECIESID_mean_fsts.txt contains the mean whole genome fst values for each pairwise comparison

SPECIESID_weighted_fsts.txt contains the weighted whole genome fst values for each pairwise comparison

### Plot IBD

# Sweep signatures
