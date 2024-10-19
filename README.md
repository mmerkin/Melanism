# Melanism

The code and files used to compare melanic haplotypes across three species of geometrid moth.

Any words in all caps correspond to variables that need to be replaced with your data, unless otherwise specified.

Conda environment yaml files can be loaded with this command:

```bash
conda env create -n NAME -f FILE.yml
```

# Long read assembly

Required files: HiFi reads in fastq.gz format, HiC data in fastq.gz format, parental short reads for trio individuals in fastq.gz format

Required installations: yak, hifiasm, gfatools, yahs, compleasm, gas, seqkit
Conda environment: assembly

### HIFI read assembly

1) Count the kmers from both parental read sets. Note that yak needs to be fed the same reads twice for paired end reads. 
```bash
yak count -b36 -t16 -o OUTPUT.yak <(cat PARENT_R1.fastq.gz PARENT_R2.fastq.gz) <(cat PARENT_R1.fastq.gz PARENT_R2.fastq.gz)
```
2) Perform triobinning assembly with hifiasm
```bash
hifiasm -o OUTPUT.asm -t32 -1 PATERNAL.yak -2 MATERNAL.yak HIFI.fastq.gz
```
3) Perform HIFI-only assembly for individuals lacking paternal reads
```bash
hifiasm -o OUTPUT.asm -t 32 HIFI.fastq.gz
```
4) Convert the output gfa to a fasta file
```bash
gfatools gfa2fa HAPLOTYPE.p_ctg.gfa > HAPLOTYPE.fa
```
5) Calculate error scores for trio samples
```bash
yak trioeval -t16 PATERNAL.yak MATERNAL.yak HAPLOTYPE.fa
```

### Reference genome construction

1) Filter and map HIC data using the [arima mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline)
2) Scaffold the contigs
```bash
yahs HAPLOTYPE.fa MAPPED_HIC.bam
```
3) Perform further scaffolding using haphic, which needs to be installed locally. This worked for the scalloped hazel with its 31 chromosomes, but the pale brindled beauty produced a single scaffold, so the scaffolds in step #2 were used instead
```bash
/PATH/TO/haphic pipeline HAPLOTYPE.fa MAPPED_HIC.bam CHROMOSOME_NUMBER --correct_nrounds 10 --threads 32
```
4) Use seqkit to extract the 31 chromosomes of the scalloped hazel
```bash
seqkit head -n 31 SCAFFOLDS.fa > GENOME.fa 
```
5) HIC-contact maps can be constructed using juicer from the yahs output, which is explained on the [yahs github](https://github.com/c-zhou/yahs?tab=readme-ov-file#generate-hic-contact-maps)

### Sequencing statistics

1) Calculate compleasm score of genome completeness (similar to BUSCO)
```bash
compleasm run -a HAPLOTYPE.fa -o HAPLOTYPE_output -t 16 -l lepidoptera
```
2) Calculate more sequencing statistics with gaas
```bash
gaas_fasta_statistics.pl -f HAPLOTYPE.fa
```

Final sequencing statistics can be found [here](assembly_info)

# Read quality control

Conda environment: filtering

1) Use fastqc to generate reports for each file
```bash
fastqc *
```
2) Use multiqc to merge these into a single report
```bash
multiqc *_fastqc.zip
```

# Pangenome-based variant calling

Required files: 

Required installations:

Conda environment:

1) Create a file 
2) Use minigraph-cactus to create the pangenom
```bash
cactus-pangenome job_store SPECIES_contigs.txt --outDir output_files --outName SPECIES_pangenome --reference reference --gbz clip filter full --viz --gfa clip full --vcf --vcfReference reference --logFile ./SPECIES_pangenome.log --consCores 8
```


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
3) Edit the script [run_generate_pairwise_fst.sh](scripts/run_generate_pairwise_fst.sh) using a text editor, such as nano
4) Replace the three variables at the top with your data, ensuring you are consistent with the files made in step #1
5) Run the script
6) Download the output files and save them into your R working directory:

SPECIESID_pairwise_comparisons.txt contains a list of all the pairwise comparisons

SPECIESID_mean_fsts.txt contains the mean whole genome fst values for each pairwise comparison

SPECIESID_weighted_fsts.txt contains the weighted whole genome fst values for each pairwise comparison

### Plot IBD

# Sweep signatures
