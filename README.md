 
Single Meiosis project: 

- Bioinformatics and statistics codes to analyze the 20 tetrads of Arabidopsis thaliana from [Ref article] 
- Supplementary methods, results and code to reproduce the figures.
- Final datasets

# Input data:

Tetrads sequencing products are available in Biostudies: 
[E-MTAB-14435](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14435)
Each tetrad's sequencing product is made of 8 files (2 paired-end for each 
4 individual). 

Example:

```bash
|_tetrad:
   ├─ tetrad_M1_s1.fq
   ├─ tetrad_M1_s2.fq
   ├─ tetrad_M2_s1.fq
   ├─ tetrad_M2_s2.fq
   ├─ tetrad_M3_s1.fq
   ├─ tetrad_M3_s2.fq
   ├─ tetrad_M4_s1.fq
   ├─ tetrad_M4_s2.fq
``` 


# Tetrad's processing:

## scripts_tetrad_processing

Contain the main scripts necessary to run tetrads sequencing products 
analysis.
 
- singlemeiosis.sh (main script)
- singlemeiosis-pipeline_genotoul_V3 (conf file)
- Launch_Tetrad_Analysis_genotoul_V3.R (sub-script)
 

## bash: 

### single-meiosis-pipeline

- Align tetrads on Col + Ler genome 
- Separate bam by parental genomes (Col and Ler) and call read depth by genome position.

input:

```bash
   |_Tetrad (*.fq)
   |_Genomes
      |_ Col.fasta
      |_ Ler.fasta 
``` 


output:

```bash
   |_Tetrad
     |_log
     |_M1
        |_ M1_Coller.bam
        |_ M1_Col.bam
        |_ M1_Ler.bam
        |_ M1_Col.vcf
        |_ M1_Ler.vcf
     |_M2
     |_M3
     |_M4
    config
```

## R: variantutils 

- Process the output of the single-meiosis bash pipeline
and prepare input for hmm-nco function.
- Call the crossover.py programm

## R: hmm-nco 

author: S. Robin 

R function to run the HMM model for Tetrad's genotyping .

# Report

## Rmd: 

Supplementary methods, results and code to reproduce the figures.

## Rmd/Data:

Final datasets


