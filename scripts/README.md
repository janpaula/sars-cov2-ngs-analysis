# Scripts

This directory contains the scripts used in the bioinformatics pipeline for the analysis of SARS-CoV-2 sequencing data (2020 dataset).

The scripts implement the main steps of a typical Next Generation Sequencing (NGS) workflow, including read alignment, variant calling, and downstream statistical analysis.

## Overview of the Pipeline

The analysis workflow includes the following steps:

1. Alignment of sequencing reads to the reference genome
2. Processing of alignment files
3. Variant detection
4. Coverage analysis
5. Variant visualization and statistical analysis

## alignment.sh

This script performs the alignment of sequencing reads against the reference genome using **BWA MEM**.

Main steps:

- Index reference genome
- Align sequencing reads to the reference genome
- Generate SAM alignment file

Example command inside the script:

Input

- `bwa mem reference.fasta reads.fastq > alignment.sam`

Output:

- `alignment.sam`

## variant_calling.sh

This script processes the alignment file and performs variant calling.

Main steps:

- Convert SAM to BAM
- Sort BAM file
- Index BAM file
- Call variants using FreeBayes

Typical commands used:

- `samtools view -bS alignment.sam > alignment.bam`
- `samtools sort alignment.bam -o alignment_sorted.bam`
- `samtools index alignment_sorted.bam`
- `freebayes -f reference.fasta alignment_sorted.bam > variants.vcf`

Outputs:

- `alignment_sorted.bam`
- `variants.vcf`

## coverage.R

This R script calculates genome coverage statistics from the alignment file.

Main tasks:

- Read coverage data
- Calculate coverage depth across the genome
- Generate coverage plots

Output:

- `coverage.txt`
- `coverage_plot.png`

## variant_analysis.R

This script analyzes the variant call file (VCF) generated during variant calling.

Main tasks:

- Import VCF data
- Extract variant positions
- Calculate allele frequencies
- Generate visualizations

Figures generated include:

- Variant distribution plot
- Manhattan plot
- Allele frequency histogram

## Dependencies

The following software and R packages are required:

### Command-line tools

- BWA
- Samtools
- FreeBayes

### R packages

- ggplot2
- VariantAnnotation
- Gviz
- kableExtra

## Notes

These scripts were developed as part of a bioinformatics analysis of SARS-CoV-2 sequencing data from 2020, demonstrating a standard pipeline for viral genome variant detection.
