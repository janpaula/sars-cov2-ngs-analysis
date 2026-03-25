# SARS-CoV-2 Variant Analysis using NGS Data

## Overview

This repository contains a bioinformatics workflow for the analysis of Next Generation Sequencing (NGS) data from SARS-CoV-2.

The dataset analyzed corresponds to viral genomic data generated in 2020, during the early phase of the COVID-19 pandemic.

The objective of this project is to identify nucleotide substitutions in the viral genome using a reference-based alignment and variant calling pipeline, as well as to evaluate sequencing coverage and allele frequency distribution across the genome.

This project demonstrates a practical application of bioinformatics tools for viral genome analysis.

## Reference Genome

The analysis uses the following reference genome:

AF086833.fasta

## Tools and Software

The analysis pipeline was performed using the following tools:

- BWA – read alignment
- Samtools – file processing and indexing
- FreeBayes – variant calling
- R – statistical analysis and visualization

R packages used:

- ggplot2
- Gviz
- VariantAnnotation
- kableExtra

## Workflow

The analysis pipeline consists of the following steps:

### 1. Reference genome indexing
The SARS-CoV-2 reference genome was indexed using BWA.

### 2. Read alignment
Sequencing reads were aligned to the reference genome using BWA MEM.

### 3. SAM/BAM processing
Alignment files were converted, sorted, and indexed using Samtools.

### 4. Variant calling
Variants were identified using FreeBayes.

### 5. Variant analysis
The VCF file was analyzed in R to evaluate:

- variant positions
- sequencing depth
- allele frequency
- variant distribution across the genome
  
### 6. Visualization

The following visualizations were generated:

- Manhattan plot of genomic variants
- Genome coverage plot
- Allele frequency distribution
 
## Repository Structure

sars-cov2-ngs-analysis/

data/
Reference genome

scripts/
Pipeline scripts and R analysis scripts

results/
Variant files and processed outputs

figures/
Generated plots

report/
PDF report containing the full analysis


## Work flow

FASTQ
 ↓
BWA alignment
 ↓
SAM/BAM processing
 ↓
FreeBayes variant calling
 ↓
VCF analysis in R
 ↓
Visualization


## Example Variant Calls

| Position | Reference | Alternative | Depth | Caller |
|--------|--------|--------|--------|--------|
| 5876 | C | T | 2 | FreeBayes |
| 5954 | T | C | 2 | FreeBayes |

## Example Figures

The repository includes graphical summaries of the analysis:

- Genome coverage plot
- Variant Manhattan plot
- Allele frequency distribution

## Future Improvements

Possible improvements for this pipeline include:

- Functional annotation of variants
- Comparison with global SARS-CoV-2 mutation databases
- Phylogenetic analysis of detected variants

## Author

Jânice Roberta de Paula

Bioinformatics researcher interested in:

- genomics
- neurological genetics
- precision medicine
