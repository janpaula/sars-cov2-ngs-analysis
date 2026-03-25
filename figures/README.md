# Figures
This directory contains the graphical outputs generated during the analysis of the SARS-CoV-2 NGS dataset.

All figures were produced in R using the packages ggplot2, Gviz, and VariantAnnotation, based on the processed variant and coverage data.

These visualizations summarize key aspects of the sequencing analysis, including variant distribution, sequencing depth, and allele frequency patterns across the viral genome.

# Figure List

## Genome Coverage Plot
This plot shows the sequencing depth across the SARS-CoV-2 genome.

It allows evaluation of:

 - regions with high or low coverage
 - uniformity of sequencing
 - potential gaps in read alignment

Coverage values were derived from the BAM alignment file using SAMtools.

## Variant Manhattan Plot

This plot displays the genomic positions of detected variants across the viral genome.

Each point represents a variant identified by the variant calling pipeline.

The Manhattan-style visualization helps identify:

 - clusters of mutations
 - genomic regions with variant accumulation
 - overall distribution of detected substitutions
 

## Allele Frequency Distribution

This histogram shows the distribution of allele frequencies for the detected variants.

Allele frequency (AF) is calculated as:

$$
AF = \frac{AO}{DP}
$$

Where:

- **AO** = number of reads supporting the alternative allele  
- **DP** = total read depth at that position

Running this script will reproduce the visualizations stored in this directory.


