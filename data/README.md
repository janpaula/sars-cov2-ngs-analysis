# Data

This directory contains the reference genome used for the analysis.

## Reference Genome

AF086833.fasta

This FASTA file contains the reference genome used for read alignment and variant detection.

The reference genome is used as the template for mapping sequencing reads and identifying nucleotide variants.

## Notes

The reference genome must be indexed before alignment using BWA:

bwa index AF086833.fasta
