# SPAGE-finder Manual (Brief version)

A multi-step computational pipeline to search for Survival associated Pairwise Gene Expression states finder (SPAGE-finder).

# Analysis pipeline summary
See main.script.R for running the analysis.

## Input format
The input includes mRNA expression matrix and patient clinical-demographic information. The TCGA dataset (filtered to Census Cancer Genes) is provided here for quick analysis. Example dataset available in 'data' folder.

## Output format
Final SPAGEs list is generated in a matrix format where each row represents a SPAGE that is annotated by a quadruple (x,y,bin,effect), where x and y are the two interacting genes, bin is a number indicating the bin annotation and effect annotated the significance level where the sign of the effect represents the direction of the interaction; Positive sign represents higher survival risk while negative sign represents lower survival risk.

## C Code
C code is being compiled by R thus compiling it before running the pipeline is not needed.
