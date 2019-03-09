# SPAGE-finder Manual (Brief version)

This repository contains code and coumentation for a multi-step computational pipeline to search for Survival associated Pairwise Gene Expression states finder (SPAGE-finder).
Previously, the repository was called EnGIne. The manuscript describing the project can be found at:
https://www.biorxiv.org/content/10.1101/253120v2

The repository can be retrieved from GitHub via the command
git clone https://github.com/asmagen/SPAGEfinder

This creates a new directory called SPAGEfinder where this README.md can be found.

SPAGEfinder has two subdirectories: 
data
R  

# Analysis pipeline summary
See main.script.R for running the analysis.

# Code files
There are seven code files:
aggregateLogRankGene.cpp
analyze.pairwise.significance.R  
calculate.base.cox.model.R  
calculate.candidates.cox.fdr.R  
main.script.functions.R  
main.script.R 
merge.pancancer.results.R

Six of these are in R and one in C++
The C++ code is compiled from within R. Therefore, one need not compile it before running the pipeline.

## Input format
The input includes mRNA expression matrix and patient clinical-demographic information which are stored in a single data source object and corresponding file. The TCGA dataset (filtered to Census Cancer Genes) is provided here for quick analysis (Example dataset available: 'data/data.mRNA.RData').

## Output format
Final SPAGEs list is generated in a matrix format where each row represents a SPAGE that is annotated by a quadruple (x,y,bin,effect), where x and y are the two interacting genes, bin is a number indicating the bin annotation and effect annotated the significance level where the sign of the effect represents the direction of the interaction; Positive sign represents higher survival risk while negative sign represents lower survival risk.

## Analysis setup

### UNIX commands
Connect to remote server [Example: ssh USER_NAME@SERVER_ADDRESS]  
Invoke R version 3.3.1 [Example (may be different across systems): module purge; module add R/3.3.1; R]  
The following commands are set and run in the R environment.

### Install R packages
install.packages(pkgs = c('Rcpp','survival','rslurm','foreach','doMC','data.table','igraph'))
Specify a repository of choice.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("survcomp", version = "1.24.0")

### Define relevant analysis paths (in R)
r.package.path = 'USER_SET_PATH' # Define path for the downloaded pipeline scripts and data  
temp = 'USER_SET_PATH' # Define path for analysis results set by user  

### Change analysis parameters if defaults are inappropriate for the dataset  
dataset = 'pancancer.drivers'# dataset/project name to be used as an analysis folder  
p.val.quantile.threshold = 0.8 # Log-Rank threshold (p value quantile)  

queues   = 'high_throughput' # SLURM HPC queue  
memory   = '8GB' # Memory allocation per job  
walltime = '8:00:00' # Time limit per job  
num.jobs = 100 # Number of jobs  

Continue analysis by executing the commands in 'R/main.script.R'

### Creating new datasets for analysis
The function preprocess.genomic.data (r.package.path) can be used to process and perform binning of a 'dataset' object located at 'data/dataset.RData') and constructed in the following format:  
mRNA field containing RSEM normalized mRNA measurements (rows corresponding to genes and columns to samples)  
scna field containing copy-number variation measurements (rows corresponding to genes and columns to samples)  
samples field containing sample IDs as factors  
type field containing cancer types as factors  
sex field containing sex annotation as factors  
race field containing race annotation as factors  
time field containing patient survival as number of days to death  
status field containing patient death/alive status as 0 or 1, respectively  
