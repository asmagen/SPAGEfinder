# SPAGE-finder Manual (*Magen et al*)

This repository contains code and documentation for a multi-step computational pipeline to search for Survival associated Pairwise Gene Expression states finder (SPAGE-finder).
Previously, the repository was called EnGIne. The manuscript describing the project can be found at:
https://www.biorxiv.org/content/10.1101/253120v2

The repository can be retrieved from GitHub via the command
*git clone https://github.com/asmagen/SPAGEfinder*

This creates a new directory called SPAGEfinder where this README.md can be found.

SPAGEfinder has two subdirectories: 
1. data
2. R  

# Code files
There are seven code files:
- aggregateLogRankGene.cpp
- analyze.pairwise.significance.R  
- calculate.base.cox.model.R  
- calculate.candidates.cox.fdr.R  
- main.script.functions.R  
- main.script.R 
- merge.pancancer.results.R

Six of these are in R and one in C++.
The C++ code is compiled from within R. Therefore, one need not compile it before running the pipeline.

## Input format
The input includes mRNA expression matrix and patient clinical-demographic information which are stored in a single data source object and corresponding file. The TCGA dataset (filtered to Census Cancer Genes) is provided here for quick analysis (Example dataset available: 'data/data.mRNA.RData').

## Output format
Final SPAGEs list is generated in a matrix format where each row represents a SPAGE that is annotated by a quadruple (x,y,bin,effect), where x and y are the two interacting genes, bin is a number indicating the bin annotation and effect annotated the significance level where the sign of the effect represents the direction of the interaction; Positive sign represents higher survival risk while negative sign represents lower survival risk.

## Analysis setup

### UNIX commands
Connect to remote server [Example: ssh USER_NAME@SERVER_ADDRESS]  
git clone https://github.com/asmagen/SPAGEfinder.git  
Invoke R version 3.3.1 [Example (may be different across systems): module purge; module add R/3.3.1; R]  
The following commands are set and run in the R environment.

### Install R packages
Install the required R packages into the default location (no need to specify where to install).  
install.packages(pkgs = c('Rcpp','RcppArmadillo','survival','rslurm','foreach','doMC','data.table','igraph','whisker','foreach'))  
Specify a repository of choice and verify successful installations.
source("https://bioconductor.org/biocLite.R"); biocLite("survcomp")
Verify successful installation.

### Define relevant analysis paths (in R)
*r.package.path* = 'USER_SET_PACKAGE_PATH' # Define path for the downloaded pipeline scripts and data  
*results.path* = 'USER_SET_ANALYSIS_DIRECTORY_PATH' # Define path for analysis results set by user  

### Assign values to additional analysis and slurm parameters

The suggested values shown below may be changed by the user.  

*p.val.quantile.threshold* = 0.8 # Log-Rank threshold (p value quantile)  

*queues*   = 'high_throughput' # SLURM HPC queue  
*memory*   = '8GB' # Memory allocation per job  
*walltime* = '8:00:00' # Time limit per job  
*num.jobs* = 100 # Number of jobs  

Continue analysis by executing the commands in 'R/main.script.R'

### Creating new datasets for analysis
The function preprocess.genomic.data (r.package.path) can be used to process and perform binning of a 'dataset' object located at 'data/dataset.RData') and constructed in the following fields:  
*mRNA* - RSEM normalized mRNA measurements (rows corresponding to genes and columns to samples)  
*scna* - copy-number variation measurements (rows corresponding to genes and columns to samples)  
*samples* - sample IDs as factors  
*type* - cancer types as factors  
*sex* - sex annotation as factors  
*race* - race annotation as factors  
*time* - patient survival as number of days to death  
*status* - patient death/alive status as 0 or 1, respectively  
