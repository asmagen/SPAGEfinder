# SPAGE-finder Manual (*Magen et al*)

This repository contains code and documentation for a multi-step computational pipeline to search for Survival associated Pairwise Gene Expression states finder (SPAGE-finder).
Previously, the repository was called EnGIne. The manuscript describing the project can be found at:
https://www.biorxiv.org/content/10.1101/253120v2

The repository can be retrieved from GitHub via the command
*git clone https://github.com/asmagen/SPAGEfinder*

This creates a new directory called SPAGEfinder where this README.md can be found. If one reads README.md as a text file, there will be asterisks around the variable names; do not copy the asterisks when copying a command to assign a value to a variable.

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
Install the required R packages into the default location (no need to specify where to install, enter 'yes' to indicate installation to personal library if asked).  
install.packages(pkgs = c('Rcpp','RcppArmadillo','survival','rslurm','foreach','doMC','data.table','igraph','whisker','foreach'))  
Specify a repository of choice and verify successful installations by loading packages (Example: library('Rcpp')).  
source("https://bioconductor.org/biocLite.R"); biocLite("survcomp")  
Verify successful installation.

### Define relevant analysis paths (in R)
*r.package.path* = 'USER_SET_PACKAGE_PATH' # Define path for the downloaded pipeline scripts and data (Example: '/USER/SPAGEfinder')  
*results.path* = 'USER_SET_ANALYSIS_DIRECTORY_PATH' # Define path for analysis results set by user (Example: '/USER/analysis/TCGA_analysis')  

### Assign values to additional analysis and slurm parameters

The suggested values shown below may be adjusted by the user as needed.  

*p.val.quantile.threshold* = 0.8 # Log-Rank threshold (p value quantile)  

The next 7 parameters may need to be adjusted based on the specification of your high-performance computing system. Run the following command to obtain the info about the available queues, memory, walltime and num.jobs (number of concurrent jobs) resources:  
*sacctmgr show qos format=name,MaxJobs,MaxWall,MaxTRES*  

	    	   Name  MaxJobs    MaxWall       MaxTRES  
			---------- ------- ----------- -------------  
				 normal                                    
				default      16    01:00:00        mem=4G  
		 throughput     125    18:00:00       mem=36G  
 high_throghput     300    08:00:00        mem=8G  
					large       5 11-00:00:00      mem=128G  
				 xlarge       1 21-00:00:00      mem=512G  
					 long      16  7-00:00:00       mem=12G  
    workstation       4  7-00:00:00       mem=48G  

Based on this information and the scope of the analysis (whole-genome or in this example case, only about 500 genes) you would define the following parameters:

*queues*   = 'throughput' # SLURM HPC queue  
*num.jobs* = 50 # Number of concurrent jobs  
*walltime* = '1:00:00' # Time limit per job  
*memory*   = '4GB' # Memory allocation per job  

And the following parameters specifically for the merge.pancancer.results function as it requires more memory than the usual:
*large.queues*   = 'throughput' # SLURM HPC queue   
*large.walltime* = '1:00:00' # Time limit per job  
*large.memory*   = '36GB' # Memory allocation per job  

The appropriate parameters for whole genome analysis (analysis of about 20k genes) are:
*num.jobs* = 120, *walltime* = '8:00:00', *memory* = '8GB', *large.queues* = 'large', *large.walltime* = '5:00:00', *large.memory* = '120GB'  

Note that in some systems there is no need to specify the queues parameter. If the queues specification above results in an error, use *queues = NA* and *large.queues = NA* to let the system choose the appropriate queues by itself.

Continue the analysis by executing the commands in 'R/main.script.R'  

### Creating new datasets for analysis
The function preprocess.genomic.data (r.package.path) can be used to process and perform binning of a 'dataset' object located at 'data/dataset.RData') and constructed in the following fields:  
- *mRNA* - RSEM normalized mRNA measurements (rows corresponding to genes and columns to samples)  
- *scna* - copy-number variation measurements (rows corresponding to genes and columns to samples)  
- *samples* - sample IDs as factors
- *type* - cancer types as factors
- *sex* - sex annotation as (two) factors
- *race* - race annotation as factors
- *time* - patient survival as number of days to death  
- *status* - patient death/alive status as 0 or 1, respectively  

The annotation or format of *samples, type, sex and race* is not important as long as the variables are converted to factors.

### Potential problems

- The following error may rarely come up in the *get_slurm_output* function:  
	'slurm_load_jobs error: Socket timed out on send/recv operation'  
	The error does not reflect the failure of the analysis but only a failure of monitoring the job execution. Simply rerun the *get_slurm_output* function.  
- The following error may come up due to insufficient resource allocation.  
	'The following files are missing: ... Check failed jobs error outputs'  

