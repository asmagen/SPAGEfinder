# EncyclopediaGeneticInteractions

A multi-step computational pipeline (Encyclopedia of Genetic Interactions - EnGIne) to search for clinically relevant functional gene interactions given their association with patient survival.

# Analysis pipeline summary
See main.GI.script.R for full details.

  ## Set up environment
  ### (Connecting to high performance computing environment and installing packages; Slurm Workload Manager required)
  ssh username@server
  
  module purge;module add R/common/3.3.1;R
  
  install.packages(pkgs = c('Rcpp','survival','rslurm','foreach','doMC','data.table','igraph','survcomp'))

  ## Preprocess genomic and clinical data
  ### (Loading mRNA and clinical data, identifying abundant clinical-demographic subgroups and perform mRNA binning)
  preprocess.genomic.data (data,temp,dataset,base.res.path,constrain.gene.set = drivers)
  
  ## Run Log-Rank analysis
  ### (Distributed C++ code for all 9 bins survival trends modelling)
  run.distributed.pairwise.significance (scripts,temp,dataset,results.path,queues,num.jobs,memory,walltime)

  ## Merge Log-Rank analysis files
  merge.clinical.results (temp,scripts,dataset,p.val.quantile.threshold,base.res.path,results.path = file.path( base.res.path,'results' ),large.queues = 'large',memory = '120GB',walltime = '3:00:00')

  ## Run Null molecular calculation
  ### (Calculating the expected number of samples per bin to adjust to bin position biases)
  calculate.null.molecular(temp,dataset,null.molecular.file,shuffle.frac = 0.5)

  ## Molecular filter and base cox testing
  ### (Performing distributed molecular depletion/enrichment analysis based on Null molecular and filtering potential interactions with inconsistent molecular-survival trends)
  calculate.base.cox.model (scripts,data,temp,dataset,results.path,queues,num.jobs,memory,walltime)

  ## Confounder controlled Cox analysis including shuffled analysis for FDR control
  ### (Performing distributed Cox Proportional hazards regression models (survival modelling) controlled for clinical-demographic confounders; Including shuffled analysis for FDR correction)
  calculate.candidates.cox.fdr (scripts,data,temp,dataset,results.path,queues,num.jobs,memory,walltime)
  
  ## Compile final GI list based on FDR quantile
  ### (Comparing Cox significance to shuffled analysis for FDR filter; Including PPI filtering option)
  states = get.final.GIs (temp,data,results.path = base.res.path,LLR.threshold = 0.99,PPI = F)
  
  write.csv(states,file = file.path(results.path,'final.GIs.csv'))
