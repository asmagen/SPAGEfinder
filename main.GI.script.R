
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Set up environment
  # ____________________________________________________________________________________________________________________________________________________________

  # Connect to HPC cluster
  ssh magen@cbcbsub00.umiacs.umd.edu
  screen -r -d GI
  module purge;module add R/common/3.3.1;R
  install.packages(pkgs = c('Rcpp','survival','rslurm','foreach','doMC','data.table','igraph','survcomp'))
  # Put script files in scripts folder 'pipeline_example'
  # Put data files (TCGA dataset 'prob.TCGA.extended.RData', Hippie file 'hippie_current.txt' and Census cancer genes 'Census_all.csv') in data folder
  # Create analysis folder ('temp')

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Configure analysis paths and parameters
  # ____________________________________________________________________________________________________________________________________________________________

  source('/cbcbhomes/magen/RSR/pan.cancer/pipeline_example/main.GI.script.functions.R')
  # /cbcbhomes/magen/RSR/pan.cancer/pipeline_example

  scripts = '/cbcbhomes/magen/RSR/pan.cancer/pipeline_example' # Path for pipeline scripts
  data    = '/cbcbhomes/magen/RSR/data' # Path for data folder
  outputs = '/cbcbhomes/magen/RSR/results' # Path for HPC job run outputs
  temp    = '/cbcb/project2-scratch/amagen/RSR/data' # Path for analysis results
  dataset                  = 'pancancer.example' # dataset/project name to be used as an analysis folder
  p.val.quantile.threshold = 0.8 # Log-Rank threshold (p value quantile)
  
  queues   = 'high_throughput' # SLURM HPC queue
  memory   = '8GB' # Memory allocation per job
  walltime = '8:00:00' # Time limit per job
  num.jobs = 100 # Number of jobs

  base.res.path = file.path(temp,dataset)
  
  create.folders( base.res.path )

  library(data.table)
  census = fread(file.path(data,'Census_all.csv'))
  census = data.table(census)
  drivers = census$"Gene Symbol"
  drivers
  
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Preprocess genomic and clinical data
  # ____________________________________________________________________________________________________________________________________________________________

  results.path=base.res.path
  preprocess.genomic.data (data,temp,dataset,base.res.path,constrain.gene.set = drivers)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Run Log-Rank analysis
  # ____________________________________________________________________________________________________________________________________________________________

  job.ids = run.distributed.pairwise.significance (scripts,temp,dataset,results.path,queues,num.jobs,memory,walltime)

  # Run calculate.null.molecular in the meanwhile - it's independent
  results.path = base.res.path
  null.molecular.file = file.path(results.path,paste('null.molecular.RData',sep='.'))
  if( !file.exists(null.molecular.file) ) {
    calculate.null.molecular(temp,dataset,null.molecular.file,shuffle.frac = 0.5)
  }

  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # Merge Log-Rank analysis files
  setwd(base.res.path)
  job.ids = merge.clinical.results (temp,scripts,dataset,p.val.quantile.threshold,base.res.path,results.path = file.path( base.res.path,'results' ),large.queues = 'large',memory = '120GB',walltime = '3:00:00')
  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Molecular filter and base cox testing
  # ____________________________________________________________________________________________________________________________________________________________
  
  molecular = F
  job.ids = calculate.base.cox.model (scripts,data,temp,dataset,results.path,queues,num.jobs,memory,walltime)
  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Confounder controlled Cox analysis including shuffled analysis for FDR control
  # ____________________________________________________________________________________________________________________________________________________________

  support = 0
  job.ids = calculate.candidates.cox.fdr (scripts,data,temp,dataset,results.path,queues,num.jobs,memory,walltime)
  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Compile final GI list based on FDR quantile
  # ____________________________________________________________________________________________________________________________________________________________

  states = get.final.GIs (temp,data,results.path = base.res.path,LLR.threshold = 0.99,PPI = F)
  head(states,20)
  dim(states)

  print(round(table(states[,3]*sign(states[,4]))/nrow(states),2))

  save(states,file = file.path(results.path,'final.GIs.RData'))
  write.csv(states,file = file.path(results.path,'final.GIs.csv'))
