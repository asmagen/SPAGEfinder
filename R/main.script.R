
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Set up environment
  # ____________________________________________________________________________________________________________________________________________________________

  # User needs to set the following environment-specific parameters:
  # USER_NAME - the user name used to login to the HPC server
  # r.package.path - Path for the downloaded pipeline scripts and data
  # temp - Path for analysis results
  # 
  # 

  # Connect to HPC cluster
  ssh USER_NAME@SERVER_ADDRESS
  screen -r -d GI
  module purge;module add R/common/3.3.1;R
  install.packages(pkgs = c('Rcpp','survival','rslurm','foreach','doMC','data.table','igraph','survcomp'))

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Configure analysis paths and parameters
  # ____________________________________________________________________________________________________________________________________________________________

  scripts = file.path(r.package.path,'scripts')
  data    = file.path(r.package.path,'data') # Path for data folder
  outputs = file.path(r.package.path,'results') # Path for HPC job run outputs
  temp   # Path for analysis results set by user
  
  dataset = 'pancancer.drivers' # dataset/project name to be used as an analysis folder
  p.val.quantile.threshold = 0.8 # Log-Rank threshold (p value quantile)
  
  queues   = 'high_throughput' # SLURM HPC queue
  memory   = '8GB' # Memory allocation per job
  walltime = '8:00:00' # Time limit per job
  num.jobs = 100 # Number of jobs

  source(file.path(scripts,'main.script.functions.R'))

  base.res.path = file.path(temp,dataset)
  
  create.folders( base.res.path )
  
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Preprocess genomic and clinical data
  # ____________________________________________________________________________________________________________________________________________________________

  # Do not use if using example dataset from repository
  # preprocess.genomic.data (data,temp,dataset,base.res.path)

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
