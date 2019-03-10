
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Set up environment
  # ____________________________________________________________________________________________________________________________________________________________

  source(file.path(r.package.path,'R','main.script.functions.R'))
  
  create.folders( results.path )
  
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Preprocess genomic and clinical data
  # ____________________________________________________________________________________________________________________________________________________________

  # Do not use if using example dataset from repository
  # preprocess.genomic.data (r.package.path)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Run Log-Rank analysis
  # ____________________________________________________________________________________________________________________________________________________________

  job.ids = run.distributed.pairwise.significance (r.package.path,results.path,queues,num.jobs,memory,walltime)

  # Run calculate.null.molecular in the meanwhile - it's independent
  
  calculate.null.molecular(r.package.path,results.path)

  setwd(job.ids$path)
  get_slurm_out(job.ids,wait = T)
  cleanup_files(job.ids)

  # Merge Log-Rank analysis files
  setwd(results.path)
  job.ids = merge.clinical.results (r.package.path,results.path,p.val.quantile.threshold,large.queues = 'large',memory = '120GB',walltime = '3:00:00')
  setwd(job.ids$path)
  get_slurm_out(job.ids,wait = T)
  cleanup_files(job.ids)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Molecular filter and base cox testing
  # ____________________________________________________________________________________________________________________________________________________________
  
  job.ids = calculate.base.cox.model (r.package.path,results.path,queues,num.jobs,memory,walltime)
  setwd(job.ids$path)
  get_slurm_out(job.ids,wait = T)
  cleanup_files(job.ids)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Confounder controlled Cox analysis including shuffled analysis for FDR control
  # ____________________________________________________________________________________________________________________________________________________________

  job.ids = calculate.candidates.cox.fdr (r.package.path,results.path,queues,num.jobs,memory,walltime)
  setwd(job.ids$path)
  get_slurm_out(job.ids,wait = T)
  cleanup_files(job.ids)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Compile final GI list based on FDR quantile
  # ____________________________________________________________________________________________________________________________________________________________

  states = get.final.GIs (temp,data,results.path,LLR.threshold = 0.99,PPI = F)
  head(states,20)
  dim(states)

  print(round(table(states[,3]*sign(states[,4]))/nrow(states),2))

  write.csv(states,file = file.path(results.path,'SPAGEs.csv'))
