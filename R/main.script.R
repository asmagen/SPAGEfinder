
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

  get_slurm_output(job.ids)
  cleanup_files(job.ids)

  # Merge Log-Rank analysis files
  job.ids = merge.clinical.results (r.package.path,results.path,p.val.quantile.threshold,large.queues,large.memory,large.walltime)
  get_slurm_output(job.ids)
  cleanup_files(job.ids)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Molecular filter and base cox testing
  # ____________________________________________________________________________________________________________________________________________________________
  
  job.ids = calculate.base.cox.model (r.package.path,results.path,queues,num.jobs,memory,walltime)
  get_slurm_output(job.ids)
  cleanup_files(job.ids)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Confounder controlled Cox analysis including shuffled analysis for FDR control
  # ____________________________________________________________________________________________________________________________________________________________

  job.ids = calculate.candidates.cox.fdr (r.package.path,results.path,queues,num.jobs,memory,walltime)
  get_slurm_output(job.ids)
  cleanup_files(job.ids)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Compile final SPAGEs list
  # ____________________________________________________________________________________________________________________________________________________________

  states = get.final.SPAGEs (r.package.path,results.path,FDR.threshold = 0.99,apply.PPI.filter = F) # Obtain the final set of SPAGEs based on FDR threshold = 0.99 and without applying the PPI filter
  head(states,20)
  dim(states)

  print(round(table(states[,3]*sign(states[,4]))/nrow(states),2)) # Summarize SPAGE type distribution

  write.csv(states,file = file.path(results.path,'SPAGEs.csv'))
