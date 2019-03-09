
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Set up environment
  # ____________________________________________________________________________________________________________________________________________________________

  source(file.path(scripts,'main.script.functions.R'))

  results.path = file.path(temp,dataset)
  
  create.folders( results.path )
  
  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Preprocess genomic and clinical data
  # ____________________________________________________________________________________________________________________________________________________________

  # Do not use if using example dataset from repository
  # preprocess.genomic.data (data,temp,dataset,results.path)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Run Log-Rank analysis
  # ____________________________________________________________________________________________________________________________________________________________

  job.ids = run.distributed.pairwise.significance (r.package.path,temp,dataset,results.path,queues,num.jobs,memory,walltime)

  # Run calculate.null.molecular in the meanwhile - it's independent
  
  null.molecular.file = file.path(results.path,'null.molecular.RData')
  if( !file.exists(null.molecular.file) ) {
    calculate.null.molecular(temp,dataset,null.molecular.file,shuffle.frac = 0.5)
  }

  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # Merge Log-Rank analysis files
  setwd(results.path)
  job.ids = merge.clinical.results (temp,r.package.path,dataset,p.val.quantile.threshold,results.path,results.path = file.path( results.path,'results' ),large.queues = 'large',memory = '120GB',walltime = '3:00:00')
  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Molecular filter and base cox testing
  # ____________________________________________________________________________________________________________________________________________________________
  
  molecular = F
  job.ids = calculate.base.cox.model (r.package.path,temp,dataset,results.path,queues,num.jobs,memory,walltime)
  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Confounder controlled Cox analysis including shuffled analysis for FDR control
  # ____________________________________________________________________________________________________________________________________________________________

  job.ids = calculate.candidates.cox.fdr (r.package.path,temp,dataset,results.path,queues,num.jobs,memory,walltime)
  setwd(job.ids$path)
  tryCatch({get_slurm_out(job.ids, outtype = 'table',wait = T);get_slurm_out(job.ids, outtype = 'table');get_slurm_out(job.ids, outtype = 'table')},error=function(v) v)
  tryCatch({cleanup_files(job.ids);cleanup_files(job.ids);cleanup_files(job.ids)},error=function(v) v)

  # ____________________________________________________________________________________________________________________________________________________________
  #
  # Compile final GI list based on FDR quantile
  # ____________________________________________________________________________________________________________________________________________________________

  states = get.final.GIs (temp,data,results.path = results.path,LLR.threshold = 0.99,PPI = T)
  head(states,20)
  dim(states)

  print(round(table(states[,3]*sign(states[,4]))/nrow(states),2))

  save(states,file = file.path(results.path,'final.GIs.RData'))
  write.csv(states,file = file.path(results.path,'final.GIs.csv'))
