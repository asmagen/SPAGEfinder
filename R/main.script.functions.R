clear.previous.results <- function( results.path ) {

  system( paste('rm -r ',results.path,'/*',sep = ''),
          intern = TRUE,
          ignore.stderr = TRUE
        )
  dir.create( results.path, 
    recursive = TRUE, showWarnings = FALSE )
}

create.folders <- function( results.path ) {

  dir.create( 
    results.path, 
    recursive = TRUE, showWarnings = FALSE )

  results = file.path(results.path,'results')
  dir.create( 
    results, 
    recursive = TRUE, showWarnings = FALSE )

  return(results)
}

preprocess.genomic.data <- function(r.package.path) {

  nbins = 3

  cat('Preprocess genomic data\n')

  suppressMessages(library(survival))

  cat('Loading raw data\n')
  load(file.path(r.package.path,'data','dataset.RData'))

  measurements = prob$mRNA
  genomic.instability = colSums(abs(prob$scna)>1,na.rm=T)/nrow(prob$scna)
  genes = prob$genes
  clinical = data.frame(samples = prob$samples,type = prob$types,sex = prob$sex,age = prob$age,race = prob$race,genomic.instability = genomic.instability,time = prob$survival[,1],status = 1-prob$survival[,2])
  
  # Filter under-represented subgroups    
  strata = strata(clinical$type,clinical$sex,clinical$race)
  ranked.strata = sort(table(strata),decreasing=T)
  represented.strata = ranked.strata[ranked.strata>=100]
  represented.strata = represented.strata[!is.na(represented.strata)]
  matches = strata %in% names(represented.strata)

  matches = matches & clinical$time>0
  measurements = measurements[,matches]
  clinical = clinical[matches,]
  clinical$race = factor(clinical$race)

  measurements = qnorm.col(measurements)

  fn = nrow(measurements)
  sn = ncol(measurements)
  
  cat('Get bin mapping\n')
  start.time <- Sys.time()
  bin.map = matrix(NA,nrow=fn,ncol=sn)
  strata = strata(clinical$type,clinical$sex,clinical$race)
  types = unique(strata)
  fgi = 1
  for( fgi in 1:fn ) {
    type = types[1]
    for( type in types ) {
      type.indices = strata == type
      c = as.double(measurements[fgi,type.indices])
      bin.map[fgi,type.indices] = get.bin.intervals(c,nbins)
    }
    cat(as.numeric(Sys.time() - start.time),fgi,
        paste(sprintf('%.f',100*fgi/fn),'%',sep=''))
    cat('\r')
  }

  cat(as.numeric(Sys.time() - start.time),'\n')

  save(bin.map,measurements,genes,clinical,file = file.path(r.package.path,'data','data.mRNA.RData'))
}

qnorm.col <- function(gene.exp.resid)
{
  gene.exp.norm <- gene.exp.resid
  for( sl in seq(ncol(gene.exp.resid ))) {
    mat = gene.exp.resid[,sl]
    is.na = is.na(mat)
    mat = rank(mat, ties.method = "average")
    mat[is.na] = NA
    mat = qnorm(mat / (length(mat)+1))
    gene.exp.norm[,sl] = mat
  }
 gene.exp.norm
}

get.bin.intervals <- function(c,nbins,recursion = 1) {
  
  binp = seq(from = 0, to = 1, by = 1/nbins)
  b   = unique(quantile(c,binp,na.rm=TRUE))
  if(length(b[which(!is.na(b))])==0) return(NA)
  add = 0
  if( b[1] == 0 ) add = 1
  b = b[b!=0]
  if( add == 1 ) b = c(0,b)
  if( length(b) < nbins-1 ) return(NA)
  bin = b[1:(length(b)-1)]
  
  intervals = findInterval(c, bin)-1
  counts    = table(intervals)

  len = length(c)
  if( max(abs(counts-(len/nbins))) > max(1,0.005*len)
          && recursion > 0 ) {
    noise = sd(c)/100
    c = c + runif(len, -noise, noise)
    return(get.bin.intervals(c,nbins,recursion-1))
  }

  return(intervals)
}

get.bin.quantile.intervals <- function(c,quantiles,recursion = 1) {
  
  binp = c(0,quantiles,1)
  b   = unique(quantile(c,binp,na.rm=TRUE))
  if(length(b[which(!is.na(b))])==0) return(NA)
  add = 0
  if( b[1] == 0 ) add = 1
  b = b[b!=0]
  if( add == 1 ) b = c(0,b)
  if( length(b) < nbins-1 ) return(NA)
  bin = b[1:(length(b)-1)]
  
  intervals = findInterval(c, bin)-1
  counts    = table(intervals)

  len = length(c)
  if( max(abs(counts-(len/nbins))) > max(1,0.005*len)
          && recursion > 0 ) {
    noise = sd(c)/100
    c = c + runif(len, -noise, noise)
    return(get.bin.intervals(c,nbins,recursion-1))
  }

  return(intervals)
}

run.distributed.pairwise.significance <- function( r.package.path,results.path,queues,num.jobs,memory,walltime ) {

  clear.previous.results(file.path(results.path,'results'))

  setwd(results.path)
  source(file.path(r.package.path,'R','analyze.pairwise.significance.R'))
  library(rslurm,quietly=T)
  params = data.frame(r.package.path,results.path,workers=num.jobs,id=seq(num.jobs));head(params)
  if (is.na(queues)) {
    sopt <- list(mem = memory, time = walltime, share = TRUE)
  } else {
    sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  }
  setwd(results.path);getwd()
  job.ids <- slurm_apply(analyze.pairwise.significance, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  source(file.path(r.package.path,'R','main.script.functions.R'))
  setwd(file.path(results.path))
  cat('\nRun calculate.null.molecular\n')
  return(job.ids)
}

merge.clinical.results <- function( r.package.path,results.path,p.val.quantile.threshold,large.queues,large.memory,large.walltime ) {

  start.time <- Sys.time()

  output.path = file.path(results.path,'output')

  setwd(results.path)
  source(file.path(r.package.path,'R','merge.pancancer.results.R'))
  library(rslurm,quietly=T)
  params = data.frame(r.package.path,results.path,num.jobs,p.val.quantile.threshold);head(params)
  if (is.na(large.queues)) {
    sopt <- list(mem = large.memory, time = large.walltime, share = TRUE)
  } else {
    sopt <- list(qos = large.queues, mem = large.memory, time = large.walltime, share = TRUE)
  }
  setwd(results.path);getwd()
  job.ids <- slurm_apply(merge.pancancer.results, params, nodes = 1, cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  source(file.path(r.package.path,'R','main.script.functions.R'))
  setwd(file.path(results.path))
  
  cat('\nMerging pancancer results\nRun get_slurm_output\n')

  return(job.ids)
}

calculate.null.molecular <- function (
  r.package.path,results.path,shuffle.frac = 0.1) {

  null.molecular.file = file.path(results.path,'null.molecular.RData')
  if( file.exists(null.molecular.file) ) {
    return('NULL Molecular already exists')
  }

  cat('\nCalculating Null Molecular\n')

  load(file = file.path(r.package.path,'data','data.mRNA.RData'))

  nshuffled = ceiling(shuffle.frac * choose(nrow(bin.map),2))

  suppressMessages(library(foreach))
  suppressMessages(library(doMC))
  registerDoMC(10)

  symmetric.bins = c(1:3,5:6,9)
  pair.sample = cbind(sample(nrow(bin.map),nshuffled,replace = T),sample(nrow(bin.map),nshuffled,replace = T))
  pair.sample = pair.sample[pair.sample[,1]!=pair.sample[,2],]
  nrows = nrow(pair.sample)

  stats <- foreach (row=seq(nrows),.combine=rbind,.inorder=T) %dopar% {
    x = bin.map[pair.sample[row,1],]
    y = bin.map[pair.sample[row,2],]
    bins = x*3 + y + 1
    obs = table(bins)/length(bins)
    obs = obs[symmetric.bins]
    obs
  }
  quantiles = apply(stats,2,function(v) quantile(v,c(0.05,0.95),na.rm=T))
  save(quantiles,stats,file = null.molecular.file)
  cat('\nFinished calculate.null.molecular\n')
}

calculate.base.cox.model <- function( r.package.path,results.path,queues,num.jobs,memory,walltime ) {

  output.path = file.path(results.path,'output')
  candidates.path = file.path(results.path,'candidates')
  dir.create( candidates.path,recursive = TRUE,showWarnings = FALSE )

  file.name = 'mRNA.candidates.base.cox.significance*'

  file.remove(list.files(
                path = candidates.path,
                pattern=file.name,
                include.dirs = FALSE,
                recursive = TRUE,
                full.names = TRUE))

  start.time <- Sys.time()

  setwd(results.path)
  source(file.path(r.package.path,'R','calculate.base.cox.model.R'))
  library(rslurm,quietly=T)
  params = data.frame(r.package.path,results.path,workers=num.jobs,id=seq(num.jobs));head(params)
  if (is.na(queues)) {
    sopt <- list(mem = memory, time = walltime, share = TRUE)
  } else {
    sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  }
  setwd(results.path);getwd()
  job.ids <- slurm_apply(calculate.base.cox.model, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  source(file.path(r.package.path,'R','main.script.functions.R'))
  setwd(file.path(results.path))

  cat('\nCalculating base cox model\nRun get_slurm_output\n')

  return(job.ids)
}

calculate.candidates.cox.fdr <- function( r.package.path,results.path,queues,num.jobs,memory,walltime ) {

  output.path = file.path(results.path,'output')
  candidates.path = file.path(results.path,'candidates')
  
  dir.create( candidates.path,recursive = TRUE,showWarnings = FALSE )

  file.name = 'mRNA.candidates.signed.delta.loglik*'

  file.remove(list.files(
                path = candidates.path,
                pattern=file.name,
                include.dirs = FALSE,
                recursive = TRUE,
                full.names = TRUE))

  setwd(results.path)
  source(file.path(r.package.path,'R','calculate.candidates.cox.fdr.R'))
  library(rslurm,quietly=T)
  params = data.frame(r.package.path,results.path,workers=num.jobs,id=seq(num.jobs));head(params)
  if (is.na(queues)) {
    sopt <- list(mem = memory, time = walltime, share = TRUE)
  } else {
    sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  }
  setwd(results.path);getwd()
  job.ids <- slurm_apply(calculate.candidates.cox.fdr, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  source(file.path(r.package.path,'R','main.script.functions.R'))
  setwd(file.path(results.path))

  cat('\nCalculating candidates cox fdr\nRun get_slurm_output\n')

  return(job.ids)
}


get.final.SPAGEs <- function( r.package.path,results.path,FDR.threshold = 0.99,apply.PPI.filter = F ) {
  
  cat('______________________________________________________________\n')
  start.time <- Sys.time()
  load(file = file.path(r.package.path,'data','data.mRNA.RData'))

  selected.functional.states = {}
  bin = 1
  bins = c(1:3,5:6,9)
  for( bin in bins ) {
    tryCatch({
      load(file.path(results.path,'candidates',paste('mRNA.signed.delta.loglik',bin,'pancancer.results.RData',sep='.')))
      significance.threshold = quantile(abs(shuffled.candidates.signed.delta.loglik[,3]),FDR.threshold,na.rm=T)
      significant = abs(candidates.signed.delta.loglik[,3])>significance.threshold
      selected.functional.states = rbind(selected.functional.states,cbind(candidates.signed.delta.loglik[significant,1:2],bin,candidates.signed.delta.loglik[significant,3]))
    },error = function(e) e )
  }
  selected.functional.states = selected.functional.states[rowSums(is.na(selected.functional.states))==0,]

  if( apply.PPI.filter ) {
    ppi.distance.file = file.path(results.path,paste(paste(min(bins),max(bins),sep='-'),'mRNA.HPRD.distances.RData',sep='.'))
    if( !file.exists(ppi.distance.file) ) {
      cat('Calculating PPI distances\n')
      pairs = cbind(genes[selected.functional.states[,1]],genes[selected.functional.states[,2]])
      library(data.table)
      ppi = fread(file.path(r.package.path,'data','hippie_current.txt'))
      ppi = cbind(sapply(ppi$V1,function(g) strsplit(g,'_')[[1]][1]),sapply(ppi$V3,function(g) strsplit(g,'_')[[1]][1]))
      rownames(ppi) = NULL
      ppi = ppi[ppi[,1]!=ppi[,2],]
      ppi = ppi[!is.na(ppi[,1])&!is.na(ppi[,2]),]
      ppi = unique(ppi)
      head(ppi)
      dim(ppi)
      length(unique(c(ppi[,1],ppi[,2])))
      suppressMessages(library(igraph))
      pi.network = graph.data.frame(ppi, directed=F)
      start.time <- Sys.time()
      ppi.distances = calculate.ppi.distances(pi.network,pairs)
      pair.ppi.distances = data.frame(selected.functional.states[,1:2],ppi.distances)
      print(Sys.time() - start.time)
      summary(ppi.distances)
      save(pair.ppi.distances,file=ppi.distance.file)
    } else {
      cat('Loading distances\n')
      load(ppi.distance.file)
    }

    if(!identical(pair.ppi.distances[,1],selected.functional.states[,1]) | !identical(pair.ppi.distances[,1],selected.functional.states[,1])) {
      gi.pairs = apply(selected.functional.states,1,function(v) paste(v[1:2],collapse='.'))
      ppi.pairs = apply(pair.ppi.distances,1,function(v) paste(v[1:2],collapse='.'))
      match = match( gi.pairs,ppi.pairs )
      if( sum(is.na(match))==0 ) {
        ppi.pairs = ppi.pairs[match]
        pair.ppi.distances = pair.ppi.distances[match,]
      } else {
        cat('PPI distances matrix doesn\'t match\n')
        return()
      }
    }
    sum(pair.ppi.distances[,3]==1,na.rm=T)
    sum(pair.ppi.distances[,3]==2,na.rm=T)
    cat('Removing PPI distances > 2 pairs\n')
    pair.ppi.distances[is.na(pair.ppi.distances)] = Inf
    selected.functional.states = selected.functional.states[pair.ppi.distances[,3]<3,]
    dim(selected.functional.states)
  }
  
  selected.functional.states = as.data.frame(selected.functional.states)
  selected.functional.states[,1] = genes[selected.functional.states[,1]]
  selected.functional.states[,2] = genes[selected.functional.states[,2]]
  selected.functional.states[,5] = selected.functional.states[,4] * sign(selected.functional.states[,4])
  selected.functional.states[,4] = -sign(selected.functional.states[,4])
  colnames(selected.functional.states)[4:5] = c('direction','effect')

  cat('States number =',nrow(selected.functional.states),'\n')

  print(Sys.time() - start.time)
  return(selected.functional.states)
}

calculate.ppi.distances <- function(
  pi.network,pair.names ) {
  column.1 = V(pi.network)[get.vertex.attribute(pi.network, "name", index=V(pi.network)) %in% pair.names[,1]]
  column.2 = V(pi.network)[get.vertex.attribute(pi.network, "name", index=V(pi.network)) %in% pair.names[,2]]
  map1 = match(pair.names[,1],names(column.1))
  map2 = match(pair.names[,2],names(column.2))
  distances = array(NA,nrow(pair.names))
  suppressMessages(require(doMC))
  suppressMessages(require(foreach))
  registerDoMC(cores = 12)
  d = foreach(p = seq(nrow(pair.names)), .inorder=F) %do%{
    pair = pair.names[p,]
    if( is.na(map1[p]) | is.na(map2[p]) ) return(NA)
    shortest.paths(pi.network, v=column.1[map1[p]], to=column.2[map2[p]])
  }
  distances = unlist(d)
  return(distances)
}

get_slurm_output <- function (sjob) {

  Sys.sleep(1)
  queued = length(system(paste('squeue -hn', sjob$jobname),intern = T)) > 0
  while(length(system(paste('squeue -hn', sjob$jobname),intern = T)) > 0) {
    Sys.sleep(1)
  }

  res_files <- paste0("results_", 0:(sjob$nodes - 1), ".RDS")
  tmpdir <- file.path(sjob$path,paste0("_rslurm_", sjob$jobname))
  missing_files <- setdiff(res_files, dir(path = tmpdir))
  
  if (length(missing_files) > 0) {
      missing_list <- paste(missing_files, collapse = ", ")
      warning(cat("The following files are missing:", missing_list,'\nCheck failed jobs error outputs\n'))
  }
  
  res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
  if (length(res_files) == 0) return(cat('Distributed job failed\n'))
  
  slurm_out <- lapply(res_files, readRDS)
  slurm_out <- do.call(c, slurm_out)
  slurm_out <- as.data.frame(do.call(rbind, slurm_out))

  cat('Distributed job finished\nRun cleanup_files and the next command\n')
}
