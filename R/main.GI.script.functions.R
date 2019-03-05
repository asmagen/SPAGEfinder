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

  dir.create( 
    file.path(results.path,'output'), 
    recursive = TRUE, showWarnings = FALSE )

  results = file.path(results.path,'results')
  dir.create( 
    results, 
    recursive = TRUE, showWarnings = FALSE )

  dir.create( 
    file.path(results.path,'plots'), 
    recursive = TRUE, showWarnings = FALSE )

  return(results)
}

preprocess.genomic.data <- function(
				data,temp,dataset,results.path,constrain.gene.set = NA) {

  nbins = 3

  cat('Preprocess genomic data\n')
  file.name = file.path(scripts,'data','data.mRNA.RData')

  suppressMessages(library(survival))
  cat('Loading raw data\n')

  load(file.path(data,'prob.TCGA.extended.RData'))

  measurements = prob$mRNA
  dim(measurements)
  summary(measurements[,1:10])

  genomic.instability = colSums(abs(prob$scna)>1,na.rm=T)/nrow(prob$scna)

  genes = prob$genes
  clinical = data.frame(samples = prob$samples,type = prob$types,sex = prob$sex,age = prob$age,race = prob$race,stage = prob$stage,genomic.instability = genomic.instability,time = prob$survival[,1],status = 1-prob$survival[,2])
  head(clinical)
  
  # Filter under-represented subgroups    
  strata = strata(clinical$type,clinical$sex,clinical$race)
  ranked.strata = sort(table(strata),decreasing=T)
  represented.strata = ranked.strata[ranked.strata>=100]
  represented.strata = represented.strata[!is.na(represented.strata)]
  write.csv(represented.strata,file=file.path(temp,'confounder.strata.csv'))
  names(represented.strata)
  matches = strata %in% names(represented.strata)
  sum(matches)
  apply(sapply(names(represented.strata),function(strat){ strsplit(strat,', ')[[1]] }),1,unique)

  matches = matches & clinical$time>0
  measurements = measurements[,matches]
  clinical = clinical[matches,]
  clinical$race = factor(clinical$race)
  summary(clinical)
  head(clinical)
  dim(clinical)
  dim(measurements)

  if (length(constrain.gene.set) > 1) {
    indices = genes %in% constrain.gene.set
    measurements = measurements[indices,]
    genes = genes[indices]
  }

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

  save(bin.map,measurements,genes,clinical,file = file.name)
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

run.distributed.pairwise.significance <- function( scripts,temp,dataset,results.path,queues,num.jobs,memory,walltime ) {

  output.path = file.path(results.path,'output')

  setwd(results.path)
  source(file.path(scripts,'analyze.pairwise.significance.R'))
  library(rslurm,quietly=T)
  params = data.frame(scripts,temp,dataset,'mRNA',workers=num.jobs,id=seq(num.jobs));head(params)
  sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(analyze.pairwise.significance, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  setwd(scripts);source('main.GI.script.functions.R')
  return(job.ids)
}

merge.clinical.results <- function( temp,scripts,dataset,p.val.quantile.threshold,base.res.path,results.path,large.queues,memory,walltime ) {

  start.time <- Sys.time()

  output.path = file.path(base.res.path,'output')

  cat('Merging pancancer results\n')

  setwd(results.path)
  source(file.path(scripts,'merge.pancancer.results.R'))
  library(rslurm,quietly=T)
  params = data.frame(scripts,temp,dataset,'mRNA',num.jobs,p.val.quantile.threshold);head(params)
  sopt <- list(qos = large.queues, mem = memory, time = walltime, share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(merge.pancancer.results, params, nodes = 1, cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  setwd(scripts);source('main.GI.script.functions.R')
  return(job.ids)

  print(Sys.time() - start.time)
}

calculate.null.molecular <- function (
  temp,dataset,null.molecular.file,shuffle.frac = 0.1) {

  load(file = file.path(scripts,'data',dataset,'data.mRNA.RData'))

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
}

calculate.base.cox.model <- function( scripts,data,temp,dataset,results.path,queues,num.jobs,memory,walltime ) {

  cat('Calculating base cox model\n')

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
  source(file.path(scripts,'calculate.base.cox.model.R'))
  library(rslurm,quietly=T)
  params = data.frame(scripts,data,temp,dataset,'mRNA',workers=num.jobs,id=seq(num.jobs));head(params)
  sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(calculate.base.cox.model, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  setwd(scripts);source('main.GI.script.functions.R')
  return(job.ids)
}

calculate.candidates.cox.fdr <- function( scripts,data,temp,dataset,results.path,queues,num.jobs,memory,walltime ) {

  cat('Calculating candidates cox fdr\n')

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
  source(file.path(scripts,'calculate.candidates.cox.fdr.R'))
  library(rslurm,quietly=T)
  params = data.frame(scripts,data,temp,dataset,'mRNA',support,workers=num.jobs,id=seq(num.jobs));head(params)
  sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(calculate.candidates.cox.fdr, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt)
  job.ids$path = getwd()
  setwd(scripts);source('main.GI.script.functions.R')
  return(job.ids)
}


get.final.GIs <- function( temp,data,results.path,LLR.threshold,PPI = F ) {
  
  cat('______________________________________________________________\n')
  start.time <- Sys.time()
  load(file = file.path(scripts,'data','data.mRNA.RData'))

  selected.functional.states = {}
  bin = 1
  for( bin in bins ) {
    tryCatch({
      load(file.path(results.path,paste('mRNA.candidates.signed.delta.loglik',bin,'pancancer.results.RData',sep='.')))
      significance.threshold = quantile(abs(shuffled.candidates.signed.delta.loglik[,3]),LLR.threshold,na.rm=T)
      significant = abs(candidates.signed.delta.loglik[,3])>significance.threshold
      selected.functional.states = rbind(selected.functional.states,cbind(candidates.signed.delta.loglik[significant,1:2],bin,candidates.signed.delta.loglik[significant,3]))
    },error = function(e) e )
  }
  selected.functional.states = selected.functional.states[rowSums(is.na(selected.functional.states))==0,]

  if( PPI ) {
    ppi.distance.file = file.path(results.path,paste(paste(min(bins),max(bins),sep='-'),interaction.direction,'mRNA.HPRD.distances.RData',sep='.'))
    if( !file.exists(ppi.distance.file) ) {
      cat('Calculating PPI distances\n')
      pairs = cbind(genes[selected.functional.states[,1]],genes[selected.functional.states[,2]])
      library(data.table)
      ppi = fread(file.path(data,'hippie_current.txt'))
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
  colnames(selected.functional.states)[4] = 'effect'

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
