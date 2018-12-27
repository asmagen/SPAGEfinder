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
				data,temp,dataset,datatype,nbins,results.path,
        sampling.pool,quantile.threshold,validation=0) {

  cat('Preprocess genomic data\n')
  file.name = file.path(results.path,
                paste(
                  ifelse(validation==0,'data','validation.data'),
                  datatype,
                  'RData',
                sep='.')
                )

  suppressMessages(library(survival))
  if( !file.exists(file = file.name) ) {
    
    cat('Loading raw data\n')
    if(validation) {
      
      load(file.path(temp,'prob.curtis.extended.RData'))
      if( datatype=='mRNA' )
      {
        measurements = prob$mRNA
      } else if( datatype=='scna' ) {
        measurements = prob$scna
      }
      genomic.instability = colSums(abs(prob$scna)>1,na.rm=T)/nrow(prob$scna)
      head(genomic.instability)
      genes = prob$genes
      clinical = data.frame(samples = prob$samples,type = 'BRCA',sex = 'FEMALE',age = prob$age,race = 'WHITE',stage = prob$stage,grade = prob$grade,subtypes = prob$subtypes,size = prob$size,p53.mut = prob$p53.mut,cellularity = prob$cellularity,treatment = prob$treatment,histological.type = prob$histological.type,er = prob$er,her2 = prob$her2,pr = prob$pr,genomic.instability=genomic.instability,time = prob$survival[,1],status = 1-prob$survival[,2])
    } else {

      load(file.path(temp,'prob.TCGA.extended.RData'))

      if( datatype=='mRNA' )
      {
        measurements = prob$mRNA
      } else if( datatype=='scna' ) {
        measurements = prob$scna
      }
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
      # clinical$status = factor(clinical$status)
      summary(clinical)
      head(clinical)
      dim(clinical)
      dim(measurements)
    }
    
    fn = nrow(measurements)
    sn = ncol(measurements)

    measurements = qnorm.col(measurements)
    
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
    
  } else {
    load(file = file.name)
    fn = nrow(bin.map)
    sn = ncol(bin.map)
  }

  if( validation == 0 ){
      save(bin.map,measurements,genes,clinical,file = file.name)
    } else {
        save(bin.map,measurements,genes,clinical,
              file = file.path(results.path,paste(
                  'validation',
                  'data',
                  datatype,
                  'RData',
                sep='.')))
    }
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

run.distributed.pairwise.significance <- function(
        scripts,temp,dataset,datatype,results.path,
        queues,jobs,batches,memory,walltime) {

  output.path = file.path(results.path,'output')

  num.jobs = sum(jobs)

  source('/cbcbhomes/magen/RSR/pan.cancer/analyze.pairwise.significance.R')
  library(rslurm,quietly=T)
  params = data.frame(scripts,temp,dataset,datatype,workers=num.jobs,id=seq(num.jobs));head(params)
  sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(analyze.pairwise.significance, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt,libPaths = "/cbcb/sw/RedHat-7-x86_64/common/local/R/common/3.3.1")#, add_objects = c('path','data.path')
  return(job.ids)
}

merge.clinical.results <- function(
      temp,scripts,dataset,datatype,p.val.quantile.threshold,base.res.path,num.jobs,results.path ) {

  start.time <- Sys.time()

  output.path = file.path(base.res.path,'output')

  cat('Merging pancancer results\n')

  source('/cbcbhomes/magen/RSR/pan.cancer/merge.pancancer.results.R')
  library(rslurm,quietly=T)
  params = data.frame(scripts,temp,dataset,datatype,num.jobs,p.val.quantile.threshold);head(params)
  sopt <- list(qos = 'large', mem = '120GB', time = '3:00:00', share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(merge.pancancer.results, params, nodes = 1, cpus_per_node = 1, submit = TRUE, slurm_options = sopt,libPaths = "/cbcb/sw/RedHat-7-x86_64/common/local/R/common/3.3.1")#, add_objects = c('path','data.path')
  return(job.ids)

  print(Sys.time() - start.time)
}

calculate.null.molecular <- function (
  temp,dataset,datatype,null.molecular.file,shuffle.frac = 0.1) {

  load(file = file.path(temp,dataset,paste('data',datatype,'RData',sep='.')))

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

calculate.candidates.cox.fdr <- function(
        scripts,data,temp,dataset,datatype,p.val.quantile.threshold,support,
        results.path,queue.name,num.jobs,memory,walltime) {

  cat('Calculating candidates cox fdr\n')

  output.path = file.path(results.path,'output')
  candidates.path = file.path(results.path,'candidates')
  
  dir.create( candidates.path,recursive = TRUE,showWarnings = FALSE )

  if( support == 1 ) {
    file.name = paste(datatype,'support.candidates.signed.delta.loglik*',sep='.')
  } else {
    file.name = paste(datatype,'candidates.signed.delta.loglik*',sep='.')
  }

  file.remove(list.files(
                path = candidates.path,
                pattern=file.name,
                include.dirs = FALSE,
                recursive = TRUE,
                full.names = TRUE))

  source('/cbcbhomes/magen/RSR/pan.cancer/calculate.candidates.cox.fdr.R')
  library(rslurm,quietly=T)
  params = data.frame(scripts,data,temp,dataset,datatype,support,workers=num.jobs,id=seq(num.jobs));head(params)
  sopt <- list(qos = queues, mem = memory, time = walltime, share = TRUE)
  setwd(results.path);getwd()
  job.ids <- slurm_apply(calculate.candidates.cox.fdr, params, nodes = nrow(params), cpus_per_node = 1, submit = TRUE, slurm_options = sopt,libPaths = "/cbcb/sw/RedHat-7-x86_64/common/local/R/common/3.3.1")#, add_objects = c('path','data.path')
  return(job.ids)
}


load.train.data <- function(
                      temp, 
                      data,
                      datatype,
                      results.path,
                      LLR.threshold,
                      use.fdr = T,
                      opposite.direction=F,
                      bins = 1:9,
                      limit = NA,
                      interaction.direction = 0,
                      shuffle = F,
                      random = F,
                      include.bins = F,
                      PPI = T,
                      functional.interactions = F,
                      random.subset = NA ) {
  
  cat('______________________________________________________________\n')
  start.time <- Sys.time()
  cat('Load train data\n')

  cat('Lodaing bins & survival data\n')
  non.na.genes = NA
  tryCatch({
    load(file = file.path(results.path,paste('validation.data',datatype,'RData',sep='.')))
    non.na.genes = which(rowSums(is.na(measurements))<ncol(measurements)/2)
  },error = function(c) c)
  load(file = 
    file.path(results.path,paste('data',datatype,'RData',sep='.')))
  load(file.path(temp,'prob.TCGA.extended.RData'))
  genomic.instability = NA
  genes = prob$genes

    selected.functional.states = {}
    significance.thresholds = {}
    bin = 1
    for( bin in bins ) {
      tryCatch({
        load(file.path(results.path,paste(datatype,'candidates.signed.delta.loglik',bin,'pancancer.results.RData',sep='.')))

        if( use.fdr ) {
          significance.threshold = quantile(abs(shuffled.candidates.signed.delta.loglik[,3]),LLR.threshold,na.rm=T)
        } else {
          significance.threshold = quantile(abs(candidates.signed.delta.loglik[,3]),LLR.threshold,na.rm=T)
        }
        significance.thresholds = c(significance.thresholds,significance.threshold)
        if( opposite.direction == T ) {
          significant = abs(candidates.signed.delta.loglik[,3])<significance.threshold
        } else {
          significant = abs(candidates.signed.delta.loglik[,3])>significance.threshold
        }
        print(sum(significant,na.rm=T))
        
        selected.functional.states = rbind(selected.functional.states,cbind(candidates.signed.delta.loglik[significant,1:2],bin,candidates.signed.delta.loglik[significant,3]))
      },error = function(e) e )
    }
    dim(selected.functional.states)
    head(selected.functional.states)

    selected.functional.states = selected.functional.states[rowSums(is.na(selected.functional.states))==0,]
    if( interaction.direction != 0 ) {
      selected.functional.states = selected.functional.states[sign(selected.functional.states[,4])==interaction.direction,]
    }

    dim(selected.functional.states)
    head(selected.functional.states)

    if( !is.na(non.na.genes) ) {
      cat('Removing pairs not existing in validation\n')
      pairs = selected.functional.states[,1:2]
      indices = pairs[,1]%in%non.na.genes & pairs[,2]%in%non.na.genes
      selected.functional.states = selected.functional.states[indices,]
    }

    cat('Removing duplicate states\n')
    df = selected.functional.states[,1:2]
    rownames(df) = seq(nrow(df))
    udf = unique(df)
    indices = as.numeric(rownames(udf))
    selected.functional.states = selected.functional.states[indices,]

  dim(selected.functional.states)

  cat('Selected bin summary')
  print(round(table(selected.functional.states[,3])/nrow(selected.functional.states),2))
  print(round(table(selected.functional.states[,3]*sign(selected.functional.states[,4]))/nrow(selected.functional.states),2))

  if( !is.na(limit) ) {
    cat('*** LIMITING TO TOP',limit,'RANKED STATES ***\n')
    selected.functional.states = selected.functional.states[order(abs(selected.functional.states[,4]),decreasing=T),][1:min(limit,nrow(selected.functional.states)),]
    cat('After limit selected bin summary')
    print(round(table(selected.functional.states[,3])/nrow(selected.functional.states),2))
  }

  if( PPI ) {
    ppi.distance.file = file.path(results.path,paste(paste(min(bins),max(bins),sep='-'),interaction.direction,datatype,'HPRD.distances.RData',sep='.'))
    if( !file.exists(ppi.distance.file) ) {
      cat('Calculating PPI distances\n')
      pairs = cbind(genes[selected.functional.states[,1]],genes[selected.functional.states[,2]])
      library(data.table)
      ppi = fread(file.path(data,'hippie_current.txt'))
      ppi = cbind(sapply(ppi$V1,function(g) strsplit(g,'_')[[1]][1]),sapply(ppi$V3,function(g) strsplit(g,'_')[[1]][1]))
      # ppi = data.frame(fread(file.path(temp,'biogrid.csv')));head(ppi)
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

  if( !is.na(random.subset) ) {
    cat('**** Get random subset of pairs ****\n')
    selected.functional.states = selected.functional.states[sample(nrow(selected.functional.states),random.subset*nrow(selected.functional.states)),] 
  }

  if( shuffle == T ) {
    cat('**** SHUFFLING STATES ****\n')
    nrows = nrow(selected.functional.states)
    if( include.bins == T ) {
      bins.row.order = sample(nrows,nrows)
    } else {
      bins.row.order = seq(nrow(selected.functional.states))
    }
    selected.functional.states[,1:3] = cbind(
                                    get.shuffled.gene.pair.list(
                                      selected.functional.states[,1:2],overall.start.time=Sys.time() ),
                                    selected.functional.states[bins.row.order,3]  )
  }

  if( random == T ) {
    cat('**** RANDOMIZING STATES ****\n')
    nrows = nrow(selected.functional.states)
    ngenes = nrow(bin.map)
    if( include.bins == T ) {
      bins = sample(9,nrows,replace=T)
    } else {
      bins = selected.functional.states[,3]
    }
    selected.functional.states[,1:3] = cbind(
      sample(ngenes,nrows,replace=T),
      sample(ngenes,nrows,replace=T),
      bins
      )
  }
  
  cat('States number =',nrow(selected.functional.states),'\n')

  s = selected.functional.states[,3]*sign(selected.functional.states[,4])
  s[abs(s)==4] = 2*sign(s[abs(s)==4])
  s[abs(s)==7] = 3*sign(s[abs(s)==7])
  s[abs(s)==8] = 6*sign(s[abs(s)==8])
  t = table(s)/length(s)

  clinical$genomic.instability=genomic.instability
  res = list(bins=bin.map,genes=genes,measurements=measurements,clinical=clinical,states=selected.functional.states,bin.abundance = t,significance.thresholds=significance.thresholds,genomic.instability)
  # str(res)
  print(Sys.time() - start.time)
  return(res)
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
