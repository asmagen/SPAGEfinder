
calculate.candidates.cox.fdr <- function(r.package.path,results.path,workers,id) {

  print('calculate.candidates.cox.fdr')

  cat(date(),'\n')
  
  print(r.package.path)
  print(results.path)
  print(workers)
  print(id)

  perform.cox <- function( mRNA1,mRNA2,bins,clinical ) {

    options(warn=2)
    bins.numeric = bins
    bins.numeric[bins.numeric==T] = 1
    dt2 = cbind(bin=bins.numeric,gene1=mRNA1,gene2=mRNA2,clinical)
    signed.delta.loglik = NA
    effect.size = NA
    tryCatch({
      cox.out = coxph(Surv(time,status)~gene1+gene2+bin+age+strata(sex,race,type), data=dt2)
      coefficients = summary(cox.out)$coefficients
      signs = sign(coefficients['bin',1]) * -1
      signed.delta.loglik = signs * -log(coefficients['bin',5])
      effect.size = coefficients['bin',1] * -1
    },error=function(x){})
    options(warn=1)
    return(list(signed.delta.loglik=signed.delta.loglik,effect.size=effect.size))
  }

  cat('Loading candidates\n')
  load(file.path(results.path,'mRNA.candidates.base.cox.significance.pancancer.results.RData'))
  head(candidates.base.cox.significance)
  nrow(candidates.base.cox.significance)
  str(shuffled.candidates)
  
  nbins = sum(sapply(shuffled.candidates,length)>0)

  cat('Calculating work indices\n')
  ids = seq(workers)
  workers = floor(workers/nbins)
  ids = ids[seq(workers * nbins)]
  bins = ids %% nbins
  bins[bins==0] = nbins
  if( id > workers * nbins ) return()
  bin = bins[id]
  map = c(1:3,5:6,9)
  mapped.bin = map[bins[id]]
  id = which(which(bins==bins[id])==id)
  candidates = candidates.base.cox.significance[candidates.base.cox.significance[,bin+2]==1,1:2]
  shuffled.candidates = shuffled.candidates[[mapped.bin]]
  fn = nrow(candidates)
  work = floor(fn/workers)
  start = (id-1)*work+1
  end = min(id*work,fn)
  if( id == workers ) end = fn
  nrows = end - start + 1

  if (start >= end) return(NA)

  cat('Loading data for cox\n')

  load(file = file.path(r.package.path,'data','data.mRNA.RData'))

  suppressMessages(library(survival))
  suppressMessages(library(survcomp))
  table(clinical$status)
  clinical$age = scale(clinical$age,center = T,scale = T)

  cat('Pairwise Cox p.vals\n')
  candidates.signed.delta.loglik = array(NA,nrows)
  effect.size = array(NA,nrows)
  shuffled.candidates.signed.delta.loglik = array(NA,nrows)
  shuffled.effect.size = array(NA,nrows)
  start.time <- Sys.time()
  index = 1
  chunk = ceiling(nrows/10)
  row.index = 1
  for(row.index in start:end) {
    if( (row.index-start) %% chunk == 0 || row.index-start == nrows) cat(round((row.index-start)/nrows, 1),'')
    pair = candidates[row.index,]
    row = pair[1]
    col = pair[2]
    bins = bin.map[row,]*3 + bin.map[col,] + 1 == mapped.bin
    res = perform.cox(
                      mRNA1 = scale(measurements[row,],center = T,scale = T),
                      mRNA2 = scale(measurements[col,],center = T,scale = T),
                      bins,
                      clinical)
    candidates.signed.delta.loglik[index] = res$signed.delta.loglik
    effect.size[index] = res$effect.size
    pair = shuffled.candidates[row.index,]
    row = pair[1]
    col = pair[2]
    bins = bin.map[row,]*3 + bin.map[col,] + 1 == mapped.bin
    res = perform.cox(
                      mRNA1 = measurements[row,],
                      mRNA2 = measurements[col,],
                      bins,
                      clinical)
    shuffled.candidates.signed.delta.loglik[index] = res$signed.delta.loglik
    shuffled.effect.size[index] = res$effect.size
    index = index + 1
  }
  candidates.signed.delta.loglik = cbind(candidates[start:end,],candidates.signed.delta.loglik,effect.size)
  shuffled.candidates.signed.delta.loglik = cbind(shuffled.candidates[start:end,],shuffled.candidates.signed.delta.loglik,shuffled.effect.size)
  print(Sys.time() - start.time)

  results.dir = file.path(results.path,'candidates')
  cat('Creating results.dir',results.dir,'\n')
  dir.create(results.dir,recursive=T)

  o.start.time <- Sys.time()
  if( id > 1 ) {
    
    prev.candidates.signed.delta.loglik.res = file.path(results.dir,paste('mRNA.signed.delta.loglik',mapped.bin,id-1,'pancancer.results.RData',sep='.'))
    
    if(!file.exists(file = prev.candidates.signed.delta.loglik.res)) {
      cat('Waiting for',id-1,'to finish\n')
      while( !file.exists(prev.candidates.signed.delta.loglik.res) ){ Sys.sleep(1)}
    }
    cat('Merging up to',id,'\n')
    candidates.signed.delta.loglik.last = candidates.signed.delta.loglik
    shuffled.candidates.signed.delta.loglik.last = shuffled.candidates.signed.delta.loglik
    effect.size.last = effect.size
    shuffled.effect.size.last = shuffled.effect.size
    original.id = id
    rm(candidates.signed.delta.loglik)
    rm(shuffled.candidates.signed.delta.loglik)
    while( !exists('shuffled.candidates.signed.delta.loglik') ) {
      tryCatch({
        load(prev.candidates.signed.delta.loglik.res)
      }, warning = function(w) {}, error = function(e) {})
      Sys.sleep(1)
    }
    id = original.id
    candidates.signed.delta.loglik = rbind(candidates.signed.delta.loglik,candidates.signed.delta.loglik.last)
    shuffled.candidates.signed.delta.loglik = rbind(shuffled.candidates.signed.delta.loglik,shuffled.candidates.signed.delta.loglik.last)
    effect.size = rbind(effect.size,effect.size.last)
    shuffled.effect.size = rbind(shuffled.effect.size,shuffled.effect.size.last)

    if( id == workers ){
      cat('Saving merged file',id,'\n')
      candidates.signed.delta.loglik.res = file.path(results.dir,paste('mRNA.signed.delta.loglik',mapped.bin,'pancancer.results.RData',sep='.'))
      
      save(candidates.signed.delta.loglik,shuffled.candidates.signed.delta.loglik,effect.size,shuffled.effect.size,file = candidates.signed.delta.loglik.res)
    }
  }
  
  cat('Saving',id,'\n')
  candidates.signed.delta.loglik.res = file.path(results.dir,paste('mRNA.signed.delta.loglik',mapped.bin,id,'pancancer.results.RData',sep='.'))
  
  save(id,bin,candidates.signed.delta.loglik,shuffled.candidates.signed.delta.loglik,file = candidates.signed.delta.loglik.res)
  print(Sys.time() - o.start.time)

  cat('End\n')
}
