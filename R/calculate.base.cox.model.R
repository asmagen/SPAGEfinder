calculate.base.cox.model <- function(scripts,data,temp,dataset,datatype,workers,id) {

  print('calculate.base.cox.model')

  cat(date(),'\n')
  
  print(scripts)
  print(data)
  print(temp)
  print(dataset)
  print(datatype)
  print(workers)
  print(id)

  perform.cox <- function( x,y,bins,clinical,molecular.quantiles ) {

    bins.binary = sapply(1:9, function(tt) ifelse(bins == tt, 1, 0))
    colnames(bins.binary) = sapply(1:9,function(v) paste('bin',v,sep=''))
    loglik = array(F,9)
    options(warn=2)
    map = c(1:3,2,4:5,3,5:6)
    px = table(x)/length(x)
    py = table(y)/length(y)
    exp = as.vector(px%*%t(py))
    bin = 1
    for( bin in 1:9 ) {
      bin.binary = bins.binary[,bin]
      dt1 = cbind(clinical,bin=bin.binary)
      if( nrow(dt1)==0 ) next
      tryCatch({
        cox.base = coxph(Surv(time,status)~bin+strata(type,sex,race), data=dt1)
        cox.p.val = summary(cox.base)$coefficients['bin',5]
        direction = -sign(summary(cox.base)$coefficients['bin',1])
        if( cox.p.val > 0.05 ) next
        loglik[bin] = T
        obs.molecular = sum(bins==bin)/length(bins)
        if( direction == 1 ) {
          molecular = obs.molecular < exp[bin] & obs.molecular < molecular.quantiles[1,map[bin]]
        } else {
          molecular = obs.molecular > exp[bin] & obs.molecular > molecular.quantiles[2,map[bin]]
        }
        loglik[bin] = molecular
      },error=function(x){})
    }
    signed.delta.loglik = loglik
    options(warn=1)
    return(signed.delta.loglik)
  }

  get.shuffled.gene.pair.list <- function(pairs,npairs=NA,overall.start.time) {
      if( is.na(npairs) ) npairs = nrow(pairs)
      if( length(pairs)==0 | npairs==0 | length(unique(pairs[,2])) == 1 | length(unique(pairs[,2])) == 1 ) {
        cat('Not enough unique values\n')
        return()
      }

      possible.values = list()
      unq1 = unique(pairs[,1])
      unq2 = unique(pairs[,2])
      for( v in unq1 ) {
        possible.values[[paste(v)]] = NA
      }

      start.time <- Sys.time()
      shuffled.pairs = array(0,c(npairs,2))
      index = 0
      chunk = ceiling(npairs/10)
      while( index < npairs ) {
        if( (index-1) %% chunk == 0 || index == npairs) cat(round(index/npairs, 1),'')
        if( difftime(Sys.time(),overall.start.time,units='hours') > 7) break
        ind1 = pairs[sample(npairs,1),1]
        vals = possible.values[[paste(ind1)]]
        if( length(vals)==1 && is.na(vals) ) {
          vals = setdiff( unq2, pairs[pairs[,1]==ind1,2] )#unq2[unq2>ind1]
        }
        len = length(vals)
        if( len == 0 ) next
        selected = vals[sample(len,1)]
        pair = c(ind1,selected)
        possible.values[[paste(ind1)]] = setdiff(vals,selected)
        index = index + 1
        shuffled.pairs[index,] = pair
      }
      # print(Sys.time() - start.time)

      return(shuffled.pairs)
  }


  results.path = file.path(temp,dataset)

  cat('Loading data\n')
  load(file = file.path(scripts,'data',paste('data',datatype,'RData',sep='.')))

  cat('Loading candidates\n')
  load(file.path(results.path,paste(datatype,'candidates.pancancer.results.RData',sep='.')))
  molecular.stats = NA

  library(data.table)
  load(file.path(temp,dataset,paste(datatype,'null.molecular.RData',sep='.')))

  molecular.quantiles = apply(stats,2,function(v) quantile(v,c(0.45,0.55),na.rm=T))
  molecular.quantiles
  rm(stats,quantiles)

  fn = nrow(candidates)

  cat('Calculating work indices\n')
  work = ceiling(fn/workers)
  start = (id-1)*work+1
  end = min(id*work,fn)
  nrows = end - start + 1
  rm(fn,work)

  suppressMessages(library(survival))
  suppressMessages(library(survcomp))
  table(clinical$status)
  clinical$age = scale(clinical$age,center = T,scale = T)
  head(clinical)

  cat('Pairwise Cox p.vals\n')
  candidates.base.cox.significance = array(NA,c(nrows,9))
  start.time <- Sys.time()
  index = 1
  chunk = ceiling(nrows/10)
  row.index=1
  for(row.index in start:end) {
    if( (row.index-start) %% chunk == 0 || row.index-start == nrows) cat(round((row.index-start)/nrows, 1),'')
    pair = candidates[row.index,]
    row = pair[1]
    col = pair[2]
    bins = bin.map[row,]*3 + bin.map[col,] + 1
    candidates.base.cox.significance[index,] = perform.cox(
                                    x = bin.map[row,],
                                    y = bin.map[col,],
                                    bins,
                                    clinical,
                                    molecular.quantiles)
    index = index + 1
  }
  candidates.base.cox.significance = cbind(candidates[start:end,],candidates.base.cox.significance)
  print(Sys.time() - start.time)

  results.dir = file.path(results.path,'candidates')
  cat('Creating results.dir',results.dir,'\n')
  dir.create(results.dir,recursive=T)

  o.start.time <- Sys.time()
  if( id > 1 ) {
    prev.candidates.base.cox.significance.res = file.path(results.dir,paste(datatype,'candidates.base.cox.significance',id-1,'pancancer.results.RData',sep='.'))
    if(!file.exists(file = prev.candidates.base.cox.significance.res)) {
      cat('Waiting for',id-1,'to finish\n')
      while( !file.exists(prev.candidates.base.cox.significance.res) ){ Sys.sleep(1) }
    }
    cat('Merging up to',id,'\n')
    candidates.base.cox.significance.last = candidates.base.cox.significance
    original.id = id
    rm(candidates.base.cox.significance)
    while( !exists('candidates.base.cox.significance') ) {
      tryCatch({
        load(prev.candidates.base.cox.significance.res)
      }, warning = function(w) {}, error = function(e) {})
      Sys.sleep(1)
    }
    id = original.id
    candidates.base.cox.significance = rbind(candidates.base.cox.significance,candidates.base.cox.significance.last)

    if( id == workers ){

      candidates.base.cox.significance[is.na(candidates.base.cox.significance)] = 0
      candidates.base.cox.significance = candidates.base.cox.significance[rowSums(candidates.base.cox.significance[,3:11])>0,]

      # symmetric correction
      base.bins = c(1,2:3,5:6,9) + 2
      symmetric = c(4,7,8) + 2
      map = c(2,3,6) + 2
      row.sums = rowSums(candidates.base.cox.significance[,symmetric])
      new.rows = candidates.base.cox.significance[row.sums > 0,]
      candidates.base.cox.significance[row.sums > 0,symmetric] = 0
      new.rows[,base.bins] = 0
      new.rows[,1:2] = new.rows[,2:1]
      new.rows[,map] = new.rows[,symmetric]
      new.rows[,symmetric] = 0
      candidates.base.cox.significance = rbind(candidates.base.cox.significance,new.rows)
      candidates.base.cox.significance = candidates.base.cox.significance[,c(1:2,base.bins)]
      colnames(candidates.base.cox.significance) = c('y','x',base.bins-2)

      cat('Calculating shuffled candidates\n')
      shuffled.candidates = list()
      for( bin in 3:ncol(candidates.base.cox.significance) ) {
        shuffled.candidates[[base.bins[bin-2]-2]] = get.shuffled.gene.pair.list(candidates.base.cox.significance[candidates.base.cox.significance[,bin]==1,1:2],overall.start.time=Sys.time())
      }

      cat('Saving merged file',id,'\n')
      candidates.base.cox.significance.res = file.path(results.path,paste(datatype,'candidates.base.cox.significance.pancancer.results.RData',sep='.'))
      save(candidates.base.cox.significance,shuffled.candidates,file = candidates.base.cox.significance.res)
    }
  }
  cat('Bin Base Significant Column Sums\n');print(colSums(candidates.base.cox.significance[,3:ncol(candidates.base.cox.significance)]))
  print(head(candidates.base.cox.significance,20))
  cat('Saving',id,'\n')
  candidates.base.cox.significance.res = file.path(results.dir,paste(datatype,'candidates.base.cox.significance',id,'pancancer.results.RData',sep='.'))
  save(id,candidates.base.cox.significance,file = candidates.base.cox.significance.res)
  print(Sys.time() - o.start.time)

  cat('End\n')
}
