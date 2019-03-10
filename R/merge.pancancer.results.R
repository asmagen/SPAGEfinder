
merge.pancancer.results <- function(r.package.path,results.path,num.jobs,p.val.quantile.threshold) {
 
 print('merge.pancancer.results')

 cat(date(),'\n')

 overall.start.time = Sys.time()

 print(r.package.path)
 print(results.path)
 print(num.jobs)
 cat('p.val.quantile.threshold =',p.val.quantile.threshold,'\n')

  res.file.path = file.path(results.path,'results')

  candidates.file = file.path(results.path,'mRNA.candidates.pancancer.results.RData')
  quantiles.file = file.path(results.path,'mRNA.quantiles.pancancer.results.RData')
  merged = file.path(results.path,'mRNA.merged.clinical.pancancer.results.RData')
  if( !file.exists(file = candidates.file) ) {

    if( !file.exists(file = quantiles.file) ) {

     if( !file.exists(file = merged) ) {

       file.list = list.files(res.file.path,pattern = 'mRNA.clinical.*.RData')
       
       file.list = file.list[
                     order(
                       as.numeric(
                         lapply(
                           strsplit(file.list,'\\.'),
                           function(x) x[[3]]
                           )
                         )
                       )
                     ]
       num.files = length(file.list)

      cat('Merging',num.files,'clinical files\n')
      start.time <- Sys.time()
      load(file.path(res.file.path,file.list[1]))
      fn = ncol(clinical.p.vals)
      merged.clinical.p.vals = array(NA,dim=c(fn,fn,9))
      loc = 1
      count = 1
      for( file in file.list ) {
        if( count%%10 == 0 ) cat(count,'')
        count = count + 1
         tryCatch({
           load(file.path(res.file.path,file))
           len = nrow(clinical.p.vals)
           indices = loc:(loc+len-1)
           merged.clinical.p.vals[indices,,] = clinical.p.vals
           loc = loc + len
         }, warning = function(w) {
         }, error = function(e) {
          print(e)
         })
       }
        cat('\n')
        merged.clinical.p.vals[merged.clinical.p.vals==-1000]=NA
        print(Sys.time() - start.time)

         if( ncol(merged.clinical.p.vals) == loc ) {
          cat('Saving\n') 
          start.time <- Sys.time()
           save(
              merged.clinical.p.vals,
             file = merged)
           print(Sys.time() - start.time)
         } else {
           print('FALSE ncol(merged.clinical.p.vals) == loc')
           cat('ncol(merged.clinical.p.vals) =',ncol(merged.clinical.p.vals),'loc =',loc)
           return()
         }
        
      } else {
        cat('Loading previously calculated merged results\n')
        start.time <- Sys.time()
        load(merged)
        fn = ncol(merged.clinical.p.vals)
        print(Sys.time() - start.time)
      }

      cat('Computing quantiles\n')
      start.time <- Sys.time()
      max.p.val.threshold = 0.05
      quantiles = matrix(NA,fn,9)

      for( gene in 1:fn ) {
        for( bin in 1:9 ) {
          p.vals = c(merged.clinical.p.vals[1:(gene-1),gene,bin],merged.clinical.p.vals[gene,min((gene+1),fn):fn,bin])
          quantiles[gene,bin] = quantile(abs(p.vals),1-p.val.quantile.threshold,na.rm=T,names=F)
        }
      }
      cat('Saving quantiles\n')
      save(quantiles,fn,file = quantiles.file)
      print(Sys.time() - start.time)
    } else {
      cat('Loading previously calculated quantiles\n')
      load(quantiles.file)
      load(merged)
    }

    cat('Calculating candidates\n')
    start.time <- Sys.time()
    index = 1
    aggregate.index = 1
    candidates = array(0,c(ceiling(fn*fn*2*(1-0.5*p.val.quantile.threshold)),11))
    matches = 0
    chunk = ceiling(fn/10)
    for(row in seq(fn-1)) {
      if( (row-1) %% chunk == 0 || row == fn) cat(round(row/fn, 1),'')
      for(col in (row+1):fn) {
        p.vals = merged.clinical.p.vals[row,col,]
        pattern = sign(p.vals) * (abs(p.vals) <= pmin(quantiles[row,],quantiles[col,],na.rm=T))
        if( sum(abs(pattern),na.rm=T) > 0 ) {
          candidates[index,] = c(row,col,pattern)
          index = index + 1
        }
      }
    }
    cat('\n')
    candidates = candidates[candidates[,1]!=0,]
    cat('Num candidates = ',nrow(candidates),'\n')

    cat('Calculating shuffled candidates\n')
    shuffled.candidates = NA

    print(Sys.time() - start.time)
    cat('Saving candidates\n')
    start.time <- Sys.time()
    patterns = candidates[,3:11]
    candidates = candidates[,1:2]
    save(candidates,patterns,file=candidates.file)
    print(Sys.time() - start.time)

  } else {
    cat('Loading previously calculated quantiles\n')
    load(candidates.file)
    cat('nrow(candidates) =',nrow(candidates))
  }

  cat('End Merge','\n')
}
