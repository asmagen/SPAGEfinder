
analyze.pairwise.significance <- function(r.package.path,results.path,workers,id) {

  print('analyze pairwise significance')

  cat(date(),'\n')

  print(r.package.path)
  print(results.path)
  print(workers)
  print(id)

  analyse.datatype.significance <- function(
      r.package.path,
      results.path,
      workers,
      id) {

    load(file = file.path(r.package.path,'data','data.mRNA.RData'))

    cancerType = clinical$type
    typeInx         = factor(cancerType)
    typeNum         = length(unique(typeInx))
    levels(typeInx) = seq(typeNum)-1
    survival = cbind(clinical$time,1-clinical$status)
    table(survival[,2])
    head(survival)

    fn = nrow(bin.map)
    indices = calculate.work.indices( fn,workers,id )
    start = indices[1]
    end = indices[2]
    rows = end-start
    if (start >= end) return(NA)

    library(Rcpp)
    sourceCpp(file.path(r.package.path,'R',"aggregateLogRankGene.cpp"))
    
    cat('Pairwise analysis',start,end-1,'\n')
    start.time <- Sys.time()
    cat('Clinical\n')
    clinical.p.vals = array(NA,dim=c(rows,fn,9))
    if( typeNum == 1 ) typeInx[T] = 0
    index = 1
    for( index in start:(end-1) ) {
      # cat('\n',index,' ')
      row = index - start + 1
      col = (index+1):fn
      pairs = cbind(
                    rep(index, each = length(col)),
                    col
                    )
      out = aggregateLogRankPairsMultiClass(pairs=pairs, scna1Mat=bin.map, scna2Mat=bin.map, survival = survival, typeInx=as.numeric(typeInx)-1, typeNum=typeNum, threads=12,classMap=1:9,numClass=9)

      clinical.p.vals[row,col,] = out$p
    }

    cat('Save clinical\n')
    save(clinical.p.vals,start,
    file = file.path(results.path,
                      'results',
                      paste('mRNA.clinical',id,'RData',
                        sep = ".")))

    print(Sys.time() - start.time)
  }

  calculate.work.indices <- function( fn,workers,id ) {

    pairs      = choose(fn,2)
    job        = rep(floor(pairs/workers),workers)
    indices    = array(0,dim = (workers+1))
    indices[1] = 1
    work       = array(0,dim = workers)
    gap        = array(0,dim = workers)
    for( i in 1:workers ) {
      indices[i+1] = indices[i]
      while( (i==workers || work[i] < job[i]) && indices[i+1] < fn) {
        gap[i]       = gap[i] + 1
        work[i]      = work[i] + fn - indices[i+1]
        indices[i+1] = indices[i+1] + 1
      }
    }
    gs = indices[id]
    gt = indices[id+1]

    return( c(gs,gt) )
  }

    analyse.datatype.significance(
      r.package.path,
      results.path,
      workers,
      id
    )
}

