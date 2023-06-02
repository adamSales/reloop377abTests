
if(full | estMain){
  print('main estimation')
  resTotal <-
    datPW%>%
    split(~problem_set+model)%>%
    LAP(.,function(x) {
      #cat(".")
      filename <- paste0('results/miniResults/',x$problem_set[1],x$model[1],'.RData')
      if(file.exists(filename)){
        load(filename)
        return(out)
      }

      out <- try(
                                allEst(
                                  Y=x$completion_target,
                                  Tr=x$Z,
                                  Z=as.matrix(x[,covNames]),
                                  yhat=as.matrix(x[,'completion_prediction']),
                                  ps=x$problem_set[1],
                                  model=x$model[1],
                                  fast=fast
                                )
                              )
      save(out,file=filename)
      out                                
    }
            )

  save(resTotal,file=paste0('results/resTotal',S,'.RData'))
  if(!keepIntermediateResults) file.remove("results/miniResults/*")

} else load(paste0('results/resTotal',S,'.RData'))
