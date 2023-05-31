
if(full | estMain){
  print('main estimation')
  resTotal <-
    datPW%>%
    split(~problem_set+model)%>%
    LAP(.,function(x) {
      #cat(".")
      try(
                                allEst(
                                  Y=x$completion_target,
                                  Tr=x$Z,
                                  Z=as.matrix(x[,covNames]),
                                  yhat=as.matrix(x[,'completion_prediction']),
                                  ps=x$problem_set[1],
                                  model=x$model[1],
                                  fast=FALSE
                                )
                              )
    }
            )

  save(resTotal,file='results/resTotalSlow.RData')

} else load('results/resTotalSlow.RData')
