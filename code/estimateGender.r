if(full | estGen){
  print('gender estimation')
  resGender <-
    datPWg%>%
    bind_rows(mutate(datPWg,male="B"))%>%
    split(~male+problem_set)%>%
    LAP(.,function(x) {
      #cat("."
      filename <- paste0('results/miniResults/',x$problem_set[1],
                          x$male[1],'.RData')
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
                                  model=x$male[1],
                                  fast=fast
                                )
                              )
      save(out,file=filename)
      out
    }
            )

  save(resGender,file='results/resGender.RData')

} else load('results/resGender.RData')
