if(full | estSub){
  print('subgroup estimation')

  resSubgroups <- LAP(
                      c(covNames,'completion_prediction')%>%setNames(.,.),
                      function(covName){
                        lapply(c(L=1/3,H=2/3), function(p) {
                          datPW$sub=if(p<.5) datPW[[covName]]<quantile(datPW[[covName]],p) else
                                                datPW[[covName]]>quantile(datPW[[covName]],p)
                          datPW%>%
                          filter(model=='combined')%>%
                          filter(sub)%>%
                          split(.,.$problem_set)%>%
                          map(function(x) if(min(table(x$Z))>9){
                            filename <- paste0('results/miniResults/', 
                                              x$problem_set[1],covName,
                                              ifelse(p<.5,"Low","High"),".RData")
                            if(file.exists(filename)){
                              load(filename)
                              return(out)
                            }
                            #print(x$problem_set[1])
                            out <- try(
                                cbind(
                                  allEst(
                                    Y=x$completion_target,
                                    Tr=x$Z,Z=as.matrix(x[,covNames]),
                                    yhat=as.matrix(x[,'completion_prediction']),
                                    ps=x$problem_set[1],
                                    model=unique(x$model),
                                    fast=fast),
                                  nt=sum(x$Z),
                                  nc=sum(1-x$Z),
                                  covName=covName,
                                  side=ifelse(p<.5,"Low","High")
                            )
                          )
                          save(out,file=filename)
                          out
                        } else data.frame(
                          nt=sum(x$Z),
                          nc=sum(1-x$Z),
                          covName=covName,
                          side=ifelse(p<.5,"Low","High"),
                          ps=x$problem_set[1]
                        )
                      )
                    }
                        )
                      }
        )

  sv <- try(save(resSubgroups,file=paste0("results/subgroupResults",S,".RData")))
  if(!keepIntermediateResults) file.remove("results/miniResults/*")

} else{
  print('skipping subgroup estimation')
  load(paste0("results/subgroupResults",S,".RData"))
}
