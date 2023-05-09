library(loop.estimator)
#library(tidyverse)
library(dplyr)
library(purrr)
library(parallel)
source('code/reloopFunctions.r')
load('data/pairwiseData.RData')


cl <- makeCluster(10)
ce <- clusterEvalQ(cl,source('code/reloopFunctions.r'))
ce <- clusterEvalQ(cl,library(loop.estimator))
ce <- clusterEvalQ(cl,library(dplyr))
ce <- clusterEvalQ(cl,library(purrr))


#### covariates
covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]
p=length(covNames)

exl=datPW%>%group_by(problem_set,model)%>%
  summarize(
    n=n(),
    nt=sum(Z),
    nc=n-nt,
    vY=var(completion_target,na.rm=TRUE),
    vY1=var(completion_target[Z==1],na.rm=TRUE),
    vY0=var(completion_target[Z==0],na.rm=TRUE),
    pval=2*(pbinom(min(nt,nc)-1,n,.5)+dbinom(min(nt,nc),n,.5)/2),
    inclReloop=min(nt,nc)>=3,
    incl=min(nt,nc)>=p+3,
    inbetween=min(nc,nt)>=11,
    excl=min(nc,nt)>=(p+2)*5+1
  )%>%
  filter(vY1>0,vY0>0)

datPW <- filter(datPW,problem_set%in%exl$problem_set)
datPW <- filter(datPW,model=='combined')


clusterExport(cl,"datPW")
clusterExport(cl,"covNames")


resSubgroups <- parLapply(cl,
                      c(covNames,'completion_prediction')%>%setNames(.,.),
                      function(covName){
                        lapply(c(L=1/3,H=2/3), function(p) {
                          datPW$sub=if(p<.5) datPW[[covName]]<quantile(datPW[[covName]],p) else
                                                datPW[[covName]]>quantile(datPW[[covName]],p)
                          datPW%>%
                          #filter(problem_set=="PSA25TAControl;Treatment 1")%>%
                          filter(sub)%>%
                          split(.,.$problem_set)%>%
                          map(~if(min(table(.x$Z))>9){
                              try(
                                cbind(
                                  allEst(
                                    Y=.x$completion_target,
                                    Tr=.x$Z,Z=as.matrix(.x[,covNames]),
                                    yhat=as.matrix(.x[,'completion_prediction']),
                                    ps=.x$problem_set[1],
                                    model=unique(.x$model)),
                                  nt=sum(.x$Z),
                                  nc=sum(1-.x$Z),
                                  covName=covName,
                                  side=ifelse(p<.5,"Low","High")
                            )
                          )
                        } else data.frame(
                          nt=sum(.x$Z),
                          nc=sum(1-.x$Z),
                          covName=covName,
                          side=ifelse(p<.5,"Low","High"),
                          ps=.x$problem_set[1]
                        )
                      ) 
                    }
                        )
                      }     
        )

save(resSubgroups,file="results/subgroupResults.RData")


resSub2=NULL
for(i in 1:10) for(j in 1:2) resSub2=c(resSub2,resSubgroups[[i]][[j]])
