library(loop.estimator)
library(tidyverse)
library(parallel)
source('../ethenExperiments/code/reloopFunctions.r')
load('pairwiseData.RData')

cl <- makeCluster(4)
ce <- clusterEvalQ(cl,source('../ethenExperiments/code/reloopFunctions.r'))
ce <- clusterEvalQ(cl,library(loop.estimator))
ce <- clusterEvalQ(cl,library(tidyverse))



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


### exclusion criteria:
### inclusive: p+3 observations in each condition, cuz with k predictors you need n>=k+1 to fit ols
### and we have p+1 predictors (including remnant predictions) and we are dropping one observation
### each time we fit so....
### exclusive: min(n_C,n_T)>=(p+2)*5+1 so you have 5 samples per predictor--popular rule of thumb(
### in-between: min(n_C,n_T)>=10+1 so you have 5 samples per predictor in ReLoop (ie just remnant)

datPW <- filter(datPW,problem_set%in%exl$problem_set)



clusterExport(cl,"datPW")
clusterExport(cl,"covNames")


##################### now for keeps

res <-
  datPW%>%
  filter(model=='combined')%>%
  split(~male+problem_set)%>%
  parLapply(cl,.,function(.x) try(
                              allEst(
                                Y=.x$completion_target,
                                Tr=.x$Z,Z=as.matrix(.x[,covNames]),
                                yhat=as.matrix(.x[,'completion_prediction']),
                                ps=.x$problem_set[1],
                                male=mean(.x$male)))
                            )

save(resTotal,file='resTotal.RData')
