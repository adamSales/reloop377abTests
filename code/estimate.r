load('processedData/pairwiseData.RData')


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

if(nclust>0){
  clusterExport(cl,"datPW")
  clusterExport(cl,"covNames")
}


#################################################################
#### Main Analysis
#################################################################
LAP <- if(nclust>0) function(X,FUN,...) parLapply(cl,X,FUN,...) else function(X,FUN,...) lapply(X,FUN,...)

estMain <- TRUE
if(file.exists('results/resTotalSlow.RData'))
 if(file.mtime('results/resTotalSlow.RData')> 
    max(
      file.mtime('code/estimate.r'),
      file.mtime('processedData/pairwiseData.RData')
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE

if(full | estMain){
  print('main estimation')
  resTotal <-
    datPW%>%
    split(~problem_set+model)%>%
    LAP(.,function(.x) try(
                                allEst(
                                  Y=.x$completion_target,
                                  Tr=.x$Z,
                                  Z=as.matrix(.x[,covNames]),
                                  yhat=as.matrix(.x[,'completion_prediction']),
                                  ps=.x$problem_set[1],
                                  model=.x$model[1],
                                  fast=FALSE
                                )
                              )
            )

  save(resTotal,file='results/resTotalSlow.RData')

} else load('results/resTotalSlow.RData')

#### hard part: we have 5 sample size-based restrictions (including "none") and 4 p-value-based restrictions, leading to 20 possibilities. DECISION (prior to seeing results! tho this is unverifiable--you have to trust me): report results for restrictive conditions (p>.1, excl) and discuss differences with other sets of restictions. Provide plots for all combinations in the appendix.

resExcl<- map(resTotal,
              function(res)
                res[
                  exl$excl[exl$model==res[[1]]$model[1]]&
                  exl$pval[exl$model==res[[1]]$model[1]]>0.1])

names(resExcl)=vapply(resExcl,function(x) x[[1]]$model[1],'a')


Contrasts=map_dfr(resExcl,makeContrast)

Contrasts$model=factor(Contrasts$model,levels=c('action','student','assignment','combined'))

save(Contrasts,file='results/contrasts.RData')


#################################################################
#### Subgroup Analysis
#################################################################

estSub <- TRUE
if(file.exists('results/subgroupResults.RData'))
 if(file.mtime('results/subgroupResults.RData')> 
    max(
      file.mtime('code/estimate.r'),
      file.mtime('processedData/pairwiseData.RData')
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE

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

} else{
  print('skipping subgroup estimation')
  load('results/subgroupResults.RData')
}

resSub2=NULL
for(i in 1:10) for(j in 1:2) resSub2=c(resSub2,resSubgroups[[i]][[j]])

resSub2 <- resSub2[map_lgl(resSub2,~.$ps[1]%in%names(resExcl[[1]]))]

resSub2 <- resSub2[vapply(resSub2,function(x) any(is.finite(x$Var))&min(x$Var,na.rm=TRUE)>1e-10,TRUE)]

ContrastsSub=makeContrast(resSub2) #map_dfr(resSub2,makeContrast)

save(ContrastsSub,file='results/contrastsSub.RData')
