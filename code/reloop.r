library(loop.estimator)
library(tidyverse)
library(parallel)
source('reloopFunctions.r')
load('pairwiseData.RData')

cl <- makeCluster(4)
ce <- clusterEvalQ(cl,source('reloopFunctions.r'))
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

resTotal <- parLapply(cl,
                      unique(datPW$model),
                      function(mod)
                        datPW%>%
                        filter(model==mod)%>%
                        split(.,.$problem_set)%>%
                        map(~try(
                              allEst(
                                Y=.x$completion_target,
                                Tr=.x$Z,Z=as.matrix(.x[,covNames]),
                                yhat=as.matrix(.x[,'completion_prediction']),
                                ps=.x$problem_set[1],
                                model=unique(.x$model))
                            )
                            )
                      )
save(resTotal,file='resTotal.RData')


#### hard part: we have 5 sample size-based restrictions (including "none") and 4 p-value-based restrictions, leading to 20 possibilities. DECISION (prior to seeing results! tho this is unverifiable--you have to trust me): report results for restrictive conditions (p>.1, excl) and discuss differences with other sets of restictions. Provide plots for all combinations in the appendix.

resExcl<- map(resTotal,
              function(res)
                res[
                  exl$excl[exl$model==res[[1]]$model[1]]&
                  exl$pval[exl$model==res[[1]]$model[1]]>0.1])

names(resExcl)=vapply(resExcl,function(x) x[[1]]$model[1],'a')

#### contrasts of interest: reloop vs simpDiff, reloopPlus vs simpDiff, reloopPlus vs Loop

### look at "combined"
makeContrast=function(res){
  cont1=sapply(res,
               function(x) x['simpDiff','Var']/x['reloopOLS','Var'])

  cont2=sapply(res,
               function(x) x['simpDiff','Var']/x['reloopPlus','Var'])

  cont3=sapply(res,
               function(x) x['loop','Var']/x['reloopPlus','Var'])
  ps=sapply(res,function(x) x$ps[1])
  
  data.frame(
    ps=rep(ps,3),
    comp=rep(c('\n$\\frac{\\varhat(\\tsd)}{\\varhat(\\trc)}$\\n',
               '\n$\\frac{\\varhat(\\tsd)}{\\varhat(\\trcpen)}$\\n',
               '\n$\\frac{\\varhat(\\tss[\\bx,\\mathrm{RF}])}{\\varhat(\\trcpen)}$\\n'),
             each=length(cont1)),
    ssMult=c(cont1,cont2,cont3),
    model=res[[1]]$model[1]
  )
}


Contrasts=map_dfr(resExcl,makeContrast)

Contrasts$model=factor(Contrasts$model,levels=c('action','student','assignment','combined'))

p<-Contrasts%>%
  ggplot(aes(ssMult))+
  geom_dotplot( method="histodot", binwidth = .02 )  +
    labs( x = "Relative Ratio of Sample Variances", y="" ) +
    geom_vline( xintercept = 1, col="red" ) +
  facet_grid(model~comp)+
    theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        text=element_text(size=12),
        strip.text=element_text(size=12,lineheight=0.5))


p2<-Contrasts%>%
  ggplot(aes(model,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+facet_wrap(~comp,nrow=1)+
  scale_y_continuous(trans="log10")

Contrasts <- Contrasts%>%
  group_by(model,comp)%>%
  mutate(q1=quantile(ssMult,.25),q3=quantile(ssMult,.75),iqr=q3-q1,
         outlier=ifelse(ssMult<q1-2*iqr,'low',
                 ifelse(ssMult>q3+2*iqr,'high','no')))%>%
  ungroup()


tikz('writeup/figure/combined.tex',width=6.4,height=2,standAlone=FALSE)
print(p)
dev.off()


############################
### performance vs within-sample
############################


datPW%>%
  group_by(problem_set)%>%
  summarize(rho=cor(completion_prediction,completion_target))%>%
  left_join(
    filter(Contrasts,model=='combined',
           comp=='\n$\\frac{\\varhat(\\tsd)}{\\varhat(\\trc)}$\\n'),
    by=c("problem_set"="ps")
  )%>%
  ggplot(aes(rho,ssMult))+geom_point()+geom_smooth(method="lm",formula=y~poly(x,2),se=FALSE)


############################
### performance vs CV
############################



res.05 <- datPW%>%
  filter(pval>0.05)%>%
  split(.,.$problem_set)%>%
  map(~try(allEst(Y=.x$completion_target,Tr=.x$Z,Z=as.matrix(.x[,covNames]),yhat=as.matrix(.x[,'completion_prediction']),ps=.x$problem_set[1])))

save(res.05,file='res.05.RData')

res.01  <- datPW%>%
  filter(pval>0.01)%>%
  split(.,.$problem_set)%>%
  map(~try(allEst(Y=.x$completion_target,Tr=.x$Z,Z=as.matrix(.x[,covNames]),yhat=as.matrix(.x[,'completion_prediction']),ps=.x$problem_set[1])))

save(res.01,file='res.01.RData')

res.05 <- datPW%>%
  filter(pval>0.05)%>%
  split(.,.$problem_set)%>%
  map(~try(allEst(Y=.x$completion_target,Tr=.x$Z,Z=as.matrix(.x[,covNames]),yhat=as.matrix(.x[,'completion_prediction']),ps=.x$problem_set[1])))

save(res.05,file='res.05.RData')

res.1 <- datPW%>%
  filter(pval>0.1)%>%
  split(.,.$problem_set)%>%
  map(~try(allEst(Y=.x$completion_target,Tr=.x$Z,Z=as.matrix(.x[,covNames]),yhat=as.matrix(.x[,'completion_prediction']),ps=.x$problem_set[1])))
save(res.1,file='res.1.RData')

## reloopOLS vs simpDiff
cont1=sapply(res.1,
             function(x) x['simpDiff','V2']/x['reloopOLS','V2'])

ns=sapply(res.1,function(x) c(nt=sum(datPW$problem_set==x$ps[1]&datPW$Z==1),nc=sum(datPW$problem_set==x$ps[1]&datPW$Z==0)))

hist(cont1,50)
hist(cont1[ns['nt',]>3&ns['nc',]>3],50)
hist(cont1[ns['nt',]>3&ns['nc',]>3],50)


plot(ns,cont1)

cont2=sapply(res.1,
             function(x) x['loop','V2']/x['reloopPlus','V2'])
hist(cont2[ns>20],50)
mean(cont2[ns>20])


cont3=sapply(res.1,
             function(x) x['simpDiff','V2']/x['reloopPlus','V2'])
hist(cont3[ns>20],50)
mean(cont3[ns>20])
