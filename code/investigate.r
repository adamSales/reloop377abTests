library(loop.estimator)
library(tidyverse)
source('reloopFunctions.r')
load('pairwiseData.RData')

#### covariates
covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]


############ trying things out
ps1=datPW$problem_set[1]
datS=filter(datPW,problem_set==ps1)

covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]

loop(Tr=datS$Z,Y=datS$completion_target,Z=as.matrix(datS[,covNames]))

loop_ols(Tr=datS$Z,Y=datS$completion_target,Z=as.matrix(datS[,covNames]))
loop_ols(Tr=datS$Z,Y=datS$completion_target,Z=as.matrix(datS[,c(covNames,'completion_prediction')]))

loop(Tr=datS$Z,Y=datS$completion_target,Z=as.matrix(datS[,c(covNames,'completion_prediction')]))

loop(Tr=datS$Z,Y=datS$completion_target,Z=as.matrix(datS[,covNames]),pred=reloop,yhat=datS[,'completion_prediction'])

#### lin-style analysis
simpDiff0=lm_robust(completion_target~Z,data=datS)
rebar0=lm_robust(completion_target-completion_prediction~Z,data=datS)
reloop0=lm_lin(completion_target~Z,covariates=~completion_prediction,data=datS)

## other covariates
form=as.formula(paste('~',paste(covNames,collapse='+')))

loop0=lm_lin(completion_target~Z,covariates=form,data=datS)
reloopPlus0=lm_lin(completion_target~Z,covariates=update(form,~.+completion_prediction),data=datS)

### lin sampling variances
linVars=map(list(sd=simpDiff0,reb=rebar0,rel=reloop0,loop=loop0,relPlus=reloopPlus0),
            ~.$std.error['Z']^2)


#### colinearity?
ccc=datS[,covNames]
cor(ccc)
pairs(ccc)
pca=prcomp(ccc,scale.=TRUE)
screeplot(pca)

loop_ols(Tr=datS$Z,Y=datS$completion_target,Z=pca$x[,1:4])

#### random forest prediction
library(randomForest)
rf1=randomForest(ccc[datS$Z==1,],datS$completion_target[datS$Z==1])
rf1a=randomForest(ccc[datS$Z==1,],factor(datS$completion_target[datS$Z==1]),importance=TRUE)
plot(predict(rf1),datS$completion_prediction[datS$Z==1])
boxplot(predict(rf1)~datS$completion_target[datS$Z==1])
## whoa!
par(mfrow=c(3,3))
bp=apply(ccc,2,function(x) boxplot(x~datS$completion_target))

par(mfrow=c(1,1))
ccc%>%
  mutate(Y=factor(datS$completion_target),Z=datS$Z)%>%
  filter(Z==1)%>%
  pivot_longer(cols=-c(Y,Z),names_to="covariate",values_to="val")%>%
  ggplot(aes(Y,val))+geom_violin()+facet_wrap(~covariate,scales="free_y")


ccc2=ccc%>%
  mutate(perCompSB=
           student_prior_completed_skill_builder_count/student_prior_started_skill_builder_count,
         perCompPS=
           student_prior_completed_problem_set_count/student_prior_started_problem_set_count,
         perCompSB=ifelse(is.finite(perCompSB),perCompSB,mean(perCompSB,na.rm=TRUE)),
         perCompPS=ifelse(is.finite(perCompPS),perCompPS,mean(perCompPS,na.rm=TRUE))
         )%>%
  select(
    -student_prior_completed_skill_builder_count,
    -student_prior_completed_problem_set_count)

summary(lm(datS$completion_target[datS$Z==1]~as.matrix(ccc2[datS$Z==1,])))
summary(lm(datS$completion_target[datS$Z==1]~as.matrix(ccc[datS$Z==1,])))

rf1b=randomForest(ccc2[datS$Z==1,],factor(datS$completion_target[datS$Z==1]),importance=TRUE)

rfall=randomForest(ccc2,factor(datS$completion_target),importance=TRUE)

library(arm)
par(mfrow=c(3,3))
for(i in 1:ncol(ccc2)) binnedplot(ccc2[[i]], datS$completion_target)
  

allEst(Y=datS$completion_target,Tr=datS$Z,Z=as.matrix(datS[,covNames]),yhat=as.matrix(datS[,'completion_prediction']),ps=datS$problem_set[1])
