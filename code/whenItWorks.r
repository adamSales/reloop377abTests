library(tidyverse)
library(tikzDevice)
library(broom)
library(estimatr)
library(randomForest)
load('data/pairwiseData.RData')

load('results/resTotalSlow.RData')
load('results/contrasts.RData')

cv0=read_csv('data/cross_validation_results.csv')
crosswalk <- read_csv('data/exp_norm_map.csv')

cvMeas=read_csv('data/cross_validation_metrics.csv')%>%
  select(-starts_with("problems_"))%>%
  rename_with(function(x) paste0(map_chr(strsplit(x,"_"),~.[2]),'CV'),starts_with('completion_'))


cv <- left_join(cv0,crosswalk,by=c("target_sequence"="normal_id"))%>%
  filter(!is.na(experiment_id))%>%
  group_by(experiment_id,model)%>%
  summarize(rhoCV=cor(completion_target,completion_prediction),
            #nCV=n(),
            #mseCV=mean((completion_target-completion_prediction)^2),
            meanYcv=mean(completion_target,na.rm=TRUE))

covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]
covForm=as.formula(paste('completion_target~',paste(covNames,collapse='+')))

datPWs <- split(datPW,datPW$problem_set)
rf <- lapply(datPWs, function(x) randomForest(covForm,data=subset(x,model=='combined')))
rfDat=data.frame(problem_set=names(datPWs),
                 rf=vapply(rf,function(x) x$rsq[length(x$rsq)],1.1))

save(rfDat,file='results/rfDat.RData')

## rfDat=datPW%>%
##   filter(model=='combined')%>%
##   nest(data=-problem_set)%>%
##   mutate(
##     rf=map(data,~randomForest(.x[,covNames],.x$completion_target)$rsq[500]))%>%
##   unnest(rf)

resDat=datPW%>%
  group_by(problem_set,model)%>%
  summarize(rhoTrue=cor(completion_prediction,completion_target),
            logn=log(n()),
            mseTrue=mean((completion_prediction-completion_target)^2),
            meanYcTrue=mean(completion_target[Z==0]),
            logitYcTrue=qlogis(meanYcTrue)
            )%>%
  left_join(
    Contrasts,
    by=c("problem_set"="ps",'model')
  )%>%
  mutate(experiment_id=map_chr(strsplit(problem_set,'Control|Treatment'),~.[1]))%>%
  left_join(cv)%>%
  left_join(filter(cvMeas,experiment_id!='None'))%>%
  left_join(rfDat)




### models

#### within-sample data

within=resDat%>%
  filter(!is.na(compSimp))%>%
  group_by(model,compSimp)%>%
  mutate(rho2=rhoTrue^2,
         across(c(rho2,logn,logitYcTrue,rf),scale))%>%
  nest()%>%
  mutate(
    fit=map(data,~lm(ssMult~rho2+logn+logitYcTrue+rf,data=.x)),
    robfit=map(data,~lm_robust(ssMult~rho2+logn+logitYcTrue+rf,data=.x,clusters=experiment_id)),
    tidied=map(robfit,tidy))%>%
  unnest(tidied)

within%>%
  ungroup()%>%
  filter(model=='combined',term!='(Intercept)')%>%
  select(-model,-data,-fit,-robfit,-c(conf.low:outcome))%>%
  print(n=Inf)


within%>%
  filter(term!='(Intercept)')%>%
  ggplot(aes(term,estimate,ymin=conf.low,ymax=conf.high))+
  geom_point()+geom_errorbar(width=0)+geom_hline(yintercept=0)+
  facet_grid(compSimp~model)


within%>%
  filter(term=='logn')%>%
  mutate(fv=map(fit,fitted),resid=map(fit,resid))%>%unnest(c(fv,resid))%>%
  ggplot(aes(fv,resid))+geom_point()+geom_hline(yintercept=0)+
  geom_smooth(se=F)+
  facet_wrap(~compSimp+model,scales='free')

within%>%
  filter(term=='logn',model=='student')%>%
  mutate(data=map(data,~select(.,rho2,logn,logitYcTrue,rf)),resid=map(fit,resid))%>%unnest(c(data,resid))%>%
  pivot_longer(cols=c(rho2,logn,logitYcTrue,rf),names_to='variable',values_to='regressor')%>%
  ggplot(aes(regressor,resid))+geom_point()+geom_hline(yintercept=0)+
  geom_smooth(se=F)+
  facet_wrap(~compSimp+variable,scales='free')



par(mfrow=c(3,4))
within%>%
  group_by(compSimp,model)%>%
  mutate(title=paste(model[1],compSimp[1]))%>%
  group_walk(~plot(.x$fit[[1]],which=1,main=.x$title[1]))


#### within-sample data
