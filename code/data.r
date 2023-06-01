#library(tidyverse)

dat <- read_csv('data/experiment_results.csv')

crosswalk <- read_csv('data/exp_norm_map.csv')

dat <- left_join(dat,crosswalk,by=c("sequence_id"="experiment_id"))

datSp <- split(dat,dat$sequence_id,drop=TRUE)

datSp <- lapply(datSp,function(dd){
    priors <- read_csv(paste0('data/covariates/',dd$sequence_id[1],'/priors.csv'))
    alog <- read_csv(paste0('data/covariates/',dd$sequence_id[1],'/exp_alogs.csv'))
    dd <- left_join(dd,priors,by=c("user_id"="student_id"))%>%
        left_join(select(alog,student_id,assigned_condition),by=c("user_id"="student_id"))
    dd
})

#sapply(datSp,ncol)

dat <- reduce(datSp,bind_rows)

save(dat,file='processedData/combinedData.RData')

dat%>%select(starts_with('student_prior'))%>%sapply(function(x) mean(is.na(x)))

########### mean imputation
dat <- dat%>%
    group_by(sequence_id,class_id)%>%
    mutate(across(starts_with('student_prior'),~ifelse(is.na(.),mean(.,na.rm=TRUE),.)))%>%
    group_by(sequence_id)%>%
    mutate(across(starts_with('student_prior'),~ifelse(is.na(.),mean(.,na.rm=TRUE),.)))%>%
    ungroup()%>%
    mutate(across(starts_with('student_prior'),~ifelse(is.na(.),mean(.,na.rm=TRUE),.)))

dat%>%select(starts_with('student_prior'))%>%sapply(function(x) mean(is.na(x)))

covNames <- names(dat)[startsWith(names(dat),'student_prior')]

dat$problem_set <- dat$sequence_id

dat <- dat%>%filter(!is.na(assigned_condition),assigned_condition!="Not Assigned")%>%
    mutate(Z=ifelse(startsWith(assigned_condition,"Treatment"),1,0),
           ctl=ifelse(startsWith(assigned_condition,"Control"),1,0))

save(dat,file='processedData/combinedData.RData')

dat%>%group_by(sequence_id)%>%
    summarize(n=n(),p=mean(Z))%>%
    ggplot(aes(n,p))+geom_point()


#### look at all pairwise-comparisons
datSp <- split(dat,dat$sequence_id)

## get pairwise
Pairs=function(cc){
    mat=outer(cc,cc,paste,sep=';')
    tt=as.vector(mat[upper.tri(mat)])
    do.call("rbind",strsplit(tt,';'))

}


pairwise=function(dd){
    pp=Pairs(unique(dd$assigned_condition))

    map_dfr(1:nrow(pp),
            function(i)
                filter(dd,assigned_condition%in%pp[i,])%>%
                mutate(
                    comparison=paste(pp[i,],collapse=';'),
                    n1=sum(assigned_condition==pp[i,1]),
                       n2=sum(assigned_condition==pp[i,2]),
                       pval=2*pbinom(min(n1,n2),n1+n2,.5)
                       ))
}

datPW=datSp%>%
    map_dfr(pairwise)

datPW$problem_set<- paste0(datPW$problem_set,datPW$comparison)

datPW <- datPW%>%
  group_by(problem_set)%>%
  mutate(cond1=assigned_condition[1],Z=ifelse(assigned_condition==cond1,1,0))%>%
  select(-cond1)%>%
  ungroup()

save(datPW,file='processedData/pairwiseData.RData')


datPW%>%group_by(problem_set,comparison)%>%
  summarize(n=n1+n2,nn=n(),prop=n1/n)%>%
  ungroup()%>%
  summarize(all(n==nn))


## plot proportions assigned between categories
datPW%>%group_by(problem_set,comparison)%>%
  summarize(n=n1+n2,nn=n(),prop=n1/n)%>%
  ggplot(aes(n,prop))+geom_point()+geom_hline(yintercept=c(1/6,1/3,1/2,2/3,5/6),linetype="dotted")+
  annotate("text",x=5000,y=c(1/6,1/3,1/2,2/3,5/6),label=c(1/6,1/3,1/2,2/3,5/6))


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
  filter(vY1>0,vY0>0,pval>0.1,excl)


### exclusion criteria:
### inclusive: p+3 observations in each condition, cuz with k predictors you need n>=k+1 to fit ols
### and we have p+1 predictors (including remnant predictions) and we are dropping one observation
### each time we fit so....
### exclusive: min(n_C,n_T)>=(p+2)*5+1 so you have 5 samples per predictor--popular rule of thumb(
### in-between: min(n_C,n_T)>=10+1 so you have 5 samples per predictor in ReLoop (ie just remnant)

#### hard part: we have 5 sample size-based restrictions (including "none") and 4 p-value-based restrictions, leading to 20 possibilities. DECISION (prior to seeing results! tho this is unverifiable--you have to trust me): report results for restrictive conditions (p>.1, excl) and discuss differences with other sets of restictions. Provide plots for all combinations in the appendix.


datPW <- filter(datPW,problem_set%in%exl$problem_set)

save(datPW,file='processedData/pairwiseData.RData')



############### inferred gender experiment

expMale <- read_csv('data/male_experiment_results.csv')
expOth <- read_csv('data/non_male_experiment_results.csv')


dim(expMale)
dim(expOth)

names(expOth)

expMale <- expMale%>%
  filter(!is.na(user_id),model=='combined')%>%
  group_by(sequence_id,user_id)%>%
  slice_head()%>%
  ungroup()%>%
  mutate(male="M")

expOth <- expOth%>%
  filter(!is.na(user_id),model=='combined')%>%
  group_by(sequence_id,user_id)%>%
  slice_head()%>%
  ungroup()%>%
  mutate(male="F")

datGen <- bind_rows(expMale,expOth)



datPWg <- inner_join(
  select(datPW,-completion_prediction),
  select(datGen,model,sequence_id,assignment_log_id,user_id,completion_prediction,male)
)

exlg=datPWg%>%group_by(problem_set,male)%>%
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


keepPSg <- intersect(
  exlg%>%filter(inbetween,pval>0.1,male=="M")%>%pull(problem_set)%>%unique(),
  exlg%>%filter(inbetween,pval>0.1,male=="F")%>%pull(problem_set)%>%unique())

datPWg <- filter(datPWg,problem_set%in%keepPSg)

save(datPWg,file='processedData/pairwiseDataGender.RData')
