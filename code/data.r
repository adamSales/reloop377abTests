library(tidyverse)

dat <- read_csv('experiment_results.csv')

crosswalk <- read_csv('exp_norm_map.csv')

dat <- left_join(dat,crosswalk,by=c("sequence_id"="experiment_id"))

datSp <- split(dat,dat$sequence_id,drop=TRUE)

datSp <- lapply(datSp,function(dd){
    priors <- read_csv(paste0('covariates/',dd$sequence_id[1],'/priors.csv'))
    alog <- read_csv(paste0('covariates/',dd$sequence_id[1],'/exp_alogs.csv'))
    dd <- left_join(dd,priors,by=c("user_id"="student_id"))%>%
        left_join(select(alog,student_id,assigned_condition),by=c("user_id"="student_id"))
    dd
})

#sapply(datSp,ncol)

dat <- reduce(datSp,bind_rows)

save(dat,file='combinedData.RData')

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

save(dat,file='combinedData.RData')

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

save(datPW,file='pairwiseData.RData')


datPW%>%group_by(problem_set,comparison)%>%
  summarize(n=n1+n2,nn=n(),prop=n1/n)%>%
  ungroup()%>%
  summarize(all(n==nn))


## plot proportions assigned between categories
datPW%>%group_by(problem_set,comparison)%>%
  summarize(n=n1+n2,nn=n(),prop=n1/n)%>%
  ggplot(aes(n,prop))+geom_point()+geom_hline(yintercept=c(1/6,1/3,1/2,2/3,5/6),linetype="dotted")+
  annotate("text",x=5000,y=c(1/6,1/3,1/2,2/3,5/6),label=c(1/6,1/3,1/2,2/3,5/6))


