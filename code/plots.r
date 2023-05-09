library(tidyverse)
library(tikzDevice)
load('data/pairwiseData.RData')

load('results/resTotalSlow.RData')
load('results/contrasts.RData')

#pdf('plots.pdf')



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
  geom_violin(outlier.shape=NA)+geom_jitter()+facet_wrap(~compSimp,nrow=1)+
  scale_y_continuous(trans="log10")


modContrast

p3<-Contrasts%>%
  filter(model=='combined')%>%
  mutate(
    compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt))
  )%>%
  ggplot(aes(model,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+facet_wrap(~compTxt,nrow=1)+
  scale_y_continuous(trans="log10")+geom_hline(yintercept=1)+ylab("Sampling Variance Ratio")+xlab(NULL)

p3<-Contrasts%>%
  filter(model=='combined')%>%
  mutate(
    compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt))
  )%>%
  ggplot(aes(model,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+facet_wrap(~compTxt,nrow=1)+
  scale_y_continuous(trans="log10")+geom_hline(yintercept=1)+ylab("Sampling Variance Ratio")+xlab(NULL)

p4<-Contrasts%>%
  filter(model=='combined')%>%
  mutate(compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)))%>%
ggplot(aes(compTxt,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+#facet_wrap(~comp,nrow=1)+
  scale_y_continuous(name= expression(SE(Alternative)^2~"/SE(ReLOOP)"^2),
                     sec.axis=sec_axis(~.*2,
                                       name=expression(#""%<-%retteB~POOLeR),#
                                         ReLOOP~Better%->%""),
                                       breaks=NULL,labels=NULL,guide=NULL))+
  theme(axis.title.y.right = element_text(angle = 90),text = element_text(size = 20))+
  geom_hline(yintercept=1)+xlab(NULL)+
  ggtitle("Ratios of Sampling Variances")

pdf('results.pdf',width=2,height=2.5)#,units='in',res=72)
print(p4)
dev.off()




options(
       tikzLatexPackages = c(
        getOption('tikzLatexPackages'),
        '\\usepackage[T1]{fontenc}',
        '\\usepackage{graphicx}',
        '\\usepackage{amsmath,amsfonts,amssymb}',
        '\\usepackage{bm}',
        '\\usepackage[utf8]{inputenc}',
        '\\input{g:/My Drive/iesJohannReloop/ethenExperiments/writeup/notation.tex}'))


tikz('aied.tex',width=4.7,height=2)
p3
dev.off()


p2tot <- ContrastsTot%>%
  ggplot(aes(model,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+facet_wrap(~comp,nrow=1)+
  scale_y_continuous(trans="log10")

print(p2)

Contrasts <- Contrasts%>%
  group_by(model,comp)%>%
  mutate(q1=quantile(ssMult,.25),q3=quantile(ssMult,.75),iqr=q3-q1,
         outlier=ifelse(ssMult<q1-2*iqr,'low',
                 ifelse(ssMult>q3+2*iqr,'high','no')))%>%
  ungroup()


#############################
## pred vs. prior percent correct
#############################
studPPC=datPW%>%
  filter(model=='combined')%>%
  group_by(user_id)%>%
  summarize(priorPercentComplete=mean(student_prior_completed_skill_builder_count/student_prior_started_skill_builder_count),pred=mean(completion_prediction),predSD=sd(completion_prediction),nExp=n_distinct(sequence_id),
            nsb=mean(student_prior_started_skill_builder_count))

with(studPPC,cor(priorPercentComplete,pred,use='pairwise'))
with(studPPC,cor(priorPercentComplete,pred,use='pairwise',method='spearman'))

hist(studPPC$predSD[studPPC$nExp>1])

studPPC%>%
  ggplot(aes(priorPercentComplete,pred))+geom_point()+geom_smooth(se=FALSE)+ylim(0,1)+geom_abline(slope=1,intercept=0)

studPPC%>%
  filter(nExp==1)%>%
  ggplot(aes(priorPercentComplete,pred))+geom_point()+geom_smooth(se=FALSE)+ylim(0,1)+geom_abline(slope=1,intercept=0)

with(filter(studPPC,nExp==1),cor(priorPercentComplete,pred,use='pairwise'))

mod1=lm(pred~priorPercentComplete,data=filter(studPPC,is.finite(priorPercentComplete),nExp==1))
mod2=lm(pred~priorPercentComplete*nsb,data=filter(studPPC,is.finite(priorPercentComplete),nExp==1))

############################
### performance vs within-sample
############################


p=datPW%>%
  filter(model=='combined')%>%
  group_by(problem_set)%>%
  summarize(rho=cor(completion_prediction,completion_target))%>%
  left_join(
    filter(Contrasts,model=='combined',
          compSimp=='reloopVsSD'),
    by=c("problem_set"="ps")
  )%>%
#  filter(rho<0.2,ssMult>1.3)%>%pull(problem_set)
  ggplot(aes(rho^2,ssMult))+geom_point()+geom_smooth(method="rlm",formula=y~poly(x,2),se=FALSE)+
  labs(y='V(SD)/V(ReloopOLS)')
print(p)


p=datPW%>%
  group_by(problem_set)%>%
  summarize(rho=cor(completion_prediction,completion_target))%>%
  left_join(
    filter(Contrasts,model=='combined',
           compSimp=='reloopPlusVsLoop'),
    by=c("problem_set"="ps")
  )%>%
#  filter(rho<0.2,ssMult>1.3)%>%pull(problem_set)
  ggplot(aes(rho^2,ssMult))+geom_point()+geom_smooth(method="lm")+#,formula=y~poly(x,2),se=FALSE)+
  labs(y='V(Loop)/V(ReloopPlus)')
print(p)


p=datPW%>%
  group_by(problem_set)%>%
  summarize(logn=log(n()))%>%
  right_join(
    filter(Contrasts,model=='combined',
          compSimp=='reloopVsSD'),
    by=c("problem_set"="ps")
  )%>%
#  filter(rho<0.2,ssMult>1.3)%>%pull(problem_set)
  ggplot(aes(logn,ssMult))+geom_point()+geom_smooth()+
  labs(y='V(SD)/V(ReloopOLS)')
print(p)


p=datPW%>%
  group_by(problem_set)%>%
  summarize(meanYc=mean(completion_target[Z==0]),
            logitYc=qlogis(meanYc))%>%
  right_join(
    filter(Contrasts,model=='combined',
          compSimp=='reloopVsSD'),
    by=c("problem_set"="ps")
  )%>%
#  filter(rho<0.2,ssMult>1.3)%>%pull(problem_set)
  ggplot(aes(logitYc,ssMult))+geom_point()+geom_smooth()+
  labs(y='V(SD)/V(ReloopOLS)')
print(p)


p=datPW%>%
  group_by(problem_set)%>%
  summarize(meanYcHat=mean(completion_prediction[Z==0]),
            logitYcHat=qlogis(meanYcHat))%>%
  right_join(
    filter(Contrasts,model=='combined',
          compSimp=='reloopVsSD'),
    by=c("problem_set"="ps")
  )%>%
#  filter(rho<0.2,ssMult>1.3)%>%pull(problem_set)
  ggplot(aes(logitYcHat,ssMult))+geom_point()+geom_smooth()+
  labs(y='V(SD)/V(ReloopOLS)')
print(p)

p=datPW%>%
  group_by(problem_set)%>%
  summarize(logn=log(n()))%>%
  right_join(
    filter(Contrasts,model=='combined',
           compSimp=='reloopPlusVsLoop'),
    by=c("problem_set"="ps")
  )%>%
#  filter(rho<0.2,ssMult>1.3)%>%pull(problem_set)
  ggplot(aes(logn,ssMult))+geom_point()+geom_smooth()+
  scale_x_continuous(breaks=6:10,labels=round(exp(6:10),-3))+
  labs(y='V(Loop)/V(ReloopPlus)')
print(p)

Contrasts%>%
  group_by(ps)%>%
  summarize(loopVsSD=ssMult[compSimp=='reloopPlusVsSD']/ssMult[compSimp=='reloopPlusVsLoop'],reloopVsSD=ssMult[compSimp=='reloopVsSD'])%>%
  ggplot(aes(loopVsSD,reloopVsSD))+geom_point()+geom_smooth()

############################
### performance vs CV
############################
cv0=read_csv('data/cross_validation_results.csv')
crosswalk <- read_csv('data/exp_norm_map.csv')



cv <- left_join(cv0,crosswalk,by=c("target_sequence"="normal_id"))%>%
  filter(!is.na(experiment_id))%>%
  group_by(experiment_id)%>%
  summarize(rho=cor(completion_target,completion_prediction),n=n(),mse=mean((completion_target-completion_prediction)^2))

p=Contrasts%>%
  filter(model=='combined',compSimp=='reloopVsSD')%>%
  mutate(experiment_id=map_chr(strsplit(ps,'Control|Treatment'),~.[1]))%>%
  group_by(experiment_id)%>%
  summarize(avgSSmult=mean(ssMult))%>%
  left_join(cv)%>%
  ggplot(aes(rho,avgSSmult))+geom_point()+geom_smooth(se=FALSE)+
  labs(y='V(Loop)/V(ReloopPlus) (mean)',x='out of sample cor(y,yhat)')
print(p)

p=Contrasts%>%
  filter(model=='combined')%>%
  mutate(experiment_id=map_chr(strsplit(ps,'Control|Treatment'),~.[1]))%>%
  group_by(experiment_id)%>%
  summarize(avgSSmult=mean(ssMult))%>%
  left_join(cv)%>%
  ggplot(aes(mse,avgSSmult))+geom_point()+geom_smooth(se=FALSE)+
  labs(y='V(Loop)/V(ReloopPlus) (mean)',x='out of sample mse(y,yhat)')
print(p)

cvMeas=read_csv('data/cross_validation_metrics.csv')

p=Contrasts%>%
  filter(model=='combined')%>%
  mutate(experiment_id=map_chr(strsplit(ps,'Control|Treatment'),~.[1]))%>%
  group_by(experiment_id)%>%
  summarize(avgSSmult=mean(ssMult))%>%
  left_join(cvMeas)%>%filter(model=='combined')%>%
  ggplot(aes(completion_mse,avgSSmult))+geom_point()+geom_smooth(se=FALSE)+
  labs(y='V(Loop)/V(ReloopPlus) (mean)',x='out of sample mse(y,yhat)')
print(p)

p=Contrasts%>%
  filter(model=='combined')%>%
  mutate(experiment_id=map_chr(strsplit(ps,'Control|Treatment'),~.[1]))%>%
  group_by(experiment_id)%>%
  summarize(avgSSmult=mean(ssMult))%>%
  left_join(cvMeas)%>%filter(model=='combined')%>%
  mutate(logSampleSize=log(sample_size))%>%
  select(avgSSmult,logSampleSize,starts_with('completion_'))%>%
  pivot_longer(-avgSSmult,names_to='meas',values_to='measure')%>%
  ggplot(aes(measure,avgSSmult))+geom_point()+geom_smooth(method='lm')+
  labs(y='V(Loop)/V(ReloopPlus) (mean)')+facet_wrap(~meas,scales="free_x")
print(p)
### out of sample vs in sample rho

p=datPW%>%
  group_by(sequence_id,user_id)%>%
  summarize(completion_prediction=mean(completion_prediction),completion_target=mean(completion_target))%>%
  #xtabs(~completion_target,data=.)
  group_by(sequence_id)%>%
  summarize(rhoInSample=cor(completion_prediction,completion_target))%>%
  left_join(rename(cv,rhoOutOfSample=rho,sequence_id=experiment_id))%>%
  mutate(cor=paste('cor=',round(cor(rhoOutOfSample,rhoInSample,use='pairwise'),3)),pval=paste('p=',round(cor.test(rhoOutOfSample,rhoInSample)$p.value,3)),x=c(.3,rep(NA,length(sequence_id)-1)),y1=c(.65,rep(NA,length(sequence_id)-1)),y2=c(.6,rep(NA,length(sequence_id)-1)))%>%
  ggplot(aes(rhoOutOfSample,rhoInSample))+geom_point()+geom_smooth(method='lm')+geom_text(aes(x=x,y=y1,label=cor),check_overlap=TRUE)+geom_text(aes(x=x,y=y2,label=pval),check_overlap=TRUE)
print(p)

dev.off()
