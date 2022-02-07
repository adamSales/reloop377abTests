library(tidyverse)
library(arm)
library(estimatr)
library(lme4)

cvResults <- read_csv('cross_validation_results.csv')

AUC = function(probs, true_Y){
  probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
  val = unlist(probsSort$x)
  idx = unlist(probsSort$ix)
  
  roc_y = true_Y[idx];
  stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
  stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)
  
  auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
  return(auc)
}

overallAUCs <- cvResults%>%
  group_by(model)%>%
  summarize(
    auc=AUC(completion_prediction,completion_target)
  )

overallOLS <- lm_robust(completion_target~completion_prediction,data=subset(cvResults,model=='combined'),fixed_effects = target_sequence)



overallOLSmods <- lm_robust(completion_target~completion_prediction*model,data=cvResults,fixed_effects = target_sequence)

mlm <- lmer(completion_target~completion_prediction+(completion_prediction|target_sequence),data=subset(cvResults,model=='combined'))

hist(ranef(mlm)$target_sequence$completion_prediction)

bySeq <- cvResults%>%
  filter(model=='combined')%>%
  group_by(target_sequence)%>%
  summarize(
    n=n(),
    np=sum(completion_target),
    nn=n-np,
    auc=AUC(completion_prediction,completion_target),
    aucSE=se_auc(auc,np,nn),
    mse=mean((completion_prediction-completion_target)^2),
    r=cor(completion_prediction,completion_target),
    r2=r^2,
    mseSE=sd((completion_prediction-completion_target)^2)/sqrt(n)
  )%>%
  arrange(auc)

bySeq%>%ungroup()%>%filter(n>100)%>%summarize(vAUC=var(auc),EaucSE=mean(aucSE^2),tau2Hat=vAUC-EaucSE)

expResults <- read_csv('experiment_results.csv')

expMapping <- read_csv('exp_norm_map.csv')

bySeqModel <- cvResults%>%
  filter(model=='combined')%>%
  group_by(target_sequence)%>%
  summarize(
    n=n(),
    auc=AUC(completion_prediction,completion_target),
    mse=mean((completion_prediction-completion_target)^2)
  )%>%
  arrange(auc)

length(intersect(cvResults$target_sequence,expMapping$normal_id))

bySeq$experiment <- bySeq$target_sequence%in%expMapping$normal_id



with(bySeq,plot(n,auc,col=ifelse(experiment,'blue','red')))
abline(h=0.77)
with(bySeq,plot(n,mse,col=ifelse(experiment,'blue','red')))

with(bySeq,plot(n,r,col=ifelse(experiment,'blue','red')))
with(bySeq,plot(n,r2,col=ifelse(experiment,'blue','red')))




expComplete <- cvResults%>%
  filter(
    model=='combined',
    target_sequence%in%expMapping$normal_id
    )%>%
  mutate(
    sqErr=(completion_prediction-completion_target)^2
    )

mod <- lmer(
  log(sqErr)~(1|target_sequence),
  data=expComplete)

summary(mod)

anova(mod,lm(log(sqErr)~1,data=expComplete))


anova(lm(log(sqErr)~as.factor(target_sequence),data=expComplete))

############
### binned resids

library(arm)

with(subset(cvResults,model=='combined'),
     binnedplot(completion_prediction,completion_target))

with(subset(cvResults,model=='combined'),
     binnedplot(
       completion_prediction,
       completion_prediction-completion_target))

br <- with(subset(cvResults,model=='combined'),
     binned.resids(
       completion_prediction,
       completion_prediction-completion_target)[[1]])


par(mfrow=c(2,2))

for(mm in unique(cvResults$model)){
  with(subset(cvResults,model==mm),
       binnedplot(
         completion_prediction,
         completion_prediction-completion_target,main=mm))
}

pdf('remnantBinnedPlots.pdf',width=6.5,height=9)
par(mfrow=c(3,2))
for(id in intersect(cvResults$target_sequence,expMapping$normal_id)){
  with(subset(cvResults,model=='combined'&target_sequence==id),
       binnedplot(
         completion_prediction,
         completion_prediction-completion_target,main=id))
}
dev.off()
