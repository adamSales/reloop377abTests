library(tidyverse)
load('results/res.RData')
load('results/resCombined.RData')
load('processedData/pairwiseData.RData')

covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]
p=length(covNames)

exl=datPW%>%group_by(problem_set,male)%>%
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


keepPS <- intersect(
  exl%>%filter(inbetween,pval>0.1,male)%>%pull(problem_set)%>%unique(),
  exl%>%filter(inbetween,pval>0.1,!male)%>%pull(problem_set)%>%unique())

res <- res[sapply(res,function(x) !inherits(x,'try-error'))]
resExcl<- res[sapply(res,function(x) x$ps[1]%in%keepPS)]

resBoth <- resCombined[sapply(resCombined,function(x) x$ps[1]%in%keepPS)]

resExcl <- c(resExcl,resBoth)

#### contrasts of interest: reloop vs simpDiff, reloopPlus vs simpDiff, reloopPlus vs Loop

### look at "combined"
makeContrast2=function(res,ols=TRUE){
  rl=if(ols) 'reloopOLS' else 'reloopRF'
  v=if('Var'%in%colnames(res[[1]])) 'Var' else 'var'
  cont1=sapply(res,
               function(x) x['simpDiff',v]/x[rl,v])

  cont2=sapply(res,
               function(x) x['simpDiff',v]/x['reloopPlus',v])

  cont3=sapply(res,
               function(x) x['loop',v]/x['reloopPlus',v])
  ps=sapply(res,function(x) x$ps[1])

  model <- sapply(res,function(x) x$model[1])

  data.frame(
    ps=rep(ps,3),
    comp=rep(c('$\\frac{V(\\hat{\\tau}^{DM})}{V(\\hat{\\tau}_{LOOP}(\\hat{y}^r))}$',
               '$\\frac{V(\\hat{\\tau}^{DM})}{V(\\hat{\\tau}_{LOOP}(\\hat{y}^r,x))}$',
               '$\\frac{V(\\hat{\\tau}_{LOOP}(x))}{V(\\hat{\\tau}_{LOOP}(\\hat{y}^r,x))}$'),
             each=length(cont1)),
    compSimp=rep(c('reloopVsSD',
               'reloopPlusVsSD',
               'reloopPlusVsLoop'),
             each=length(cont1)),
    ssMult=c(cont1,cont2,cont3),
    model=c(model,model,model)

  )
}

ContrastsMO <- makeContrast2(resExcl)

save(ContrastsMO,file='contrasts.RData')


load('../../ethenExperiments/results/contrasts.RData')



p <- Contrasts%>%
  filter(model=='combined')%>%
  bind_rows(ContrastsMO)%>%
  mutate(
    compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt))
  )%>%
  ggplot(aes(model,ssMult))+
  geom_jitter()+ #geom_violin(draw_quantiles=.5,alpha=0.5,size=2)
  geom_boxplot(outlier.shape=NA,width=.5)+
  facet_wrap(~compTxt,nrow=1)+
  geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10")+ylab("Sampling Variance Ratio")+xlab(NULL)
ggsave('../../JEDM/figure/lak1.jpg', plot=p,   width=6,height=4)




ContrastsMO%>%
  pivot_wider(names_from='model',values_from='ssMult')%>%
  ggplot(mapping=aes(M,O))+geom_point()+geom_abline(aes(intercept=0,slope=1))+xlim(.96,1.65)+ylim(.96,1.65)+facet_wrap(~compSimp)



ContrastsMO <- exl%>%
  mutate(model=ifelse(male,'M','O'))%>%
  select(problem_set,model,n,nt,nc,pval)%>%
  right_join(ContrastsMO,by=c("problem_set"="ps","model"))

ContrastsMO%>%
  ggplot(aes(n,ssMult))+geom_point()+geom_smooth()+
  facet_grid(model~compSimp)

ContrastsMO%>%
  mutate(ESS=2/(1/nc+1/nt))%>%
  ggplot(aes(ESS,ssMult))+geom_point()+geom_smooth()+geom_hline(yintercept=1)+
  facet_grid(model~compSimp)

ContrastsMO%>%filter(ssMult<.985)%>%mutate(sid=substr(problem_set,1, 7))%>%group_by(sid)%>%arrange(ssMult)%>%
  select(-comp)%>%as.data.frame()


#nms <- read.csv('~/Downloads/names.csv')

#pmale <- mean(nms[,2]=='Male')
 
## from Ethan's email
inferredGender=c(
unknown = 24486,
male = 22397,
female = 18394,
mostly_male = 3138,
mostly_female = 1591,
andy = 1131)

pmale=sum(inferredGender[c('male','mostly_male')])/sum(inferredGender)

PS <- sapply(unique(ContrastsMO$problem_set),
             function(ps){
               resM <- res[[paste0('TRUE.',ps)]]
               resF <- res[[paste0('FALSE.',ps)]]
               Est=pmale*resM[,'Est']+(1-pmale)*resF[,'Est']
               Var=pmale^2*resM[,'Var']+(1-pmale)^2*resF[,'Var']
               cbind(Est=Est,Var=Var,resM[,c('method','ps')],model='PS')},
             simplify=FALSE)

ContrastsPS <- makeContrast2(PS)

psPlot <- ContrastsPS%>%
mutate(
  dumb="",
    compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt))
  )%>%
  ggplot(aes(dumb,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA,width=.5)+facet_wrap(~compTxt,nrow=1)+
  scale_y_continuous(trans="log10")+geom_hline(yintercept=1)+ylab("Sampling Variance Ratio")+xlab(NULL)#geom_violin(draw_quantiles=.5,alpha=0.5,size=2)+facet_wrap(~comp,nrow=1)+
  #geom_hline(yintercept=1)+
  #scale_y_continuous(trans="log10")

ggsave('../../JEDM/figure/postStrat.jpg', plot=psPlot,   width=6,height=4)


#tikz('ps.tex', width=4, height=3,standAlone=TRUE)
#print(psPlot)
#dev.off()
