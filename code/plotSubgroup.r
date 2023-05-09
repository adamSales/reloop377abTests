library(xtable)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggeffects)
library(clubSandwich)

source('code/reloopFunctions.r')
load('data/pairwiseData.RData')

load("results/subgroupResults.RData")
load("results/resTotalSlow.RData")

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

resSub2=NULL
for(i in 1:10) for(j in 1:2) resSub2=c(resSub2,resSubgroups[[i]][[j]])


resExcl=map(resTotal,
              function(res)
                res[
                  exl$excl[exl$model==res[[1]]$model[1]]&
                  exl$pval[exl$model==res[[1]]$model[1]]>0.1])

resSub2 <- resSub2[map_lgl(resSub2,~.$ps[1]%in%names(resExcl[[1]]))]

resSub2 <- resSub2[vapply(resSub2,function(x) any(is.finite(x$Var))&min(x$Var,na.rm=TRUE)>1e-10,TRUE)]

makeContrast=function(res,ols=TRUE){
  rl=if(ols) 'reloopOLS' else 'reloopRF'
  v=if('Var'%in%colnames(res[[1]])) 'Var' else 'var'
  cont1=sapply(res,
               function(x) if(nrow(x)<3) NA else x['simpDiff',v]/x[rl,v])

  cont2=sapply(res,
               function(x) if(nrow(x)<3) NA else x['simpDiff',v]/x['reloopPlus',v])

  cont3=sapply(res,
               function(x) if(nrow(x)<3) NA else x['loop',v]/x['reloopPlus',v])

  id=map_dfr(res,~.[1,c('ps','covName','side','nt','nc')])
  id$PS=map_chr(strsplit(id$ps,"(Control)|(Treatment)"),~.[1])
  id$num=1:length(res)


  ssMult=c(cont1,cont2,cont3)

  cbind(rbind(id,id,id),
    data.frame(
      comp=rep(c('/n$/\frac{\\varhat(\\tsd)}{\\varhat(\\trc)}$\\n',
                 '\n$\\frac{\\varhat(\\tsd)}{\\varhat(\\trcpen)}$\\n',
                 '\n$\\frac{\\varhat(\\tss[\\bx,\\mathrm{RF}])}{\\varhat(\\trcpen)}$\\n'),
               each=length(cont1)),
      compSimp=rep(c('reloopVsSD',
                 'reloopPlusVsSD',
                 'reloopPlusVsLoop'),
               each=length(cont1)),
      ssMult=ssMult,
      model=res[[1]]$model[1]   
    )
  )
}

ContrastsSub=makeContrast(resSub2) #map_dfr(resSub2,makeContrast)

save(ContrastsSub,file='results/contrastsSub.RData')

p1<-ContrastsSub%>%
  filter(model=='combined')%>%
  mutate(compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)))%>%
ggplot(aes(compTxt,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+
    geom_hline(yintercept=1)+xlab(NULL)+
  scale_y_continuous(trans="log10")

ContrastsSub%>%
mutate(compTxt=
              c(reloopVsSD='ReLOOP vs. T-Test',
                reloopPlusVsSD='ReLOOP+ vs. T-Test',
                reloopPlusVsLoop='ReLOOP+ vs. LOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)))%>%
ggplot(aes(ssMult))+
  geom_histogram(bins=100,fill=NA,color="black",boundary=1)+
    geom_vline(xintercept=1)+facet_wrap(~compTxt,ncol=1)+
  scale_x_continuous(trans="log10",breaks=c(.5,.75,1,1.25,1.5,1.75,2,2.5,3,5))+
  xlab("Sampling Variance Ratio")+
  theme(axis.text.x = element_text(angle = 45,hjust=1)
  )
ggsave("../JEDM/figure/subgroupHistograms.jpg",width=4,height=3)

ssMultC=cut(ContrastsSub$ssMult,c(0,.75,.9,1,1.1,1.25,1.5,2,10))
table(ssMultC,ContrastsSub$compSimp)

round(table(ssMultC,ContrastsSub$compSimp)/4311,3)

ContrastsSub%>%group_by(compSimp)%>%summarize(mean(ssMult<1))


rng=ContrastsSub%>%group_by(compSimp)%>%group_map(~boxplot.stats(.$ssMult)$stats)%>%vapply(range,c(1.1,2.2))%>%range()

ContrastsSub%>%
  filter(model=='combined')%>%
  filter(compSimp!='reloopPlusVsSD')%>%
  mutate(compTxt=
              c(reloopVsSD='ReLOOP vs.\nT-Test',
                reloopPlusVsSD='ReLOOP+ vs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+ vs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)),
         covName=sub("student_prior_","",covName,fixed=TRUE),
         covName=sub("problem_set","PS",covName),
         covName=sub("average","avg",covName),
         covName=sub("skill_builder","SB",covName),
         covName=sub("problem","prob.",covName),
         covName=sub("median","med.",covName))%>%
ggplot(aes(fct_reorder(covName, ssMult, median),ssMult))+
  geom_boxplot(outlier.shape=NA)+
    geom_hline(yintercept=1)+xlab(NULL)+
  coord_cartesian(ylim=rng)+
facet_wrap(~compTxt+side,ncol=4)+scale_y_continuous(name="Samping Variance Ratio",trans="log10",breaks=c(.75,.9,1,1.1,1.25,1.5,2))+
guides(x =  guide_axis(angle = 90))
ggsave("../JEDM/figure/subgroupBoxplots2.jpg",width=6.5,height=4)


p4b<-ContrastsSub%>%
  filter(model=='combined')%>%
  mutate(compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)))%>%
ggplot(aes(compTxt,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA)+
    geom_hline(yintercept=1)+xlab(NULL)+
facet_wrap(~covName+side)+scale_y_continuous(trans="log10",limits = c(NA,2))

modelDat=ContrastsSub%>%
transmute(logn=log(nt+nc),#hmeanN=log(2/(1/nt+1/nc)),
logY=log(ssMult),PS=PS)%>%
filter(is.finite(logY),is.finite(logn))

model=lm(logY~splines::bs(logn,4),data=modelDat)

pred=ggpredict(
  model, vcov.fun = "vcovCR", vcov.type = "CR0", 
  vcov.args = list(cluster = modelDat$PS)
)

Wald_test(model,vcov="CR0",cluster=modelDat$PS,constraints=constrain_zero("logn",reg_ex=TRUE))

ContrastsSub$n=ContrastsSub$nc+ContrastsSub$nt
modTTest=lm_robust(log(ssMult)~splines::bs(log(n)),data=subset(ContrastsSub,compSimp=='reloopVsSD'),cluster=PS)
predTTest=predict(modTTest,subset(ContrastsSub,compSimp=='reloopVsSD'),interval='confidence')
modLOOP=lm_robust(log(ssMult)~splines::bs(log(n)),data=subset(ContrastsSub,compSimp=='reloopPlusVsLoop'),cluster=PS)
predLOOP=predict(modTTest,subset(ContrastsSub,compSimp=='reloopPlusVsLoop'),interval='confidence')

ContrastsSub$ypred=NA
ContrastsSub$lwr=NA
ContrastsSub$upr=NA

ContrastsSub[ContrastsSub$compSimp=='reloopVsSD',c('ypred','lwr','upr')]=predTTest$fit
ContrastsSub[ContrastsSub$compSimp=='reloopPlusVsLoop',c('ypred','lwr','upr')]=predLOOP$fit


#p5 <- 

ContrastsSub%>%
  filter(compSimp!='reloopPlusVsSD')%>%
  mutate(compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)))%>%
ggplot(aes(log(nt+nc),#2/(1/nt+1/nc),
      log(ssMult),clust=PS))+
#scale_x_log10()+#limits=c(10,NA) ) +
 # scale_y_log10()+
geom_point(alpha=0.3)+
geom_ribbon(mapping=aes(ymin = lwr,ymax = upr), alpha = 1,fill='red')+
geom_line(mapping=aes(x=log(nt+nc),y=ypred),color='black')+#,linewidth=2)+
#geom_smooth(method="lm_robust",formula=y~splines::bs(x,4),method.args=list(clusters=ContrastsSub$PS[ContrastsSub$compSimp!='reloopPlusVsSD']))#+
facet_wrap(~compTxt)+#,scales="free_y")+
geom_hline(yintercept = 0)+
scale_y_continuous("Sampling Variance Ratio",
                  breaks=log(c(.75,.9,1,1.1,1.25,1.5,2,3,4,5,6,7)),
                  labels = c(.75,.9,1,1.1,1.25,1.5,2,3,4,5,6,7))+
scale_x_continuous("Total Sample Size",breaks=log(c(50,100,200,500,1000,5000,10000,50000)),labels=c(50,100,200,500,1000,5000,10000,50000))
ggsave("../JEDM/figure/subgroupSampleSize.jpg",width=6.5,height=3)

  #facet_wrap(~comp,nrow=1)+
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


#######################################
## p-values

names(resSub2)=1:length(resSub2)

resSub2=lapply(1:length(resSub2),function(i) cbind(resSub2[[i]],num=i))

estimates=function(res,estimator){
  if(nrow(res)<3) return(data.frame(Est=NA,SE=NA,p=NA))
  if(any(is.na(res$Var))) return(data.frame(Est=NA,SE=NA,p=NA))
  if(min(res$Var)<1e-10) return(data.frame(Est=NA,SE=NA,p=NA))
  Est=res[ estimator,'Est']
  SE=sqrt(res[ estimator,'Var'])
  if("nt"%in%names(res)){
    nt=res$nt[1]
    nc=res$nc[2]
    df=(1/nt+1/nc)^2/(1/nt^2/(nt-1)+1/nc^2/(nc-1))
  } else df=Inf

  data.frame(num=res$num[1],Est=Est,SE=SE,p=2*pt(-abs(Est/SE),df=df))
}

pvals=map(c( 'simpDiff','loop','reloopOLS','reloopPlus')%>%setNames(.,.),
  ~map_dfr(resSub2,estimates,estimator=.))

par(mfrow=c(4,1))
lapply(pvals,function(ee) hist(ee$p))  

pvals=map(pvals,~cbind(.,pBH=p.adjust(.$p,method="fdr")))
pvals=map(pvals,~cbind(.,pBY=p.adjust(.$p,method="BY")))

significance=rbind(
Unadjusted=map_dbl(pvals,~sum(.$p<0.05)),
`Benjamini-Hochberg`=map_dbl(pvals,~sum(.$pBH<0.05)),
`Benjamini-Yekutieli`=map_dbl(pvals,~sum(.$pBY<0.05)))

colnames(significance)=c("T-Test","LOOP","ReLOOP","ReLOOP+")

print(xtable(significance),file="subgroupSignificance.tex")


map(pvals,~sum(.$pBY<0.05))
map(pvals,~mean(.$pBY<0.05))


resExcl[[4]]=lapply(1:length(resExcl[[4]]),function(i) cbind(resExcl[[4]][[i]],num=i))


pvalsFull=map(c( 'simpDiff','loop','reloopOLS','reloopPlus')%>%setNames(.,.),
  ~map_dfr(resExcl[[4]],estimates,estimator=.))

pvalsFull=map(pvalsFull,~cbind(.,pBH=p.adjust(.$p,method="fdr"),pBY=p.adjust(.$p,method="BY")))

significanceFull=rbind(
Unadjusted=map_dbl(pvalsFull,~sum(.$p<0.05)),
`Benjamini-Hochberg`=map_dbl(pvalsFull,~sum(.$pBH<0.05)),
`Benjamini-Yekutieli`=map_dbl(pvalsFull,~sum(.$pBY<0.05)))
colnames(significanceFull)=c("T-Test","LOOP","ReLOOP","ReLOOP+")
print(xtable(significanceFull),file="../JEDM/fullSignificance.tex")



map( pvalsFull,~sum(.$p<0.05))
map( pvalsFull,~mean(.$p<0.05))

map( pvalsFull,~sum(.$pBH<0.05))
map( pvalsFull,~mean(.$pBH<0.05))

map( pvalsFull,~sum(.$pBY<0.05))
map( pvalsFull,~mean(.$pBY<0.05))


load('results/contrasts.RData')

p3<-Contrasts%>%
  filter(model=='combined')%>%
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
  scale_y_continuous(trans="log10")+geom_hline(yintercept=1)+ylab("Sampling Variance Ratio")+xlab(NULL)

ggsave('../JEDM/figure/mainResults.jpg', plot=p3,   width=6,height=4)


p4<-Contrasts%>%
  mutate(
    compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt))
  )%>%
  ggplot(aes(model,ssMult))+
  geom_jitter()+geom_boxplot(outlier.shape=NA,width=.5)+facet_wrap(~compTxt,nrow=1)+
  scale_y_continuous(trans="log10")+geom_hline(yintercept=1)+ylab("Sampling Variance Ratio")+xlab(NULL)+theme(axis.text.x = element_text(angle = 45,hjust=1))

ggsave('../JEDM/figure/modelResults.jpg', plot=p4,   width=6,height=4)

