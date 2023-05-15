######################
## main analysis plots
######################

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

ggsave('figure/mainResults.jpg', plot=p3,   width=6,height=4)


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

ggsave('figure/modelResults.jpg', plot=p4,   width=6,height=4)

######################
## subgroup analysis plots
######################

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


ContrastsSub%>%
  filter(compSimp!='reloopPlusVsSD')%>%
  mutate(compTxt=
              c(reloopVsSD='ReLOOP\nvs.\nT-Test',
                reloopPlusVsSD='ReLOOP+\nvs.\nT-Test',
                reloopPlusVsLoop='ReLOOP+\nvs.\nLOOP')[compSimp],
         compTxt=factor(compTxt,levels=unique(compTxt)))%>%
  select(compSimp,ssMult,n)%>%arrange(compSimp,n,ssMult)%>%filter(ssMult<0.8)

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
