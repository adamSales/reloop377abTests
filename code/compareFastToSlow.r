library(tidyverse)

print(load('results/resTotalSlow.RData'))
slowTot <- resTotal

print(load('results/resTotalFast.RData'))
fastTot <- resTotal


length(slow)
length(fast)

slow <- map_dfr(slowTot,~.['reloopPlus',])%>%
  select(-Est)%>%rename(slow=Var)

fast <- map_dfr(fastTot,~.['reloopPlus',])%>%
  select(-Est)%>%rename(fast=Var)


comp <- full_join(slow,fast)

comp%>%ggplot(aes(slow,fast))+geom_point()+geom_abline(slope=1,intercept=0)
comp%>%ggplot(aes(sqrt(slow),sqrt(fast)))+geom_point()+geom_abline(slope=1,intercept=0)

comp%>%ggplot(aes(slow,fast-slow))+geom_point()+geom_abline(slope=0,intercept=0)

comp%>%ggplot(aes(fast-slow))+geom_histogram(fill=NA,color="black")

comp%>%ggplot(aes(slow/fast))+geom_histogram(fill=NA,color="black")+scale_x_continuous(trans="log10")

with(comp,t.test(fast-slow))


comp%>%ggplot(aes(slow,fast))+geom_point()+geom_abline(slope=1,intercept=0)+facet_wrap(~model)
comp%>%ggplot(aes(sqrt(slow),sqrt(fast)))+geom_point()+geom_abline(slope=1,intercept=0)+facet_wrap(~model)
ggsave('figure/slowVSfastScatter.jpg')


comp%>%ggplot(aes(slow,fast-slow))+geom_point()+geom_abline(slope=0,intercept=0)+facet_wrap(~model)

comp%>%ggplot(aes(fast-slow))+geom_histogram(fill=NA,color="black")+facet_wrap(~model)

comp%>%ggplot(aes(slow/fast))+geom_histogram(fill=NA,color="black")+facet_wrap(~model)+scale_x_continuous(trans="log10")+geom_vline(xintercept=1)
ggsave('figure/slowVSfastSSmult.jpg')



comp%>%group_by(model)%>%group_map(~t.test(.$fast-.$slow))

comp%>%group_by(model)%>%group_map(~wilcox.test(.$fast,.$slow,paired=TRUE))
