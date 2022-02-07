simpDiff <- function(Y,Tr,...){
  est=mean(Y[Tr==1])-mean(Y[Tr==0])
  v=var(Y[Tr==1])/sum(Tr)+var(Y[Tr==0])/sum(Tr==0)
  c(est,v)
}

rebar <- function(Y,Tr,yhat,...)
  simpDiff(Y-yhat,Tr)

allEst <- function(Y,Tr,Z,yhat,ps="ps",model="mod",fast=TRUE){

  stopifnot(length(model)==1)
  stopifnot(length(ps)==1)
  out <- as.data.frame(
    rbind(
      simpDiff=simpDiff(Y=Y,Tr=Tr),
      rebar=rebar(Y=Y,Tr=Tr,yhat=yhat),
      reloopOLS=loop(Y=Y,Tr=Tr,Z=cbind(yhat=yhat),pred=loop_ols),
      reloopRF=loop(Y=Y,Tr=Tr,Z=cbind(yhat=yhat),pred=loop_rf),
      loop=loop(Y=Y,Tr=Tr,Z=Z),
      reloopPlus=loop(Y=Y,Tr=Tr,Z=Z,pred=if(fast) reloop else reloop_slow,yhat=yhat)
    ))
  names(out)=c('Est','Var')
  out$method=rownames(out)
  out$ps=ps
  out$model=model
#  attr(out,"ps") <- ps
  out
}

