simpDiff <- function(Y,Tr,...){
  est=mean(Y[Tr==1])-mean(Y[Tr==0])
  v=var(Y[Tr==1])/sum(Tr)+var(Y[Tr==0])/sum(Tr==0)
  c(est,v)
}

rebar <- function(Y,Tr,yhat,...)
  simpDiff(Y-yhat,Tr)

allEst <- function(Y,Tr,Z,yhat,ps="ps",model="mod",condNs=FALSE,fast=TRUE){

  stopifnot(length(model)==1)
  stopifnot(length(ps)==1)
  out <- as.data.frame(
    rbind(
      simpDiff=simpDiff(Y=Y,Tr=Tr),
      rebar=rebar(Y=Y,Tr=Tr,yhat=yhat),
      reloopOLS=loop(Y=Y,Tr=Tr,Z=cbind(yhat=yhat),pred=loop_ols,p=if(condNs) mean(Tr) else .5),
      reloopRF=loop(Y=Y,Tr=Tr,Z=cbind(yhat=yhat),pred=loop_rf,p=if(condNs) mean(Tr) else .5),
      loop=loop(Y=Y,Tr=Tr,Z=Z,p=if(condNs) mean(Tr) else .5),
      reloopPlus=loop(Y=Y,Tr=Tr,Z=Z,pred=if(fast) reloop else reloop_slow,yhat=yhat,p=if(condNs) mean(Tr) else .5)
    ))
  names(out)=c('Est','Var')
  out$method=rownames(out)
  out$ps=ps
  out$model=model
#  attr(out,"ps") <- ps
  out
}


allEstCE <- function(Y,Tr,Z,yhat,ps="ps",model="mod"){
  stopifnot(length(model)==1)
  stopifnot(length(ps)==1)
  out <- as.data.frame(
    rbind(
      simpDiff=simpDiff(Y=Y,Tr=Tr),
      rebar=rebar(Y=Y,Tr=Tr,yhat=yhat),
      #reloopOLS=unlist(ate.glmnet(X=cbind(yhat),Y=Y,W=Tr))[1:2],
      reloopRF=unlist(ate.randomForest(X=cbind(yhat,int=1),Y=Y,W=Tr))[1:2],
      loop=unlist(ate.randomForest(X=Z,Y=Y,W=Tr))[1:2],
      reloopPlus=unlist(ate.randomForest(X=cbind(Z,yhat),Y=Y,W=Tr))[1:2]#,
    ))
 # names(out)=c('Est','Var')
  out$method=rownames(out)
  out$ps=ps
  out$model=model
#  attr(out,"ps") <- ps
  out
}
