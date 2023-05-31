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
      reloopOLS=loop(Y=Y,Tr=Tr,Z=cbind(yhat=yhat),pred=loop_ols,p= .5),
      reloopRF=loop(Y=Y,Tr=Tr,Z=cbind(yhat=yhat),pred=loop_rf,p= .5),
      loop=loop(Y=Y,Tr=Tr,Z=Z,p= .5),
      reloopPlus=loop(Y=Y,Tr=Tr,Z=Z,pred=if(fast) reloop else reloop_slow,yhat=yhat,p= .5)
    ))
  names(out)=c('Est','Var')
  out$method=rownames(out)
  out$ps=ps
  out$model=model
  cat('.')
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


### this is a function for comparing SEs between estimates
makeContrast=function(res,ols=TRUE){
  rl=if(ols) 'reloopOLS' else 'reloopRF'
  v=if('Var'%in%colnames(res[[1]])) 'Var' else 'var'
  cont1=sapply(res,
               function(x) if(nrow(x)<3) NA else x['simpDiff',v]/x[rl,v])

  cont2=sapply(res,
               function(x) if(nrow(x)<3) NA else x['simpDiff',v]/x['reloopPlus',v])

  cont3=sapply(res,
               function(x) if(nrow(x)<3) NA else x['loop',v]/x['reloopPlus',v])

  ssMult=c(cont1,cont2,cont3)


  out <- data.frame(
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

  if('covName'%in%colnames(res[[1]])){
    id=map_dfr(res,~.[1,c('ps','covName','side','nt','nc')])
    id$PS=map_chr(strsplit(id$ps,"(Control)|(Treatment)"),~.[1])
    id$num=1:length(res)
    out <- cbind(rbind(id,id,id),out)
  } else{
    ps=sapply(res,function(x) x$ps[1])
    out <- cbind(ps=rep(ps,3),out)
  }

  out

}
