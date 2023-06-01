load('processedData/pairwiseData.RData')


#### covariates
covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]
p=length(covNames)

S <- ifelse(fast,"Fast","Slow")


if(nclust>0) if(sock){
  clusterExport(cl,"datPW")
  clusterExport(cl,"covNames")
}


#################################################################
#### Main Analysis
#################################################################
LAP <- if(nclust>0){
  if(sock){
     function(X,FUN,...) parLapply(cl,X,FUN,...)
  } else {
    function(X,FUN,...) mclapply(X,FUN,mc.cores=nclust,mc.preschedule=FALSE,...)
  }
} else function(X,FUN,...) lapply(X,FUN,...)

estMain <- TRUE
if(file.exists(paste0('results/resTotal',S,'.Rdata')))
 if(file.mtime(paste0('results/resTotal',S,'.Rdata'))>
    max(
      file.mtime('code/estimateMain.r'),
      file.mtime('processedData/pairwiseData.RData'),
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE


source('code/estimateMain.r')


Contrasts=makeContrast(resTotal)

Contrasts$model=factor(Contrasts$model,levels=c('action','student','assignment','combined'))

save(Contrasts,file=paste0('results/contrasts',S,'.Rdata'))


#################################################################
#### Subgroup Analysis
#################################################################

estSub <- TRUE
if(file.exists(paste0('results/subgroupResults',S,'.Rdata')))
 if(file.mtime(paste0('results/subgroupResults',S,'.Rdata'))>
    max(
      file.mtime('code/estimateSubgroup.r'),
      file.mtime('processedData/pairwiseData.RData'),
      file.mtime('code/reloopFunctions.r')))
        estSub <- FALSE

source('code/estimateSubgroup.r')

resSub2=NULL
for(i in 1:10) for(j in 1:2) resSub2=c(resSub2,resSubgroups[[i]][[j]])

rtps=map_chr(resTotal,~.$ps[1])
resSub2 <- resSub2[map_lgl(resSub2,~.$ps[1]%in%rtps)]

resSub2 <- resSub2[vapply(resSub2,function(x) any(is.finite(x$Var))&min(x$Var,na.rm=TRUE)>1e-10,TRUE)]

save(resSub2,file=paste0('results/subgroupResults2',S,'.Rdata'))

ContrastsSub=makeContrast(resSub2) #map_dfr(resSub2,makeContrast)

save(ContrastsSub,file=paste0('results/contrastsSub',S,'.Rdata'))


#################################################################
#### Gender Experiment
#################################################################
load('processedData/pairwiseDataGender.RData')


estGen <- TRUE
if(file.exists(paste0('results/resGender',S,'.Rdata')))
 if(file.mtime(paste0('results/resGender',S,'.Rdata'))>
    max(
      file.mtime('code/estimateGender.r'),
      file.mtime('processedData/pairwiseDataGender.RData'),
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE

source('code/estimateGender.r')


ContrastsGender=makeContrast(resGender)

save(ContrastsGender,file=paste0('results/contrastsGender',S,'.Rdata'))



##################
### Post Stratification
#################

## from Ethan's email
inferredGender=c(
unknown = 24486,
male = 22397,
female = 18394,
mostly_male = 3138,
mostly_female = 1591,
andy = 1131)

pmale=sum(inferredGender[c('male','mostly_male')])/sum(inferredGender)

PS <- sapply(unique(ContrastsGender$ps),
             function(ps){
               resM <- resGender[[paste0('M.',ps)]]
               resF <- resGender[[paste0('F.',ps)]]
               Est=pmale*resM[,'Est']+(1-pmale)*resF[,'Est']
               Var=pmale^2*resM[,'Var']+(1-pmale)^2*resF[,'Var']
               cbind(Est=Est,Var=Var,resM[,c('method','ps')],model='PS')},
             simplify=FALSE)

ContrastsPS <- makeContrast(PS)

save(PS,ContrastsPS,file=paste0('results/PostStratification',S,'.Rdata'))
