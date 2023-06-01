load('processedData/pairwiseData.RData')


#### covariates
covNames <- names(datPW)[startsWith(names(datPW),'student_prior')]
p=length(covNames)




if(nclust>0 & sock){
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
if(file.exists('results/resTotalSlow.RData'))
 if(file.mtime('results/resTotalSlow.RData')> 
    max(
      file.mtime('code/estimateMain.r'),
      file.mtime('processedData/pairwiseData.RData'),
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE


source('code/estimateMain.r')


Contrasts=map_dfr(resTotal,makeContrast)

Contrasts$model=factor(Contrasts$model,levels=c('action','student','assignment','combined'))

save(Contrasts,file='results/contrasts.RData')


#################################################################
#### Subgroup Analysis
#################################################################

estSub <- TRUE
if(file.exists('results/subgroupResults.RData'))
 if(file.mtime('results/subgroupResults.RData')> 
    max(
      file.mtime('code/estimateSubgroup.r'),
      file.mtime('processedData/pairwiseData.RData'),
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE

source('code/estimateSubgroup.r')

resSub2=NULL
for(i in 1:10) for(j in 1:2) resSub2=c(resSub2,resSubgroups[[i]][[j]])

resSub2 <- resSub2[map_lgl(resSub2,~.$ps[1]%in%names(resTotal[[1]]))]

resSub2 <- resSub2[vapply(resSub2,function(x) any(is.finite(x$Var))&min(x$Var,na.rm=TRUE)>1e-10,TRUE)]

save(resSub2,file='results/subgroupResults2.RData')

ContrastsSub=makeContrast(resSub2) #map_dfr(resSub2,makeContrast)

save(ContrastsSub,file='results/contrastsSub.RData')


#################################################################
#### Gender Experiment
#################################################################

estGen <- TRUE
if(file.exists('results/resGender.RData'))
 if(file.mtime('results/resGender.RData')> 
    max(
      file.mtime('code/estimateGender.r'),
      file.mtime('processedData/pairwiseDataGender.RData'),
      file.mtime('code/reloopFunctions.r')))
        estMain <- FALSE

source('code/estimateGender.r')


ContrastsGender=makeContrast(resGender)

save(ContrastsGender,file='results/contrastsGender.RData')



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

PS <- sapply(unique(ContrastsGender$problem_set),
             function(ps){
               resM <- resGender[[paste0('M.',ps)]]
               resF <- resGender[[paste0('F.',ps)]]
               Est=pmale*resM[,'Est']+(1-pmale)*resF[,'Est']
               Var=pmale^2*resM[,'Var']+(1-pmale)^2*resF[,'Var']
               cbind(Est=Est,Var=Var,resM[,c('method','ps')],model='PS')},
             simplify=FALSE)

ContrastsPS <- makeContrast(PS)

save(PS,ContrastsPS,file='results/PostStratification.RData')