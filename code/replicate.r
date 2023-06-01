library(loop.estimator)
library(parallel)
library(tikzDevice)
library(dplyr)
library(purrr)
library(ggplot2)
library(readr)
library(xtable)
library(ggeffects)
library(clubSandwich)
library(estimatr)
library(forcats)
library(ggpmisc)

nclust <- commandArgs(TRUE)

full <- FALSE
fast <- TRUE

if(length(nclust)==0) nclust <- 0 else nclust <- as.numeric(nclust)

print(nclust)

if(nclust>0){
  if(.Platform$OS.type=='windows'){
    sock <- TRUE
    cl <- makeCluster(nclust, outfile="")
    ce <- clusterEvalQ(cl,source('code/reloopFunctions.r'))
    ce <- clusterEvalQ(cl,library(loop.estimator))
    ce <- clusterEvalQ(cl,library(dplyr))
    ce <- clusterEvalQ(cl,library(readr))
    ce <- clusterEvalQ(cl,library(purrr))
  } else sock <- FALSE
}


source('code/dirs.r')


source('code/reloopFunctions.r')

#############################
###  main analysis
#############################

### make pairwise data
## merge remnant imputations w treatment assignment, outcome, and covariate files
## mean-imputation for covariates
## format the data into pairwise comparisons
## create processedData/datPW.RData
runData <- TRUE
if(file.exists('processedData/pairwiseData.RData'))
  if(file.mtime('processedData/pairwiseData.RData')>file.mtime('code/data.r'))
    runData <- FALSE

if(runData | full){
  print('processing data')
  source('code/data.r')
} else print('skipping data')

### estimate effects & standard errors
print('estimation')
source('code/estimate.r')

print('plots')
source('code/plots.r')
