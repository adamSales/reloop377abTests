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


nclust <- commandArgs(TRUE)

if(length(nclust)==0) nclust <- 0 else nclust <- as.numeric(nclust)

if(nclust>0){
  cl <- makeCluster(nclust)
  ce <- clusterEvalQ(cl,source('code/reloopFunctions.r'))
  ce <- clusterEvalQ(cl,library(loop.estimator))
  ce <- clusterEvalQ(cl,library(dplyr))
  ce <- clusterEvalQ(cl,library(readr))
  ce <- clusterEvalQ(cl,library(purrr))

}



source('code/reloopFunctions.r')

#############################
###  main analysis
#############################

### make pairwise data
## merge remnant imputations w treatment assignment, outcome, and covariate files
## mean-imputation for covariates
## format the data into pairwise comparisons
## create processedData/datPW.RData
source('code/data.r')

### estimate effects & standard errors
source('code/estimate.r')

source('code/plots.r')
