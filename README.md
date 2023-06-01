# Replication code for Sales, Prihar, Gagnon-Bartch, Heffernan (2023). "Using Auxiliary Data to Boost Precision in the Analysis of A/B Tests on an Online Educational Platform: New Data and New
Results" forthcoming in the Journal of Educational Data Mining

This code replicates the effect and standard error estimation in the paper, but does not replicate the neural net model fit in the remnant. 
For raw data, neural net code (in Python), and the resulting predictions, please see https://osf.io/k8ph9/

The replication is inexact due to randomess in Random Forest algorithms. Results should be qualitatively similar to those reported in the paper. 

Follow these steps to replicate results:

1. Download data and Neural Net results from https://osf.io/k8ph9/:
    -  data/covariates.zip
    -  data/exp_norm_map.csv
    -  results/experiment_results.csv
    -  results/male_experiment_results.csv
    -  results/non_male_experiment_results.csv
2. Place _all_ of these files in the `data/` directory for this project
3. Un-zip covariates.zip. By this point, the project "data" folder should contain containing four .csv files and a sub-folder called "covariates"
4. Run the analysis from the command line:
    - sequentially: `Rscript code/replicate.r`
    - in parallel: `Rscript code/replicate.r [n.cores]` where `[n.cores]` is the number of cores you want to devote to the project.

Note: currently, all of the results use the `fast=TRUE` option of the `loop()` function. To use `fast=FALSE`, edit line 18 of [/code/replicate.r](/code/replicate.r)


```
R version 4.3.0 (2023-04-21)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 9.2 (Plow)

Matrix products: default
BLAS/LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggpmisc_0.5.2        ggpp_0.5.2           forcats_1.0.0       
 [4] estimatr_1.0.0       clubSandwich_0.5.8   ggeffects_1.2.2     
 [7] xtable_1.8-4         readr_2.1.4          ggplot2_3.4.2       
[10] purrr_0.3.5          dplyr_1.1.2          tikzDevice_0.12.4   
[13] loop.estimator_1.0.0

loaded via a namespace (and not attached):
 [1] sandwich_3.0-2     utf8_1.2.3         generics_0.1.3     lattice_0.21-8    
 [5] hms_1.1.3          digest_0.6.30      magrittr_2.0.3     grid_4.3.0        
 [9] filehash_2.4-5     jsonlite_1.8.3     Matrix_1.5-4       Formula_1.2-5     
[13] survival_3.5-5     gridExtra_2.3      fansi_1.0.4        scales_1.2.1      
[17] textshaping_0.3.6  cli_3.6.1          rlang_1.1.1        munsell_0.5.0     
[21] splines_4.3.0      withr_2.5.0        tools_4.3.0        SparseM_1.81      
[25] polynom_1.4-1      tzdb_0.4.0         MatrixModels_0.5-1 colorspace_2.1-0  
[29] vctrs_0.6.2        R6_2.5.1           zoo_1.8-12         lifecycle_1.0.3   
[33] MASS_7.3-58.4      insight_0.19.2     ragg_1.2.5         pkgconfig_2.0.3   
[37] pillar_1.9.0       gtable_0.3.3       glue_1.6.2         Rcpp_1.0.10       
[41] systemfonts_1.0.4  tibble_3.2.1       tidyselect_1.2.0   farver_2.1.1      
[45] labeling_0.4.2     compiler_4.3.0     quantreg_5.95     
```
