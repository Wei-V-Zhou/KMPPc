
dataprocess <- function(){
  
  ## 0. prepare environment and load libraries
  # rm(list = ls())
  # gc()
  set.seed(12345)
  graphics.off()
  options(stringsAsFactors = FALSE)
  # load packages
  pkgs <- c("ggplot2", "stringr", "ggpubr", "ggstatsplot", "export")
  # installpkgs <- function(pkgs){
  #   new.pkgs <- pkgs[!(pkgs %in% installed.packages()[ , "Package"])]
  #   if (length(new.pkgs))
  #    BiocManager::install(new.pkgs, ask = F, update = F)
  #   sapply(pkgs, require, character.only = T)
  # }
  # installpkgs(pkgs)
  lapply(pkgs, library, character.only = T)
  
  ## 2. data preproceeding
  rm(list = ls())
  gc()
  options(stringsAsFactors = F)
  # load data
  load("survival_inputdata.Rdata")
  # view the clinical data
  clinicaldata_view <- as.matrix(colnames(myclinicaldata))
  # read the clinical information
  choose_columns = c("DFS_MONTHS", "DFS_STATUS")
  choose_clinicaldata = myclinicaldata[ , choose_columns]
  dat1 <- choose_clinicaldata[!is.na(choose_clinicaldata$DFS_MONTHS), ]
  dat2 <- cbind(dat1, exprSet[rownames(dat1), ])
  colnames(dat2)[3] <- "JAG2"
  write.csv(dat2, "JAG2_TCGA_expr_pancreas1.csv")
  # expressed genes plot
  p <- ggboxplot(dat2, x="DFS_STATUS", y="JAG2", color="DFS_STATUS", palette="jco", add="jitter")
  p + stat_compare_means(method = "t.test")
  dat2$JAG2_group = ifelse(dat2$JAG2 > median(dat2$JAG2), 'high', 'low')
  # dat2$JAG2_group = ifelse(dat2$JAG2 > quantile(dat2$JAG2)[4], 'high', 'low')
  ggbetweenstats(data = dat2, x = JAG2_group, y = JAG2)
  
}


## 3. K-M plot analysis
# survival plot
attach(dat2)
{
  table(JAG2_group)
  my.surv <- Surv(DFS_MONTHS, DFS_STATUS == 'Recurred/Progressed')
  kmfit <- survfit(my.surv~JAG2_group, data = dat2)
  plot(kmfit, col = c("red", "blue"))
  ggsurvplot(kmfit, palette=c("#E7B800","#2E9FDF"),
             conf.int = TRUE, pval = TRUE, xlab = "Time / Month",
             ggtheme = theme_light(), risk.table = TRUE, ncensor.plot = TRUE)
}
detach(dat2)


sessionInfo()
# R version 3.6.2 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.936 
# [2] LC_CTYPE=Chinese (Simplified)_China.936   
# [3] LC_MONETARY=Chinese (Simplified)_China.936
# [4] LC_NUMERIC=C                              
# [5] LC_TIME=Chinese (Simplified)_China.936    
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
# [1] ggstatsplot_0.1.3   survminer_0.4.6     survival_3.1-8     
# [4] ggpubr_0.2.4        magrittr_1.5        cgdsr_1.3.0        
# [7] export_0.2.2        gplots_3.0.1.1      readxl_1.3.1       
# [10] ggfortify_0.4.8     pheatmap_1.0.12     reshape2_1.4.3     
# [13] stringr_1.4.0       ggplot2_3.2.1       Biobase_2.44.0     
# [16] BiocGenerics_0.30.0
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3          SparseM_1.77             
# [3] R.methodsS3_1.7.1         coda_0.19-3              
# [5] tidyr_1.0.0               knitr_1.26               
# [7] multcomp_1.4-10           data.table_1.12.6        
# [9] inline_0.3.15             generics_0.0.2           
# [11] callr_3.3.2               cowplot_1.0.0            
# [13] TH.data_1.0-10            future_1.15.1            
# [15] webshot_0.5.2             xml2_1.2.2               
# [17] httpuv_1.5.2              ggsci_2.9                
# [19] StanHeaders_2.19.0        assertthat_0.2.1         
# [21] WRS2_1.0-0                xfun_0.11                
# [23] hms_0.5.2                 evaluate_0.14            
# [25] promises_1.1.0            DEoptimR_1.0-8           
# [27] caTools_1.17.1.2          km.ci_0.5-2              
# [29] htmlwidgets_1.5.1         mcmc_0.9-6               
# [31] reshape_0.8.8             stats4_3.6.2             
# [33] paletteer_0.2.1           purrr_0.3.3              
# [35] ellipsis_0.3.0            crosstalk_1.0.0          
# [37] rcompanion_2.3.7          dplyr_0.8.3              
# [39] backports_1.1.5           insight_0.7.0            
# [41] ggcorrplot_0.1.3          MCMCpack_1.4-4           
# [43] libcoin_1.0-5             jmvcore_1.0.8            
# [45] vctrs_0.2.0               quantreg_5.52            
# [47] sjlabelled_1.1.1          abind_1.4-5              
# [49] withr_2.1.2               metaBMA_0.6.2            
# [51] robustbase_0.93-5         emmeans_1.4.2            
# [53] prettyunits_1.0.2         mnormt_1.5-5             
# [55] cluster_2.1.0             lazyeval_0.2.2           
# [57] crayon_1.3.4              pkgconfig_2.0.3          
# [59] labeling_0.3              nlme_3.1-142             
# [61] statsExpressions_0.1.1    rlang_0.4.1              
# [63] globals_0.12.4            lifecycle_0.1.0          
# [65] miniUI_0.1.1.1            groupedstats_0.1.0       
# [67] skimr_2.0.1               LaplacesDemon_16.1.1     
# [69] MatrixModels_0.4-1        sandwich_2.5-1           
# [71] EMT_1.1                   modelr_0.1.5             
# [73] cellranger_1.1.0          matrixStats_0.54.0       
# [75] broomExtra_0.0.6          lmtest_0.9-37            
# [77] flextable_0.5.6           Matrix_1.2-18            
# [79] loo_2.1.0                 mc2d_0.1-18              
# [81] KMsurv_0.1-5              carData_3.0-3            
# [83] boot_1.3-23               zoo_1.8-6                
# [85] base64enc_0.1-3           processx_3.4.1           
# [87] rjson_0.2.20              parameters_0.3.0         
# [89] bitops_1.0-6              R.oo_1.23.0              
# [91] KernSmooth_2.23-16        ggExtra_0.9              
# [93] rgl_0.100.30              multcompView_0.1-7       
# [95] manipulateWidget_0.10.0   coin_1.3-1               
# [97] robust_0.4-18.1           ggsignif_0.6.0           
# [99] scales_1.1.0              plyr_1.8.4               
# [101] gdata_2.18.0              compiler_3.6.2           
# [103] rstantools_2.0.0          RColorBrewer_1.1-2       
# [105] lme4_1.1-21               rrcov_1.4-7              
# [107] cli_1.1.0                 listenv_0.7.0            
# [109] pbapply_1.4-2             ps_1.3.0                 
# [111] TMB_1.7.15                Brobdingnag_1.2-6        
# [113] MASS_7.3-51.4             mgcv_1.8-31              
# [115] tidyselect_0.2.5          stringi_1.4.3            
# [117] forcats_0.4.0             ggrepel_0.8.1            
# [119] bridgesampling_0.7-2      survMisc_0.5.5           
# [121] grid_3.6.2                tools_3.6.2              
# [123] rio_0.5.16                rvg_0.2.2                
# [125] rstudioapi_0.10           uuid_0.1-2               
# [127] foreign_0.8-72            gridExtra_2.3            
# [129] stargazer_5.2.2           pairwiseComparisons_0.1.2
# [131] farver_2.0.1              digest_0.6.22            
# [133] shiny_1.4.0               nortest_1.0-4            
# [135] jmv_1.0.8                 Rcpp_1.0.3               
# [137] car_3.0-5                 broom_0.5.2              
# [139] metafor_2.1-0             ez_4.4-0                 
# [141] BayesFactor_0.9.12-4.2    performance_0.4.0        
# [143] later_1.0.0               httr_1.4.1               
# [145] gdtools_0.2.1             psych_1.8.12             
# [147] sjstats_0.17.7            colorspace_1.4-1         
# [149] splines_3.6.2             expm_0.999-4             
# [151] systemfonts_0.1.1         fit.models_0.5-14        
# [153] xtable_1.8-4              jsonlite_1.6             
# [155] nloptr_1.2.1              rstan_2.19.2             
# [157] zeallot_0.1.0             modeltools_0.2-22        
# [159] R6_2.4.1                  broom.mixed_0.2.4        
# [161] pillar_1.4.2              htmltools_0.4.0          
# [163] mime_0.7                  glue_1.3.1               
# [165] fastmap_1.0.1             minqa_1.2.4              
# [167] codetools_0.2-16          pkgbuild_1.0.6           
# [169] pcaPP_1.9-73              mvtnorm_1.0-11           
# [171] furrr_0.1.0               lattice_0.20-38          
# [173] tibble_2.1.3              curl_4.2                 
# [175] DescTools_0.99.30         gtools_3.8.1             
# [177] officer_0.3.6             logspline_2.1.15         
# [179] zip_2.0.4                 openxlsx_4.1.3           
# [181] rmarkdown_1.17            repr_1.0.1               
# [183] munsell_0.5.0             rsample_0.0.5            
# [185] sjmisc_2.8.2              haven_2.2.0              
# [187] gtable_0.3.0              bayestestR_0.4.0 



#============================#
#       Musician: Resonance  #
#           Date: 2020/03/16 #
# Revised author: Resonance  #
#       1st Time: 2020/03/17 #
#       2nd Time: 2020/05/27 #
#============================#