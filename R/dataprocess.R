
dataprocess <- function(clinicaldata = NULL, exprSet = NULL, x_axis = NULL, y_axis = NULL) {
  
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
  # load data
  if (is.null(clinicaldata) || is.null(exprSet)) {
    load("../data/survival_inputdata.Rdata")
  } else {
    print("Please guarantee your two files are choose_clinicaldata and exprSet, respectively!")
    cat("***Notation: Or you will load the survival_inputdata.Rdata on your own!")
  }
  
  # view the clinical data
  clinicaldata_view <- as.matrix(colnames(myclinicaldata))
  choose_columns <- c(clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DFS_MONTHS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DFS_STATUS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DSS_MONTHS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DSS_STATUS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "OS_MONTHS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "OS_STATUS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "PFS_MONTHS")))],
                      clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "PFS_STATUS")))])
  cat("This clinical dataset only has the following data: \n***Notation: x_axis presents time, y_axis presents status; please choose the corresponding data!\n")
  for (i in 1:length(choose_columns)) {
    cat(c(i, ": ", choose_columns[i], "\n"))
  }
  
  # read the clinical information
  # choose x_axis
  if (is.null(x_axis)) {
    if(interactive()){
      repeat{
        ANSWER <- readline("Please input the x_axis number (time): ")
        num <- as.numeric(ANSWER)
        if(is.na(num) || num <= 0 || num > length(choose_columns)){
          print("Please type right number format!")
        } else if (num > 0 && num <= length(choose_columns)) {
          break;
        }
      }
      x_axis <- choose_columns[num]
    }
  } else {
    if (! x_axis %in% choose_columns) {
      cat("Warning! The input of x_axis is not in the clinical dataset!\n***Please type right format!\n")
      x_axis = NULL
      stop("Input again! Or you can eliminate x_axis and follow the tips.")
    }
  }
  # choose y_axis
  if (is.null(y_axis)) {
    if(interactive()){
      repeat{
        ANSWER <- readline("Please input the y_axis number (status): ")
        num <- as.numeric(ANSWER)
        if(is.na(num) || num <= 0 || num > length(choose_columns)){
          print("Please type right number format!")
        } else if (num > 0 && num <= length(choose_columns)) {
          if (x_axis != choose_columns[num]) {
            break;
          } else {
            print("y_axis can't be the same with x_axis!")
          }
        }
      }
      y_axis <- choose_columns[num]
    }
  } else if (x_axis == y_axis) {
    stop("x_axis can't be the same with y_axis!")
  } else {
    if (! y_axis %in% choose_columns) {
      cat("Warning! The input of y_axis is not in the clinical dataset!\n***Please type right format!\n")
      x_axis = NULL
      stop("Input again! Or you can eliminate y_axis and follow the tips.")
    }
  }
  
  # filter the data
  choose_column = c(x_axis, y_axis)
  choose_clinicaldata = myclinicaldata[ , choose_column]
  dat1 <- choose_clinicaldata[!is.na(choose_clinicaldata[ , 1]), ]
  dat2 <- cbind(dat1, exprSet[rownames(dat1), ])
  geneName <- names(exprSet)
  colnames(dat2)[3] <- geneName
  # save the filtered data
  if(interactive()){
    repeat{
      ANSWER <- readline("Save the clinical data and related exprset data[y/n]: ")
      num <- trimws(tolower(ANSWER), which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      if (num == "y" || num == "yes") {
        print("Save successfully!")
        filename = paste0(geneName, "_clinical_expr_data_", Sys.Date(), ".csv")
        write.csv(dat2, file = filename)
        break;
      } else if (num == "n" || num == "no") {
        break;
      } else {
        print("Wrong format! Please type y or n!")
      }
    }
  }
  
  # expressed genes plot
  p <- ggboxplot(dat2, x = y_axis, y = geneName, color = y_axis, palette = "jco", add = "jitter")
  p + stat_compare_means(method = "t.test")
  group <- paste0(geneName, "_group")
  dat2$group = ifelse(dat2[ , geneName] > median(dat2[ , names(exprSet)]), 'high', 'low')
  # dat2$group = ifelse(dat2[ , names(exprSet)] > quantile(dat2[ , names(exprSet)])[4], 'high', 'low')
  # save the data
  save(dat2, file = "../data/kmplotdata.Rdata")
  # plot the expression img
  dat2$geneName = dat2$HTRA1
  ggbetweenstats(data = dat2, x = group, y = geneName, xlab = "Patient group", ylab = paste0(geneName, "_expression"))
  # save the image
  # if(interactive()){
  #   repeat{
  #     ANSWER <- readline("Save the expression picture?[y/n]: ")
  #     num <- trimws(tolower(ANSWER), which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  #     if (num == "y" || num == "yes") {
  #       print("Save successfully!")
  #       filename = paste0(geneName, "_clinical_expr_data_", Sys.Date(), ".pptx")
  #       graph2ppt(file = filename, width = 7, height = 5)
  #       break;
  #     } else if (num == "n" || num == "no") {
  #       break;
  #     } else {
  #       print("Wrong format! Please type y or n!")
  #     }
  #   }
  # }
  
}

sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
# [3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
# [5] LC_TIME=Chinese (Simplified)_China.936    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] export_0.2.2      ggstatsplot_0.1.3 ggpubr_0.2.4      magrittr_1.5      stringr_1.4.0    
# [6] ggplot2_3.2.1    
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_0.2.5          lme4_1.1-21               robust_0.4-18.1          
# [4] htmlwidgets_1.5.1         grid_3.6.3                munsell_0.5.0            
# [7] codetools_0.2-16          future_1.15.1             miniUI_0.1.1.1           
# [10] withr_2.1.2               Brobdingnag_1.2-6         metaBMA_0.6.2            
# [13] colorspace_1.4-1          uuid_0.1-2                knitr_1.26               
# [16] rstudioapi_0.10           stats4_3.6.3              DescTools_0.99.30        
# [19] robustbase_0.93-5         officer_0.3.6             ggsignif_0.6.0           
# [22] rcompanion_2.3.7          listenv_0.7.0             labeling_0.3             
# [25] emmeans_1.4.2             rstan_2.19.2              repr_1.0.1               
# [28] mnormt_1.5-5              MCMCpack_1.4-4            farver_2.0.1             
# [31] bridgesampling_0.7-2      coda_0.19-3               vctrs_0.2.0              
# [34] generics_0.0.2            TH.data_1.0-10            metafor_2.1-0            
# [37] xfun_0.11                 R6_2.4.1                  BayesFactor_0.9.12-4.2   
# [40] manipulateWidget_0.10.0   reshape_0.8.8             logspline_2.1.15         
# [43] assertthat_0.2.1          promises_1.1.0            scales_1.1.0             
# [46] multcomp_1.4-10           ggExtra_0.9               gtable_0.3.0             
# [49] multcompView_0.1-7        globals_0.12.4            processx_3.4.1           
# [52] mcmc_0.9-6                sandwich_2.5-1            rlang_0.4.1              
# [55] MatrixModels_0.4-1        EMT_1.1                   zeallot_0.1.0            
# [58] systemfonts_0.1.1         splines_3.6.3             TMB_1.7.16               
# [61] lazyeval_0.2.2            broom_0.5.2               inline_0.3.15            
# [64] rgl_0.100.30              reshape2_1.4.3            abind_1.4-5              
# [67] modelr_0.1.5              crosstalk_1.0.0           backports_1.1.5          
# [70] httpuv_1.5.2              rsconnect_0.8.16          tools_3.6.3              
# [73] psych_1.8.12              ellipsis_0.3.0            stargazer_5.2.2          
# [76] WRS2_1.0-0                ez_4.4-0                  Rcpp_1.0.3               
# [79] plyr_1.8.4                base64enc_0.1-3           jmvcore_1.0.8            
# [82] purrr_0.3.3               ps_1.3.0                  prettyunits_1.0.2        
# [85] pbapply_1.4-2             cowplot_1.0.0             zoo_1.8-6                
# [88] LaplacesDemon_16.1.1      haven_2.2.0               ggrepel_0.8.1            
# [91] cluster_2.1.0             furrr_0.1.0               data.table_1.12.6        
# [94] openxlsx_4.1.3            flextable_0.5.6           SparseM_1.77             
# [97] lmtest_0.9-37             mvtnorm_1.0-11            broomExtra_0.0.6         
# [100] sjmisc_2.8.2              matrixStats_0.54.0        evaluate_0.14            
# [103] hms_0.5.2                 mime_0.9                  xtable_1.8-4             
# [106] rio_0.5.16                sjstats_0.17.7            pairwiseComparisons_0.1.2
# [109] broom.mixed_0.2.4         readxl_1.3.1              gridExtra_2.3            
# [112] rstantools_2.0.0          compiler_3.6.3            tibble_2.1.3             
# [115] crayon_1.3.4              minqa_1.2.4               StanHeaders_2.19.0       
# [118] htmltools_0.4.0           mgcv_1.8-31               mc2d_0.1-18              
# [121] pcaPP_1.9-73              later_1.0.0               tidyr_1.0.0              
# [124] libcoin_1.0-5             rrcov_1.4-7               expm_0.999-4             
# [127] sjlabelled_1.1.1          jmv_1.0.8                 MASS_7.3-51.4            
# [130] boot_1.3-23               Matrix_1.2-18             car_3.0-5                
# [133] cli_1.1.0                 parallel_3.6.3            insight_0.7.0            
# [136] forcats_0.4.0             pkgconfig_2.0.3           fit.models_0.5-14        
# [139] statsExpressions_0.1.1    coin_1.3-1                foreign_0.8-72           
# [142] skimr_2.0.1               xml2_1.2.2                paletteer_0.2.1          
# [145] ggcorrplot_0.1.3          webshot_0.5.2             rvg_0.2.2                
# [148] estimability_1.3          callr_3.3.2               digest_0.6.25            
# [151] parameters_0.3.0          rmarkdown_1.17            cellranger_1.1.0         
# [154] nortest_1.0-4             gdtools_0.2.1             curl_4.2                 
# [157] shiny_1.4.0.2             gtools_3.8.1              quantreg_5.52            
# [160] modeltools_0.2-22         rjson_0.2.20              nloptr_1.2.1             
# [163] lifecycle_0.1.0           nlme_3.1-142              jsonlite_1.6.1           
# [166] carData_3.0-3             groupedstats_0.1.0        pillar_1.4.2             
# [169] ggsci_2.9                 lattice_0.20-38           loo_2.1.0                
# [172] fastmap_1.0.1             DEoptimR_1.0-8            pkgbuild_1.0.6           
# [175] survival_3.1-8            glue_1.3.1                bayestestR_0.4.0         
# [178] zip_2.0.4                 stringi_1.4.3             performance_0.4.0        
# [181] rsample_0.0.5             dplyr_0.8.3



#============================#
#       Musician: Resonance  #
#           Date: 2020/03/16 #
# Revised author: Resonance  #
#       1st Time: 2020/03/17 #
#       2nd Time: 2020/05/27 #
#       3rd Time: 2020/06/03 #
#============================#