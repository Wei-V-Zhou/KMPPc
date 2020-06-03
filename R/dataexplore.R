#' dataexplore
#' 
#' dataexplore for TCGA
#'
#' By default, this function first input cancerType such as "BreastCancer" to explore the dataset whatever you like.
#'
#' @param cancerType The cancertype you'd like to explore, including "BreastCancer", "PancreaticCancer", "Glioma"
#' @param studyId This argument you can explore by yourself, because there is too many researches or database you can take advantage of.
#' @param dataType It provides "mrna", "CNA", "mutation" for more information to dig.
#' @return There is a RData file (survival_inputdata.Rdata) exported to the data folder.
#' @export
#' @import cgdsr
#' @author Wei Zhou <247328181@@qq.com>
#' @examples
#' ## Explore the interested dataset by input the specified cancertype
#' dataexplore(cancerType = "BreastCancer", dataType = "mrna")
#' 
#' ## If you have interested datatype, you can input, such as mrna, CNA, mutation, etc
#' dataexplore(cancerType = "BreastCancer", studyId = NULL, dataType = "mrna")

dataexplore <- function(cancerType, studyId = NULL, dataType = NULL) {

  ## 0. prepare environment and load libraries
  # rm(list = ls())
  # gc()
  # set.seed(12345)
  # graphics.off()
  # options(stringsAsFactors = FALSE)
  # load packages
  # pkgs <- c("cgdsr", "ggplot2", "stringr", "ggpubr", "survival",
  #           "survminer", "readxl", "ggstatsplot", "export")
  # installpkgs <- function(pkgs){
  #   new.pkgs <- pkgs[!(pkgs %in% installed.packages()[ , "Package"])]
  #   if (length(new.pkgs))
  #    BiocManager::install(new.pkgs, ask = F, update = F)
  #   sapply(pkgs, require, character.only = T)
  # }
  # installpkgs(pkgs)
  # lapply(pkgs, library, character.only = T)

  ## 1. connect the TCGA database and obtain the data
  # create a CGDS connection objext
  mycgds <- CGDS("http://www.cbioportal.org/")
  # set verbose options to debug and troubleshoot issues
  setVerbose(mycgds, TRUE)
  # guarantee the unauthorized accessment
  mysecurecgds <- CGDS("http://cbioportal.mskcc.org/",
                       token = "fd0522cb-7972-40d0-9d83-cb4c14e8a337")
  # get list of cancer studies at server for view to choose
  cancerstudy <- getCancerStudies(mycgds)
  if(cancerType == "BreastCancer"){
    cancerid <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "brca"))), 1]
    studydoi <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "brca"))), 2]
  } else if (cancerType == "PancreaticCancer") {
    cancerid <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "paad"))), 1]
    studydoi <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "paad"))), 2]
  } else if (cancerType == "Glioma") {
    cancerid <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "glioma|lgg"))), 1]
    studydoi <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "glioma|lgg"))), 2]
  }

  # get available case lists for a given cancer study
  if(is.null(studyId)){
    MainstrId <- cbind(cancerid, studydoi)
    print.table(studydoi, right = F, justify = "centre")
    if(interactive()){
      repeat{
        ANSWER <- readline("Please input the dataset number: ")
        num <- as.numeric(ANSWER)
        if(is.na(num) || num <= 0 || num > length(studydoi)){
          print("Please type right number format!")
        } else if (num > 0 && num <= length(studydoi)) {
          break;
        }
      }
    }
    MainstrDat <- cancerid[num]
  } else if (is.numeric(studyId)) {
    MainstrDat <- cancerid[as.numeric(studyId)]
  } else {
    stop("Please input right format studyId!")
  }

  # get available genetic profiles
  mygeneticprofile <- getGeneticProfiles(mycgds, MainstrDat)
  if (dataType == "mrna") {
    # mrna data
    mrnaFile <- grep("mrna", mygeneticprofile[ , 1])
    print.table(mygeneticprofile[mrnaFile, 3], right = F, justify = "left")
    repeat{
      ANSWER <- readline("Please input your interested genetic file number: ")
      mrnaNum <- as.numeric(ANSWER)
      if (is.na(mrnaNum) || mrnaNum <= 0 || mrnaNum > length(mrnaFile)) {
        print("Please type right number format!")
      } else if (mrnaNum > 0 && mrnaNum <= length(mrnaFile)) {
        break;
      }
    }
    dataTypeFile = mrnaFile[mrnaNum]
  } else if (dataType == "methylation") {
    # methylation data
    methylationFile  <- grep("methylation", mygeneticprofile[ , 1])
    print.table(mygeneticprofile[methylationFile, 3], right = F, justify = "left")
    repeat{
      ANSWER <- readline("Please input your interested genetic file number: ")
      methylationNum <- as.numeric(ANSWER)
      if (is.na(methylationNum) || methylationNum <= 0 || methylationNum > length(methylationFile)) {
        print("Please type right number format!")
      } else if (methylationNum > 0 && methylationNum <= length(methylationFile)) {
        break;
      }
    }
    dataTypeFile = methylationFile[methylationNum]
  } else if (dataType == "CNA") {
    # CNA data
    cnaFile  <- grep("CNA", mygeneticprofile[ , 1])
    print.table(mygeneticprofile[cnaFile, 3], right = F, justify = "left")
    repeat{
      ANSWER <- readline("Please input your interested genetic file number: ")
      cnaNum <- as.numeric(ANSWER)
      if (is.na(cnaNum) || cnaNum <= 0 || cnaNum > length(cnaFile)) {
        print("Please type right number format!")
      } else if (cnaNum > 0 && cnaNum <= length(cnaFile)) {
        break;
      }
    }
    dataTypeFile = cnaFile[cnaNum]
  } else if (dataType == "mutation") {
    # mutation
    mutationFile <- grep("mutation", mygeneticprofile[ , 1])
    print.table(mygeneticprofile[mutationFile, 3], right = F, justify = "left")
    repeat{
      ANSWER <- readline("Please input your interested number: ")
      mutationNum <- as.numeric(ANSWER)
      if (is.na(mutationNum) || mutationNum <= 0 || mutationNum > length(mutationFile)) {
        print("Please type right number format!")
      } else if (mutationNum > 0 && mutationNum <= length(mutationFile)) {
        break;
      }
    }
    dataTypeFile = mutationFile[mutationNum]
  } else {
    stop("Please input right format dataType!")
  }
  geneticFile <- getGeneticProfiles(mycgds, MainstrDat)[dataTypeFile, 1]

  # get available genecase list
  mycaselist <- getCaseLists(mycgds, MainstrDat)
  if (dataType == "mrna") {
    # mRNA data
    mrnaList  <- grep("mrna", mycaselist[ , 1])
    print.table(mycaselist[mrnaList, 3], right = F, justify = "centre")
    repeat{
      ANSWER <- readline("Please input your interested case list number: ")
      mrnaNum <- as.numeric(ANSWER)
      if (is.na(mrnaNum) || mrnaNum <= 0 || mrnaNum > length(mrnaList)) {
        print("Please type right number format!")
      } else if (mrnaNum > 0 && mrnaNum <= length(mrnaList)) {
        break;
      }
    }
    dataTypeList = mrnaList[mrnaNum]
  } else if (dataType == "methylation") {
    # methylation data
    methylationList  <- grep("methylation", mycaselist[ , 1])
    print.table(mycaselist[methylationList, 3], right = F, justify = "centre")
    repeat{
      ANSWER <- readline("Please input your interested case list number: ")
      methylationNum <- as.numeric(ANSWER)
      if (is.na(methylationNum) || methylationNum <= 0 || methylationNum > length(methylationList)) {
        print("Please type right number format!")
      } else if (methylationNum > 0 && methylationNum <= length(methylationList)) {
        break;
      }
    }
    dataTypeList = methylationList[methylationNum]
  } else if (dataType == "CNA") {
    # CNA data
    cnaList  <- grep("CNA", mycaselist[ , 2])
    print.table(mycaselist[cnaList, 3], right = F, justify = "centre")
    repeat{
      ANSWER <- readline("Please input your interested case list number: ")
      cnaNum <- as.numeric(ANSWER)
      if (is.na(cnaNum) || cnaNum <= 0 || cnaNum > length(cnaList)) {
        print("Please type right number format!")
      } else if (cnaNum > 0 && cnaNum <= length(cnaList)) {
        break;
      }
    }
    dataTypeList = cnaList[cnaNum]
  } else if (dataType == "mutation") {
    # mutation
    mutationList <- grep("mutation", mycaselist[ , 2])
    print.table(mycaselist[mutationList, 3], right = F, justify = "centre")
    repeat{
      ANSWER <- readline("Please input your interested case list number: ")
      mutationNum <- as.numeric(ANSWER)
      if (is.na(mutationNum) || mutationNum <= 0 || mutationNum > length(mutationList)) {
        print("Please type right number format!")
      } else if (mutationNum > 0 && mutationNum <= length(mutationList)) {
        break;
      }
    }
    dataTypeList = mutationList[mutationNum]
  } else {
    stop("Please input right format dataType!")
  }
  caseList <- getCaseLists(mycgds, MainstrDat)[dataTypeList, 1]

  # get data slices for a specified list of genes, genetic profiles and case list
  repeat{
    ANSWER <- readline("Please input your interested gene name: ")
    geneName <- toupper(ANSWER)
    if (geneName == "") {
      print("Please type right number format!")
    } else {
      break;
    }
  }
  exprSet = getProfileData(mycgds, geneName, geneticFile, caseList)

  # get clinical data for the case list
  myclinicaldata <- getClinicalData(mycgds, caseList)
  # save the data
  save(exprSet, myclinicaldata, file="../data/survival_inputdata.Rdata")
}


# sessionInfo()
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
# [1] export_0.2.2      ggstatsplot_0.1.3 readxl_1.3.1      survminer_0.4.6   survival_3.1-8
# [6] ggpubr_0.2.4      magrittr_1.5      stringr_1.4.0     ggplot2_3.2.1     cgdsr_1.3.0
#
# loaded via a namespace (and not attached):
# [1] tidyselect_0.2.5          lme4_1.1-21               robust_0.4-18.1
# [4] htmlwidgets_1.5.1         grid_3.6.3                munsell_0.5.0
# [7] codetools_0.2-16          future_1.15.1             miniUI_0.1.1.1
# [10] withr_2.1.2               Brobdingnag_1.2-6         metaBMA_0.6.2
# [13] colorspace_1.4-1          uuid_0.1-2                knitr_1.26
# [16] rstudioapi_0.10           stats4_3.6.3              DescTools_0.99.30
# [19] robustbase_0.93-5         officer_0.3.6             ggsignif_0.6.0
# [22] rcompanion_2.3.7          listenv_0.7.0             emmeans_1.4.2
# [25] rstan_2.19.2              repr_1.0.1                mnormt_1.5-5
# [28] KMsurv_0.1-5              MCMCpack_1.4-4            bridgesampling_0.7-2
# [31] coda_0.19-3               vctrs_0.2.0               generics_0.0.2
# [34] TH.data_1.0-10            metafor_2.1-0             xfun_0.11
# [37] R6_2.4.1                  BayesFactor_0.9.12-4.2    manipulateWidget_0.10.0
# [40] reshape_0.8.8             logspline_2.1.15          assertthat_0.2.1
# [43] promises_1.1.0            scales_1.1.0              multcomp_1.4-10
# [46] ggExtra_0.9               gtable_0.3.0              multcompView_0.1-7
# [49] globals_0.12.4            processx_3.4.1            mcmc_0.9-6
# [52] sandwich_2.5-1            rlang_0.4.1               MatrixModels_0.4-1
# [55] EMT_1.1                   zeallot_0.1.0             systemfonts_0.1.1
# [58] splines_3.6.3             TMB_1.7.16                lazyeval_0.2.2
# [61] broom_0.5.2               inline_0.3.15             rgl_0.100.30
# [64] reshape2_1.4.3            abind_1.4-5               modelr_0.1.5
# [67] crosstalk_1.0.0           backports_1.1.5           httpuv_1.5.2
# [70] rsconnect_0.8.16          tools_3.6.3               psych_1.8.12
# [73] ellipsis_0.3.0            stargazer_5.2.2           WRS2_1.0-0
# [76] ez_4.4-0                  Rcpp_1.0.3                plyr_1.8.4
# [79] base64enc_0.1-3           jmvcore_1.0.8             purrr_0.3.3
# [82] ps_1.3.0                  prettyunits_1.0.2         pbapply_1.4-2
# [85] cowplot_1.0.0             zoo_1.8-6                 LaplacesDemon_16.1.1
# [88] haven_2.2.0               ggrepel_0.8.1             cluster_2.1.0
# [91] furrr_0.1.0               data.table_1.12.6         openxlsx_4.1.3
# [94] flextable_0.5.6           SparseM_1.77              lmtest_0.9-37
# [97] mvtnorm_1.0-11            broomExtra_0.0.6          sjmisc_2.8.2
# [100] matrixStats_0.54.0        evaluate_0.14             hms_0.5.2
# [103] mime_0.9                  xtable_1.8-4              rio_0.5.16
# [106] sjstats_0.17.7            pairwiseComparisons_0.1.2 broom.mixed_0.2.4
# [109] gridExtra_2.3             rstantools_2.0.0          compiler_3.6.3
# [112] tibble_2.1.3              crayon_1.3.4              minqa_1.2.4
# [115] R.oo_1.23.0               StanHeaders_2.19.0        htmltools_0.4.0
# [118] mgcv_1.8-31               mc2d_0.1-18               pcaPP_1.9-73
# [121] later_1.0.0               libcoin_1.0-5             tidyr_1.0.0
# [124] rrcov_1.4-7               expm_0.999-4              sjlabelled_1.1.1
# [127] jmv_1.0.8                 MASS_7.3-51.4             boot_1.3-23
# [130] Matrix_1.2-18             car_3.0-5                 cli_1.1.0
# [133] R.methodsS3_1.7.1         parallel_3.6.3            insight_0.7.0
# [136] forcats_0.4.0             pkgconfig_2.0.3           km.ci_0.5-2
# [139] fit.models_0.5-14         statsExpressions_0.1.1    coin_1.3-1
# [142] foreign_0.8-72            skimr_2.0.1               xml2_1.2.2
# [145] paletteer_0.2.1           ggcorrplot_0.1.3          webshot_0.5.2
# [148] rvg_0.2.2                 estimability_1.3          callr_3.3.2
# [151] digest_0.6.25             parameters_0.3.0          rmarkdown_1.17
# [154] cellranger_1.1.0          survMisc_0.5.5            nortest_1.0-4
# [157] gdtools_0.2.1             curl_4.2                  modeltools_0.2-22
# [160] shiny_1.4.0.2             gtools_3.8.1              quantreg_5.52
# [163] rjson_0.2.20              nloptr_1.2.1              lifecycle_0.1.0
# [166] nlme_3.1-142              jsonlite_1.6.1            carData_3.0-3
# [169] groupedstats_0.1.0        pillar_1.4.2              lattice_0.20-38
# [172] loo_2.1.0                 fastmap_1.0.1             DEoptimR_1.0-8
# [175] pkgbuild_1.0.6            glue_1.3.1                bayestestR_0.4.0
# [178] zip_2.0.4                 stringi_1.4.3             performance_0.4.0
# [181] rsample_0.0.5             dplyr_0.8.3



#============================#
#       Musician: Resonance  #
#           Date: 2020/03/16 #
# Revised author: Resonance  #
#       1st Time: 2020/03/17 #
#       2nd Time: 2020/05/27 #
#============================#
