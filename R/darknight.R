if (!requireNamespace("cBioPortalData", quietly = TRUE))
  BiocManager::install("cBioPortalData")

# setwd("C:/Users/asus1/AppData/Local/Temp/RtmpGsLI18/")
# filename <- dir()
# for (i in 1: length(filename)) {
#   install.packages(filename[i], repos = NULL, type = "win.binary")
# }

library("cBioPortalData")
cbio <- cBioPortal()

gbm <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "gbm_tcga",
                      genePanelId = "IMPACT341",
                      molecularProfileIds = c("gbm_tcga_rppa", "gbm_tcga_mrna")
)

laml <- cBioDataPack("laml_tcga")
