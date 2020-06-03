if (!requireNamespace("cBioPortalData", quietly = TRUE))
  BiocManager::install("cBioPortalData")

# setwd("C:/Users/asus1/AppData/Local/Temp/RtmpGsLI18/")
# filename <- dir()
# for (i in 1: length(filename)) {
#   install.packages(filename[i], repos = NULL, type = "win.binary")
# }

library("cBioPortalData")
library("rapiclient")

cbio <- cBioPortal()
client <- get_api(url = "https://pedcbioportal.kidsfirstdrc.org")

dipg <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "dipg",
                      genePanelId = "IMPACT341")


dipg <- cBioDataPack("dipg_cbttc")
