
kmplot <- function(dat2) {
  
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
  
  ## 3. K-M plot analysis
  load("kmplotdata.RData")
  # survival plot
  attach(dat2)
  {
    table(group)
    my.surv <- Surv(dat2[ , 1], dat2[ , 2] == 'Recurred/Progressed')
    kmfit <- survfit(my.surv~JAG2_group, data = dat2)
    plot(kmfit, col = c("red", "blue"))
    ggsurvplot(kmfit, palette=c("#E7B800","#2E9FDF"),
               conf.int = TRUE, pval = TRUE, xlab = "Time / Month",
               ggtheme = theme_light(), risk.table = TRUE, ncensor.plot = TRUE)
  }
  detach(dat2)
  
}