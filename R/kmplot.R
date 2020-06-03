rm(list = ls())
gc()

kmplot <- function(data = NULL, plotType = Typical) {
  
  ## 0. prepare environment and load libraries
  # rm(list = ls())
  # gc()
  set.seed(12345)
  graphics.off()
  options(stringsAsFactors = FALSE)
  # load packages
  pkgs <- c("ggplot2", "stringr", "ggpubr", "survival", "survminer", "export")
  # installpkgs <- function(pkgs){
  #   new.pkgs <- pkgs[!(pkgs %in% installed.packages()[ , "Package"])]
  #   if (length(new.pkgs))
  #    BiocManager::install(new.pkgs, ask = F, update = F)
  #   sapply(pkgs, require, character.only = T)
  # }
  # installpkgs(pkgs)
  lapply(pkgs, library, character.only = T)
  
  ## 3. K-M plot analysis
  if (is.null(data)) {
    load("kmplotdata.Rdata")
    dat2 <<- dat2
  } else {
    print("Please guarantee your file is dat2!")
    cat("***Notation: Or you will load the kmplotdata.Rdata on your own!")
  }
  
  # survival plot
  factor <- names(table(dat2[ , 2]))
  my.surv <<- Surv(dat2[ , 1], dat2[ , 2] == 'Recurred/Progressed')
  kmfit <<- survfit(my.surv~group, data = dat2)
  plot(kmfit, col = c("red", "blue"))
  img <- ggsurvplot(kmfit, palette = c("#E7B800","#2E9FDF"),
                    conf.int = TRUE, pval = TRUE, xlab = "Time / Month",
                    ggtheme = theme_light(), risk.table = TRUE, ncensor.plot = TRUE)
  print(img)
  label <- data.frame(kmfit$strata)
  labs <- cbind(rownames(label), label)
  colnames(labs) <- c("gene_group", "n")
  gene_high <- paste0(labs[1, 1] , "\n (n = ", labs[1, 2], ")")
  gene_low  <- paste0(labs[2, 1] , "\n (n = ", labs[2, 2], ")")
  names(kmfit$strata) <- c(gene_high, gene_low)
  g <- theme_survminer(font.main = c(16, "bold", "darkblue"),
                       font.x = c(14, "bold", "black"),
                       font.y = c(14, "bold", "black"),
                       font.tickslab = c(14, "bold", "black")) +
    theme(axis.line.x = element_line(size = 1),
          axis.line.y = element_line(size = 1),
          axis.ticks.x = element_line(size = 1),
          axis.ticks.y = element_line(size = 1),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12, face = "bold"))
  p <- ggsurvplot(kmfit, palette = c("red","blue"),
                  conf.int = F, pval = T, pval.method = T, surv.scale = "percent",
                  xlab ="Months" , ylab = "% Disease-free", 
                  ggtheme = g, risk.table = F, ncensor.plot = F) 
  
  print(p)
  
}
