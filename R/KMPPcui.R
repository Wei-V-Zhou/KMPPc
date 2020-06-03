#' KMPPcui
#' 
#' Launch the KMPPc user interface in local machine
#'
#' This function will automatically launch the KMPPc user interface in a web browser. 
#' The user interface provides many powerful functions which is not available by command line programming.
#' It also provides a much easier and more convenient way to quickly explore TCGA data and construct kmplot you are interested in.
#' The user interface can also be accessed by https://zhouwei.shinyapps.io/KMPPc/. Neither R nor any packages are required in this online version.
#' However, it is highly recommended that the user interface be launched locally for faster running speed.
#' 
#' @export
#' @import shiny cgdsr ggplot2 stringr ggpubr survival survminer
#' @author Wei Zhou <247328181@@qq.com>
#' @examples
#' \dontrun{
#'    KMPPcui()
#' }

KMPPcui <- function() {
      shiny::runApp(system.file("shiny", package = "KMPPc"))
}
