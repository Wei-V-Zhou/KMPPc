######################################################
##                      KMPPc                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##               Author:Wei Zhou, xxx               ##
##       Maintainer:Wei Zhou (247328181@qq.com)     ##
######################################################

library("shiny")

# Define UI for K-M Plot viewer application
shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("KMPPc: K-M Plot of Pan-cancer"),
    
    # Sidebar with controls to select a dataset and processing
    sidebarPanel(
        
        # Define main menu
        wellPanel(
            helpText(a("Github help docs", href = "https://github.com/Wei-V-Zhou/KMPPc", target = "_blank")),
            radioButtons("MainMenu", "Main Menu",
                         list("Reading in dataset" = "Input",
                              "Preprocessing" = "Preprocess",
                              "K-M Plot" = "kmplot",
                              "Differential gene analysis (optional)" = "Difftest",
                              "Filtering targets (optional)" = "Filter"
                              ))
            ),
                
        
        # Define secondary menu
        conditionalPanel(condition = "input.MainMenu == 'Input'",
                         wellPanel(
                             selectInput("CancerType", "Choose Cancer Type:",
                                         choices = c("Breast Cancer", "Pancreatic Cancer", "Glioma")),
                             submitButton("Update View"),
                             helpText("Note: Click this when network is not very good"),
                             numericInput("Datasets", "Choose an dataset:", value = 10, min = 1),
                             # conditionalPanel(condition = "input.CancerType == 'Pancreatic Cancer'",
                             #                  numericInput("Datasets", "Choose a dataset:",
                             #                               value = 3, min = 1, max = 5)),
                             # conditionalPanel(condition = "input.CancerType == 'Breast Cancer'",
                             #                  numericInput("Datasets", "Choose a dataset:",
                             #                               value = 10, min = 1, max = 12)),
                             # conditionalPanel(condition = "input.CancerType == 'Glioma'",
                             #                  numericInput("Datasets", "Choose a dataset:",
                             #                               value = 2, min = 1, max = 6)),
                             helpText("Note: The value is invalid when bigger than the view items"),
                             textInput("GeneName", "Choose a target gene:"),
                             sliderInput("obs", "Choose browser numbers", min = 0, max = 1000, value = 10),
                             submitButton("Update View")),
                         wellPanel(
                             helpText("Save raw clinicaldata"),
                             selectInput("Clinicaltabletype","Choose file type",
                                         choices = c("csv", "txt")),
                             downloadButton("Clinicaltable", "Download")
                             ),
                         wellPanel(
                             helpText("Save raw exprsetdata"),
                             selectInput("Exprsettabletype","Choose file type",
                                         choices = c("csv", "txt")),
                             downloadButton("Exprsettable", "Download")
                         )
                         ),
        
        conditionalPanel(condition = "input.MainMenu == 'Preprocess'",
                         wellPanel(
                             fileInput("file1", "Choose Clinicaldata",
                                       accept = c("text/csv", "text/comma-separated-values, text/plain", 
                                                  ".csv", "text/txt", ".txt")),
                             fileInput("file2", "Choose Exprsetdata",
                                       accept = c("text/csv", "text/comma-separated-values, text/plain", 
                                                  ".csv", "text/txt", ".txt")),
                             sliderInput("obser", "Choose browser numbers", min = 0, max = 100, value = 10),
                             # p(actionButton("Inputreadin","Read in")),
                             selectInput("TimeType", "Choose Survival Time Type:",
                                         choices = c("DFS_MONTHS", "DSS_MONTHS", "OS_MONTHS", "PFS_MONTHS")),
                             selectInput("StatusType", "Choose Survival Status Type:",
                                         choices = c("DFS_STATUS", "DSS_STATUS", "OS_STATUS", "PFS_STATUS")),
                             submitButton("Update View")
                             # tags$hr()
                         )),
        
        conditionalPanel(condition = "input.MainMenu == 'Difftest'",
                         wellPanel(
                             helpText("Save differential expression data"),
                             selectInput("Difftestplottype", "Choose plot type", choices = c("pdf", "tiff", "ps")),
                             textInput("Difftestfilewidth",  "Enter plot width",  12),
                             textInput("Difftestfileheight", "Enter plot height", 12),
                             ###### ? uioutput
                             uiOutput("Difftestfileheightui"),                                                                
                             downloadButton("Difftestsaveplot"),
                             tags$hr(),
                             submitButton("Update view")
                         )),
        
        conditionalPanel(condition = "input.MainMenu == 'kmplot'",
                         wellPanel(
                             helpText("Save KM plot data"),
                             selectInput("KMplottype", "Choose plot type", choices = c("pdf", "tiff", "ps")),
                             textInput("KMfilewidth",  "Enter plot width",  12),
                             textInput("KMfileheight", "Enter plot height", 12),
                             # downloadButton("kmplotsaveplot"),
                             tags$hr(),
                             submitButton("Update view")
                         ))
        
        
    ),
    
    # Show the K-M Plot of Pan-cancer
    mainPanel(
        
        conditionalPanel(condition = "input.MainMenu == 'kmplot'",
                         tabsetPanel(
                             tabPanel("KMTable", tableOutput("KMtable")),
                             tabPanel("KMPlot", plotOutput("KMplot"))
                         )
        ),

        conditionalPanel(condition = "input.MainMenu == 'Difftest'",
                         tabsetPanel(
                             tabPanel("DifftestPlot", plotOutput("Difftestplot"))
                         )
        ),
        
        conditionalPanel(condition = "input.MainMenu == 'Preprocess'",
                         tabsetPanel(
                             tabPanel("Datareadin",   tableOutput("readin")),
                             tabPanel("Datafastview", tableOutput("fastview")),
                             tabPanel("Clinicaldata", tableOutput("clinical")),
                             tabPanel("Exprsetdata",  tableOutput("exprset"))
                         )
                         # tableOutput("clinicaldat")
                         # uiOutput("Inputshowinstructionui"),
                         # uiOutput("Inputshowsummaryui")
                         ),
        
        conditionalPanel(condition = "input.MainMenu == 'Input'",
                         checkboxInput("Inputshowinstructiontf", "Show instructions", value = T),
                         uiOutput("Inputshowinstructionui"),
                         tabsetPanel(
                             tabPanel("Datasets", tableOutput("datasets")),
                             tabPanel("Clinicaldata", tableOutput("clinicaldataview")),
                             tabPanel("Summary", verbatimTextOutput("summary"), tableOutput("exprsetview"))
                         ))
        
        )
))


# setwd("I:/R-3.6.3/shinyapp/KMPPc")
# 
# # STEP 1 ¨C INSTALL RSCONNECT
# library("rsconnect")
# 
# # STEP 2 ¨C AUTHORIZE ACCOUNT
# rsconnect::setAccountInfo(name = 'zhouwei',
#                           token = 'A0106D487B1820ED76BDD00FE9ABC3B0',
#                           secret = 'KrIwYDBltOEzvbfwtMGu5NjfIuOAZQICGFPam2vw')
# 
# # STEP 3 ¨C DEPLOY
# shiny::runApp()
# rsconnect::deployApp(server = "shinyapps.io")





