######################################################
##                      KMPPc                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##               Author:Wei Zhou, xxx               ##
##       Maintainer:Wei Zhou (247328181@qq.com)     ##
######################################################

library("shiny")
library("cgdsr")
library("ggplot2")
library("stringr")
library("ggpubr")
library("survival")
library("survminer")
# library("ggstatsplot")
options(shiny.maxRequestSize = 30000*1024^2)


# Load data
########### ignore the link failure
mycgds <- CGDS("http://www.cbioportal.org/")
setVerbose(mycgds, TRUE)
mysecurecgds <- CGDS("http://cbioportal.mskcc.org/",
                     token = "fd0522cb-7972-40d0-9d83-cb4c14e8a337")
cancerstudy <- getCancerStudies(mycgds)

# mycancerstudy <- getCancerStudies(mycgds)[189, 1]
# mycaselist <- getCaseLists(mycgds, mycancerstudy)[6, 1]
# mygeneticprofile <- getGeneticProfiles(mycgds, mycancerstudy)[4, 1]

time <- as.character()
tmpstr <- str_extract_all(Sys.time(), "[0-9]")[[1]]
for (i in 1:length(tmpstr)) {
    time <- paste(time, tmpstr[i])
}
time <- gsub(" ", "", time)


# Define server logic 
shinyServer(function(input, output, session) {

    Mainstr <- reactiveValues()
    Maindat <- reactiveValues()
    
    observe({
        if(input$CancerType == "Breast Cancer"){
            cancerid <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "brca"))), 1]
            studydoi <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "brca"))), 2]
        } else if(input$CancerType == "Pancreatic Cancer"){
            cancerid <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "paad"))), 1]
            studydoi <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "paad"))), 2]
        } else if(input$CancerType == "Glioma"){
            cancerid <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "glioma|lgg"))), 1]
            studydoi <- cancerstudy[which(!is.na(str_extract(cancerstudy[ , 1], "glioma|lgg"))), 2]
        }
        Mainstr$id      <- cbind(cancerid, studydoi)
        Mainstr$dataset <- cancerid[input$Datasets]
    })
    observe({
        mycaselist       <- getCaseLists(mycgds, Mainstr$dataset)
        Mainstr$mrnanum  <- grep("mrna", mycaselist[ , 1])
        mygeneticprofile <- getGeneticProfiles(mycgds, Mainstr$dataset)
        Mainstr$mrnaNum  <- grep("mrna", mygeneticprofile[ , 1])[1] # (?!xx|xx)
    })
    
    ###########  If state to go or stop
    ### Input mainPanel
    # Input
    studyidInput <- reactive({
        as.character(Mainstr$dataset)
    })
    caselistInput <- reactive({
        as.character(getCaseLists(mycgds, studyidInput())[Mainstr$mrnanum, 1])
    })
    geneticprofileInput <- reactive({
        as.character(getGeneticProfiles(mycgds, studyidInput())[Mainstr$mrnaNum, 1])
    })
    exprsetInput <- reactive({
        if(is.null(input$GeneName)) return(NULL)
        getProfileData(mycgds, toupper(input$GeneName), geneticprofileInput(), caselistInput())
    })
    clinicalInput <- reactive({
        getClinicalData(mycgds, caselistInput())
    })
    
    # Output
    output$Inputshowinstructionui <- renderUI({
        if (input$Inputshowinstructiontf == T) {
            wellPanel(
                h5("Instructions:"),
                p("K-M Plot data should be prepared with clinical data and expression data. "),
                p("Step1-Reading in dataset: Choose a database to download clinical data and expression data (input your genename). "),
                p("Step2-Preprocessing: Upload the local files that you have downloaded, respectively. "),
                p("Step3-K-M Plot: Refresh the page and automatically show the plot."),
                h5("KMPPc BETA version:"),
                p("A K-M Plot tool for database discovery..."),
                p("This webpage is all constructed by R, and any code will be shared if you are interested."),
                p("Note: The numbers of patients in selected database must be as many as possible. To avoid 404 and operate locally, please download the clinical data and expression data.")                      
            )
        }
    })
    output$datasets <- renderTable(
        data.matrix(Mainstr$id), 
        striped = T, hover = T, borderd = T, spacing = "s", align = "l", rownames = T, colnames = F
    )
    output$summary <- renderPrint({
        summary(exprsetInput())
    })
    output$exprsetview <- renderTable({
        head(exprsetInput(), n = input$obs)
    })
    output$clinicaldataview <- renderTable(
        head(clinicalInput(), n = input$obs),
        striped = T, hover = T, bordered = T, spacing = "xs", align = "l"
    )
    output$Clinicaltable <- downloadHandler(
        filename = function() {paste0("Clinicaldata", time, ".", input$Clinicaltabletype)},
        content = function(file) {
            if(input$Clinicaltabletype == "txt") {
                tmpdat <- write.table(clinicalInput(), file, row.names = T, quote = F)
                tmpdat[!duplicated(rownames(tmpdat)), ]
            } else if(input$Clinicaltabletype == "csv") {
                tmpdat <- write.csv(clinicalInput(), file, row.names = T, quote = F)
                tmpdat[!duplicated(rownames(tmpdat)), ]
            }
        }
    )
    output$Exprsettable <- downloadHandler(
        filename = function() {paste0("Exprsetdata", time, ".", input$Exprsettabletype)},
        content = function(file) {
            if(input$Exprsettabletype == "txt") {
                tmpdat <- write.table(exprsetInput(), file, row.names = T, quote = F)
                tmpdat[!duplicated(rownames(tmpdat)), ]
            } else if(input$Exprsettabletype == "csv") {
                tmpdat <- write.csv(exprsetInput(), file, row.names = T, quote = F)
                tmpdat[!duplicated(rownames(tmpdat)), ]
            }
        }
    )
    
    
    ### Preprocess mainPanel
    # Input
    observe({
        inFile1 <- input$file1
        inFile2 <- input$file2
        ifelse(is.null(inFile1), return(NULL), myclinicaldat <- read.csv(inFile1$datapath))
        ifelse(is.null(inFile2), return(NULL), myexprsetdat  <- read.csv(inFile2$datapath))
        myclinicaldata <- as.matrix(myclinicaldat)
        myexprsetdata  <- as.matrix(myexprsetdat)
        rownames(myexprsetdata) <- myexprsetdata[ , 1]
        clinicaldata_view <- as.matrix(colnames(myclinicaldata))
        choose_columns <- c(clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DFS_MONTHS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DFS_STATUS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DSS_MONTHS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "DSS_STATUS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "OS_MONTHS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "OS_STATUS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "PFS_MONTHS")))],
                            clinicaldata_view[which(!is.na(str_extract(clinicaldata_view, "PFS_STATUS")))])
        choose_clinicaldata <- myclinicaldata[ , choose_columns]
        rownames(choose_clinicaldata) <- myclinicaldata[ , 1]
        Maindat$label <- choose_clinicaldata
        timetype <- input$TimeType
        if(timetype == "DFS_MONTHS") {
            timestatus = "DFS_STATUS"
        } else if(timetype == "DSS_MONTHS") {
            timestatus = "DSS_STATUS"
        } else if(timetype == "OS_MONTHS") {
            timestatus = "OS_STATUS"
        } else if(timetype == "PFS_MONTHS") {
            timestatus = "PFS_STATUS"
        }
        Maindat$status <- timestatus
        dat1 <- choose_clinicaldata[!is.na(choose_clinicaldata[ , timetype]), ]
        ########### timetype bug when changing to other types
        dat2 <- myexprsetdata[rownames(dat1), ]  
        colnames(dat2)[1] <- "PatientID"
        Maindat$rawdata <- dat2
        Maindat$genename <- colnames(dat2)[2]
        Maindat$clinicaldat <- cbind(dat2, dat1[ , c(timetype, timestatus)])
    })

    # Output
    output$readin <- renderTable({
        if(is.null(input$file1)) return(NULL)
        data.frame(Maindat$label)[1:input$obser, ]},
        striped = T, hover = T, bordered = T, spacing = "xs", align = "l", digits = 2)
    output$fastview <- renderTable({
        if(is.null(input$file1)) return(NULL)
        if(is.null(input$file2)) return(NULL)
        data.frame(Maindat$clinicaldat)[1:input$obser, ]},
        striped = T, hover = T, bordered = T, spacing = "xs", align = "l", digits = 2)
    output$clinical <- renderTable({
        inFile <- input$file1
        if(is.null(inFile)) return(NULL)
        
        read.csv(inFile$datapath)},
        striped = T, hover = T, bordered = T, spacing = "xs", align = "l", digits = 2)
    output$exprset <- renderTable({
        inFile <- input$file2
        if(is.null(inFile)) return(NULL)
        
        read.csv(inFile$datapath)},
        striped = T, hover = T, bordered = T, spacing = "xs", align = "l", digits = 2)
    
    
    ### Difftest mainPanel
    # Input
    statusInput  <- reactive({
        as.character(Maindat$status)
    })
    geneidInput  <- reactive({
        as.character(Maindat$genename)
    })
    
    # output
    output$Difftestplot <- renderPlot({
        dat <- data.frame(Maindat$rawdata)
        if(!is.null(dat)) {
            gene <- geneidInput()
            gene_group <- paste0(gene, "_group")
            genes = as.numeric(dat[ , 2])
            gene_groups = ifelse(genes > median(genes), 'high', 'low')
            # Maindat$plot <- ggboxplot(dat2, 
            #                           x = "DFS_STATUS", y = "JAG1", 
            #                           color = "DFS_STATUS", palette = "jco", add = "jitter")
            dat$genes <- genes; dat$gene_groups <- gene_groups
            Maindat$plot <- ggbetweenstats(data = dat,
                                           x = gene_groups,
                                           y = genes,

                                           xlab = gene,
                                           ylab = gene_group)
            ggbetweenstats(data = dat,
                           x = gene_groups,
                           y = genes,

                           xlab = gene,
                           ylab = gene_group)
        }
    })
    output$Difftestsaveplot <- downloadHandler(
        filename = function() { paste('Diff_expression.', input$Difftestplottype) },
        content = function(file) {
            if (input$Difftestplottype == "pdf") {
                pdf(file, width = as.numeric(input$Difftestfilewidth), height = as.numeric(input$Difftestfileheight))
            } else if (input$Difftestplottype == "tiff") {
                tiff(file, width = 700, height = 800)
            } else if (input$Difftestplottype == "ps") {
                postscript(file, width = as.numeric(input$Difftestfilewidth), height = as.numeric(input$Difftestfileheight), paper = "special")
            }                  
            plot(Maindat$plot)
            dev.off()
        }
    )
    
    
    ### K-M Plot mainPanel
    # Input
    observe({
        data <- Maindat$clinicaldat
        dat <- data.frame(data)
        if(!is.null(dat)) {
            # time <- input$TimeType; status <- Maindat$status
            gene <- geneidInput()
            gene_group <- paste0(gene, "_group")
            genes = as.numeric(dat[ , gene])
            gene_groups = ifelse(genes > median(genes), 'high', 'low')
            highnum <- 
            geneAndgroup <- cbind(genes, gene_groups)
            dat <- cbind(dat, geneAndgroup)
            dat1 <- data.frame(dat)
            dat$DFS_MONTHS <- as.numeric(dat$DFS_MONTHS)
            my.surv <<- Surv(dat$DFS_MONTHS, dat$DFS_STATUS == 'Recurred/Progressed')
            kmfit <<- survfit(my.surv~gene_groups, data = dat)
        }
        Maindat$kmfit <- kmfit
        Maindat$kmdat <- dat
    })
    kmfitInput <- reactive({
        Maindat$kmfit
    })
    kmdatInput <- reactive({
        Maindat$kmdat
    })
    
    # output
    output$KMtable <- renderTable({
        surv_summary(kmfitInput(), data = Maindat$kmdat)
    })
    output$KMplot  <- renderPlot({
        kmfit <- kmfitInput()
        names(kmfit$strata) <- c(paste0(geneidInput(), " high"), paste0(geneidInput(), " low"))
        label <- data.frame(kmfit$strata)
        labs <- cbind(rownames(label), label)
        colnames(labs) <- c("gene_group", "n")
        gene_high <- paste0(labs[1, 1] , "\n (n = ", labs[1, 2], ")")
        gene_low  <- paste0(labs[2, 1] , "\n (n = ", labs[2, 2], ")")
        names(kmfit$strata) <- c(gene_high, gene_low)
        g <- theme_survminer(font.main = c(24, "bold", "darkblue"),
                             font.x = c(18, "bold", "black"),
                             font.y = c(18, "bold", "black"),
                             font.tickslab = c(18, "bold", "black")) +
            theme(axis.line.x = element_line(size = 1.2),
                  axis.line.y = element_line(size = 1.2),
                  axis.ticks.x = element_line(size = 1.2),
                  axis.ticks.y = element_line(size = 1.2),
                  plot.title = element_text(hjust = 0.5),
                  legend.title = element_blank(),
                  legend.direction = "horizontal",
                  legend.text = element_text(size = 16, face = "bold"))
        Maindat$kmplot <- ggsurvplot(kmfit, data = Maindat$kmdat, palette = c("red", "blue"),
                                     conf.int = F, pval = T, pval.method = T, surv.scale = "percent",
                                     title = paste0("K-M plot of ", geneidInput(), " : ", studyidInput()),
                                     xlab ="Months" , ylab = "Disease-free", censor = F,
                                     ggtheme = g, risk.table = F, ncensor.plot = F)
        ggsurvplot(kmfit, data = Maindat$kmdat, palette = c("red", "blue"),
                   conf.int = F, pval = T, pval.method = T, surv.scale = "percent",
                   title = paste0("K-M plot of ", geneidInput(), " : ", studyidInput()),
                   xlab ="Months" , ylab = "Disease-free", censor = F,
                   ggtheme = g, risk.table = F, ncensor.plot = F)
    })
    output$kmplotsaveplot <- downloadHandler(
        filename = function() { paste('KM_Plot.', input$KMplottype) },
        content = function(file) {
            if (input$KMplottype == "pdf") {
                pdf(file, width = as.numeric(input$KMfilewidth), height = as.numeric(input$KMfileheight))
            } else if (input$KMplottype == "tiff") {
                tiff(file, width = 700, height = 800)
            } else if (input$KMplottype == "ps") {
                postscript(file, width = as.numeric(input$KMfilewidth), height = as.numeric(input$KMfileheight), paper = "special")
            }
            plot(Maindat$kmplot)
            dev.off()
        }
    )
    
})
