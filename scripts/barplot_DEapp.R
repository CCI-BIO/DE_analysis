## barplot_DEapp.R
## Script containing functions to allow the construction of enrichment bar plots for the DE analysis application.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

# List files
barPlotOutFileList <- eventReactive(c(input$OutputDir, input$OutputDirNew, input$tabs, input$fileChoice2, outputLocation$output),
                                    {
                                      if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                                        if(!is.null(input$fileChoice2)){
                                          compareChoice <- input$fileChoice2
                                          compareChoice <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", compareChoice)
                                          compareChoice <- gsub(pattern = "_DEanalysis", replacement = "", compareChoice)
                                          outDir <- outputLocation$output
                                          fileList <- list.files(path = outDir, pattern = ".txt", full.names = T)
                                          fileList <- fileList[grep(pattern = compareChoice, x = fileList)]
                                          fileList <- fileList[grep(pattern = "enrichment_DEanalysis.txt", x = fileList)]
                                          if(length(fileList) == 0){
                                            return(NULL)
                                          } else {
                                            tryCatch({barplot.data <- read.delim(fileList[1], header = T, stringsAsFactors = F, sep = "\t")
                                            return(fileList)},
                                            error=function(cond){
                                              UPDOWNcheck <<- FALSE
                                              return(NULL)})
                                          }
                                        }
                                      }
                                    }
)

observeEvent(c(input$OutputDir, input$OutputDirNew, input$tabs, input$fileChoice2, outputLocation$output), {
  if(!is.null(barPlotOutFileList())){
    fileList <- barPlotOutFileList()
    checkForDirection <- length(grep("DOWN", fileList))
    if(checkForDirection == 0){
      UPDOWNcheck <<- FALSE
    } else {
      UPDOWNcheck <<- TRUE
    }
  }
  if(UPDOWNcheck == TRUE){
    output$chooseDOWNUP <- renderUI({
      radioButtons("chooseDOWNUP2", h4("Plot UP or DOWN regulated results:"),
                   choices = list("UP-regulated" = "UP",
                                  "DOWN-regulated" = "DOWN"))
    })
  } else {
    output$chooseDOWNUP <- renderUI({})
  }
})

observe({
  if(length(input$combineGOKEGG) == 1){
    barPlotGroupLastSelected <<- input$combineGOKEGG
  }
  if(length(input$combineGOKEGG) < 1){
    updateCheckboxGroupInput(session, "combineGOKEGG", selected = barPlotGroupLastSelected)
  }
})

shinyjs::disable("refresh4")

observe({
  fileList <- barPlotOutFileList()
  if(!is.null(fileList)){
    tryCatch({read.delim(fileList[1], header = T, sep = "\t", stringsAsFactors = F)
      shinyjs::enable("refresh4")},
      error = function(cond){shinyjs::disable("refresh4")})
  } else {
    shinyjs::disable("refresh4")
  }
})

shinyjs::disable("downloadbarPlot")

observeEvent(c(input$barPlotFilename, input$refresh4),{
  if(input$barPlotFilename == ""){
    disable("downloadbarPlot")
  }
  if(barPlotrefresh == TRUE){
    if(input$barPlotFilename != ""){
      enable("downloadbarPlot")
    } else {
      disable("downloadbarPlot")
    }
  }
})

observe({
  fileList <- barPlotOutFileList()
  if(!is.null(fileList)){
    if(length(grep("enrichR", fileList)) > 0){
      output$rownameChoice <-renderUI({})
    } else {
      output$rownameChoice <- renderUI({radioButtons("rownameChoice2", h4("Row name style:"),
                                                     choices = list("Term ID + Description" = "termdesc",
                                                                    "Term ID only" = "term",
                                                                    "Description only" = "desc"))})
    }
  } else {
    output$rownameChoice <-renderUI({})
  }
  if(input$plotMethod == "filter"){
    output$termChoice <- renderUI({
      radioButtons("termChoice2", h4("Statistic to filter terms:"),
                   choices = list("P-value" = "pvalue",
                                  "Adjusted P-value" = "p.adjust"))
    })
    output$termThreshold <- renderUI({
      tryCatch({
        if(input$termChoice2 == "pvalue"){
          termTitle <- "P-value threshold:"
        } else {
          termTitle <- "Adj. P-value threshold:"
        }
        numericInput("termThreshold2", h4(termTitle), value = 0.05, min = 0, max = 1)
      },
      error = function(cond){return(NULL)}
      )
    })
    output$barPlotTermListOption <- renderUI({})
    output$barPlotTermList <- renderUI({})
  } else {
    output$barPlotTermListOption <- renderUI({
      radioButtons("barPlotTermListOption2", h4("Method to select terms:"),
                   choices = list("Top significantly enriched terms" = "topTerms",
                                  "Select specified terms" = "selectInput"))
    })
    output$barPlotTermList <- renderUI({
      if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
        tryCatch({
          if(input$barPlotTermListOption2 == "topTerms"){
            numericInput("topTermsNumber2",
                         h4("Number of top terms to plot:"),
                         min = 1, max = max.row, step = 1,
                         value = 10)
          } else if(input$barPlotTermListOption2 == "selectInput"){
            selectizeInput("selectInputTerms2", h4("Select terms to plot:"), choices = NULL, multiple = TRUE)
          }
        },
        error = function(cond){return(NULL)}
        )
      }
    })
    output$termChoice <- renderUI({
      tryCatch({
        if(input$barPlotTermListOption2 == "topTerms" | input$barPlotTermListOption2 == "selectInput"){
          if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
            radioButtons("termChoice2", h4("Statistic to filter terms:"),
                         choices = list("P-value" = "pvalue",
                                        "Adjusted P-value" = "p.adjust"))
          }
        }
      },
      error = function(cond){return(NULL)}
      )
    })
    output$termThreshold <- renderUI({})
  }
})

observeEvent(c(input$fileChoice2, input$barPlotTermListOption2, input$chooseDOWNUP2, input$combineGOKEGG),
             {
               fileList <- barPlotOutFileList()
               if(!is.null(fileList)){
                 if(length(grep("enrichR", fileList)) > 0){
                   column.names <- c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value","Odds.Ratio","Combined.Score","Genes","termType")
                 } else {
                   column.names <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","termType")
                 }
                 if(UPDOWNcheck == TRUE){
                   tryCatch({
                     directionChoice <- as.character(input$chooseDOWNUP2)
                     fileList <- fileList[grep(pattern = directionChoice, x = fileList)]
                     KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
                     GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
                     GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
                     GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
                     termsToInclude <- input$combineGOKEGG
                     barPlotTable <- matrix(ncol = 10)
                     colnames(barPlotTable) <- column.names
                     if("kegg" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, KEGG_table)
                     }
                     if("gobp" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, GO_BP_table)
                     }
                     if("gocc" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, GO_CC_table)
                     }
                     if("gomf" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, GO_MF_table)
                     }
                     barPlotTable <- barPlotTable[-1,]
                     
                     if(length(grep("enrichR", fileList)) > 0){
                       term_list <- barPlotTable$Term
                       names(term_list) <- barPlotTable$Term
                       term_list <- as.list(term_list)
                       updateSelectizeInput(session, "selectInputTerms2", choices = term_list, server = TRUE)
                     } else {
                       term_list <- barPlotTable$ID
                       names(term_list) <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
                       term_list <- as.list(term_list)
                       updateSelectizeInput(session, "selectInputTerms2", choices = term_list, server = TRUE)
                     }
                   },
                   error = function(cond){return(NULL)}
                   )
                 } else {
                   tryCatch({
                     KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
                     GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
                     GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
                     GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
                     GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
                     termsToInclude <- input$combineGOKEGG
                     barPlotTable <- matrix(ncol = 10)
                     colnames(barPlotTable) <- column.names
                     if("kegg" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, KEGG_table)
                     }
                     if("gobp" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, GO_BP_table)
                     }
                     if("gocc" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, GO_CC_table)
                     }
                     if("gomf" %in% termsToInclude){
                       barPlotTable <- rbind(barPlotTable, GO_MF_table)
                     }
                     barPlotTable <- barPlotTable[-1,]
                     if(length(grep("enrichR", fileList)) > 0){
                       term_list <- barPlotTable$Term
                       names(term_list) <- barPlotTable$Term
                       term_list <- as.list(term_list)
                       updateSelectizeInput(session, "selectInputTerms2", choices = term_list, server = TRUE)
                     } else {
                       term_list <- barPlotTable$ID
                       names(term_list) <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
                       term_list <- as.list(term_list)
                       updateSelectizeInput(session, "selectInputTerms2", choices = term_list, server = TRUE)
                     }
                   },
                   error = function(cond){return(NULL)}
                   )
                 }
               }
             })


# Generate Bar Plot
observeEvent({input$refresh4}, {
  barPlotrefresh <<- TRUE
  fileList <- barPlotOutFileList()
  plotMethod <- input$plotMethod
  if(length(grep("enrichR", fileList)) > 0){
    column.names <- c("Term","Overlap","P.value","Adjusted.P.value","Old.P.value","Old.Adjusted.P.value","Odds.Ratio","Combined.Score","Genes","termType")
  } else {
    column.names <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","termType")
  }
  addOverlapCount <- input$addOverlapCount
  reverseOrder <- input$reverseOrder
  if(plotMethod == "filter"){
    termChoice <- input$termChoice2
    if(length(grep("enrichR", fileList)) > 0){
      if(termChoice == "pvalue"){
        termChoice <- "P.value"
      }
      if(termChoice == "p.adjust"){
        termChoice <- "Adjusted.P.value"
      }
    }
    termThreshold <- input$termThreshold2
    if(!is.null(fileList)){
      if(UPDOWNcheck == TRUE){
        directionChoice <- as.character(input$chooseDOWNUP2)
        fileList <- fileList[grep(pattern = directionChoice, x = fileList)]
        KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
        GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
        GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
        GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
        termsToInclude <- input$combineGOKEGG
        barPlotTable <- matrix(ncol = 10)
        colnames(barPlotTable) <- column.names
        if("kegg" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, KEGG_table)
        }
        if("gobp" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, GO_BP_table)
        }
        if("gocc" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, GO_CC_table)
        }
        if("gomf" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, GO_MF_table)
        }
        barPlotTable <- barPlotTable[-1,]
        barPlotTable <- barPlotTable[order(barPlotTable[,termChoice], decreasing = F),]
        barPlotTable <- barPlotTable[which(barPlotTable[,termChoice] <= termThreshold),]
        
        if(length(grep("enrichR", fileList)) > 0){
          term_vector <- barPlotTable$Term
        } else {
          geneNum <- strsplit(barPlotTable$GeneRatio, split = "/")
          bgNum <- strsplit(barPlotTable$BgRatio, split = "/")
          newOverlap <- c()
          for(i in 1:length(geneNum)){
            newOverlap[i] <- paste(geneNum[[i]][1], "/", bgNum[[i]][1], sep = "")
          }
          barPlotTable$Overlap <- newOverlap
          if(input$rownameChoice2 == "termdesc"){
            term_vector <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
          }
          if(input$rownameChoice2 == "term"){
            term_vector <- barPlotTable$ID
          }
          if(input$rownameChoice2 == "desc"){
            term_vector <- barPlotTable$Description
            if(any(duplicated(term_vector))){
              idx <- which(duplicated(term_vector))
              terms <- term_vector[idx]
              for(i in 1:length(terms)){
                idx2 <- which(term_vector == terms[i])
                termTypes <- barPlotTable$termType[idx2]
                for(j in 1:length(termTypes)){
                  if(termTypes[j] == "KEGG Pathways"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(KEGG)")
                  }
                  if(termTypes[j] == "Biological Processes"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOBP)")
                  }
                  if(termTypes[j] == "Cellular Components"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOCC)")
                  }
                  if(termTypes[j] == "Molecular Function"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOMF)")
                  }
                }
              }
            }
          }
        }
        
        barPlotData <- data.frame("Category" = barPlotTable$termType, "Term" = term_vector, "Statistic" = log(barPlotTable[,termChoice]), "Overlap" = barPlotTable$Overlap)
        if(termChoice == "pvalue" | termChoice == "P.value"){
          ylabel <- "log(P-value)"
        }
        if(termChoice == "p.adjust" | termChoice == "Adjusted.P.value"){
          ylabel <- "log(Adj. P-value)"
        }
        if(reverseOrder){
          p <- ggplot(data = barPlotData, aes(x = reorder(Term, Statistic), y = Statistic, fill = Category)) +
            geom_bar(stat = "identity", width = 0.8)+
            theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
            ylab(ylabel)+xlab("Term")
          p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
          if(addOverlapCount){
            p <- p + geom_text(data = barPlotData, 
                               aes(x = reorder(Term, Statistic), y = Statistic, label = Overlap), 
                               hjust = +1.1, 
                               size = 4,
                               inherit.aes = TRUE)
          }
        } else {
          p <- ggplot(data = barPlotData, aes(x = reorder(Term, -Statistic), y = Statistic, fill = Category)) +
            geom_bar(stat = "identity", width = 0.8)+
            theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
            ylab(ylabel)+xlab("Term")
          p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
          if(addOverlapCount){
            p <- p + geom_text(data = barPlotData, 
                               aes(x = reorder(Term, -Statistic), y = Statistic, label = Overlap), 
                               hjust = +1.1, 
                               size = 4,
                               inherit.aes = TRUE)
          }
        }
        output$barPlot <- renderPlot({p}, height = input$barPlotHeight, width = input$barPlotWidth)
        barPlotHeight <- input$barPlotHeight*(2.54/100)
        barPlotWidth <- input$barPlotWidth*(2.54/100)
        output$downloadbarPlot <- downloadHandler(
          filename = function(){paste(input$barPlotFilename, '.png', sep = "")},
          content = function(file){ggsave(file, plot = p, device = "png", height = barPlotHeight, width = barPlotWidth, dpi = 100, units = "cm")}
        )
      } else {
        KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
        GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
        GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
        GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
        GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
        termsToInclude <- input$combineGOKEGG
        barPlotTable <- matrix(ncol = 10)
        colnames(barPlotTable) <- column.names
        if("kegg" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, KEGG_table)
        }
        if("gobp" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, GO_BP_table)
        }
        if("gocc" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, GO_CC_table)
        }
        if("gomf" %in% termsToInclude){
          barPlotTable <- rbind(barPlotTable, GO_MF_table)
        }
        barPlotTable <- barPlotTable[-1,]
        barPlotTable <- barPlotTable[order(barPlotTable[,termChoice], decreasing = F),]
        barPlotTable <- barPlotTable[which(barPlotTable[,termChoice] <= termThreshold),]
        
        if(length(grep("enrichR", fileList)) > 0){
          term_vector <- barPlotTable$Term
        } else {
          geneNum <- strsplit(barPlotTable$GeneRatio, split = "/")
          bgNum <- strsplit(barPlotTable$BgRatio, split = "/")
          newOverlap <- c()
          for(i in 1:length(geneNum)){
            newOverlap[i] <- paste(geneNum[[i]][1], "/", bgNum[[i]][1], sep = "")
          }
          barPlotTable$Overlap <- newOverlap
          if(input$rownameChoice2 == "termdesc"){
            term_vector <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
          }
          if(input$rownameChoice2 == "term"){
            term_vector <- barPlotTable$ID
          }
          if(input$rownameChoice2 == "desc"){
            term_vector <- barPlotTable$Description
            if(any(duplicated(term_vector))){
              idx <- which(duplicated(term_vector))
              terms <- term_vector[idx]
              for(i in 1:length(terms)){
                idx2 <- which(term_vector == terms[i])
                termTypes <- barPlotTable$termType[idx2]
                for(j in 1:length(termTypes)){
                  if(termTypes[j] == "KEGG Pathways"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(KEGG)")
                  }
                  if(termTypes[j] == "Biological Processes"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOBP)")
                  }
                  if(termTypes[j] == "Cellular Components"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOCC)")
                  }
                  if(termTypes[j] == "Molecular Function"){
                    term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOMF)")
                  }
                }
              }
            }
          }
        }
        
        barPlotData <- data.frame("Category" = barPlotTable$termType, "Term" = term_vector, "Statistic" = log(barPlotTable[,termChoice]), "Overlap" = barPlotTable$Overlap)
        if(termChoice == "pvalue" | termChoice == "P.value"){
          ylabel <- "log(P-value)"
        }
        if(termChoice == "p.adjust" | termChoice == "Adjusted.P.value"){
          ylabel <- "log(Adj. P-value)"
        }
        if(reverseOrder){
          p <- ggplot(data = barPlotData, aes(x = reorder(Term, Statistic), y = Statistic, fill = Category)) +
            geom_bar(stat = "identity", width = 0.8)+
            theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
            ylab(ylabel)+xlab("Term")
          p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
          if(addOverlapCount){
            p <- p + geom_text(data = barPlotData, 
                               aes(x = reorder(Term, Statistic), y = Statistic, label = Overlap), 
                               hjust = +1.1, 
                               size = 4,
                               inherit.aes = TRUE)
          }
        } else {
          p <- ggplot(data = barPlotData, aes(x = reorder(Term, -Statistic), y = Statistic, fill = Category)) +
            geom_bar(stat = "identity", width = 0.8)+
            theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
            ylab(ylabel)+xlab("Term")
          p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
          if(addOverlapCount){
            p <- p + geom_text(data = barPlotData, 
                               aes(x = reorder(Term, -Statistic), y = Statistic, label = Overlap), 
                               hjust = +1.1, 
                               size = 4,
                               inherit.aes = TRUE)
          }
        }
        output$barPlot <- renderPlot({p}, height = input$barPlotHeight, width = input$barPlotWidth)
        barPlotHeight <- input$barPlotHeight*(2.54/100)
        barPlotWidth <- input$barPlotWidth*(2.54/100)
        output$downloadbarPlot <- downloadHandler(
          filename = function(){paste(input$barPlotFilename, '.png', sep = "")},
          content = function(file){ggsave(file, plot = p, device = "png", height = barPlotHeight, width = barPlotWidth, dpi = 100, units = "cm")}
        )
      }
    }
  } else if(plotMethod == "select"){
    if(input$barPlotTermListOption2 == "topTerms"){
      termChoice <- input$termChoice2
      if(length(grep("enrichR", fileList)) > 0){
        if(termChoice == "pvalue"){
          termChoice <- "P.value"
        }
        if(termChoice == "p.adjust"){
          termChoice <- "Adjusted.P.value"
        }
      }
      if(!is.null(fileList)){
        if(UPDOWNcheck == TRUE){
          directionChoice <- as.character(input$chooseDOWNUP2)
          fileList <- fileList[grep(pattern = directionChoice, x = fileList)]
          KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
          GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
          GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
          GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
          termsToInclude <- input$combineGOKEGG
          barPlotTable <- matrix(ncol = 10)
          colnames(barPlotTable) <- column.names
          if("kegg" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, KEGG_table)
          }
          if("gobp" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_BP_table)
          }
          if("gocc" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_CC_table)
          }
          if("gomf" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_MF_table)
          }
          barPlotTable <- barPlotTable[-1,]
          barPlotTable <- barPlotTable[order(barPlotTable[,termChoice], decreasing = F),]
          numberOfTerms <- input$topTermsNumber2
          if(numberOfTerms > nrow(barPlotTable)){
            numberOfTerms <- nrow(barPlotTable)
          }
          barPlotTable <- barPlotTable[1:numberOfTerms,,drop=FALSE]
          
          if(length(grep("enrichR", fileList)) > 0){
            term_vector <- barPlotTable$Term
          } else {
            geneNum <- strsplit(barPlotTable$GeneRatio, split = "/")
            bgNum <- strsplit(barPlotTable$BgRatio, split = "/")
            newOverlap <- c()
            for(i in 1:length(geneNum)){
              newOverlap[i] <- paste(geneNum[[i]][1], "/", bgNum[[i]][1], sep = "")
            }
            barPlotTable$Overlap <- newOverlap
            if(input$rownameChoice2 == "termdesc"){
              term_vector <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
            }
            if(input$rownameChoice2 == "term"){
              term_vector <- barPlotTable$ID
            }
            if(input$rownameChoice2 == "desc"){
              term_vector <- barPlotTable$Description
              if(any(duplicated(term_vector))){
                idx <- which(duplicated(term_vector))
                terms <- term_vector[idx]
                for(i in 1:length(terms)){
                  idx2 <- which(term_vector == terms[i])
                  termTypes <- barPlotTable$termType[idx2]
                  for(j in 1:length(termTypes)){
                    if(termTypes[j] == "KEGG Pathways"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(KEGG)")
                    }
                    if(termTypes[j] == "Biological Processes"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOBP)")
                    }
                    if(termTypes[j] == "Cellular Components"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOCC)")
                    }
                    if(termTypes[j] == "Molecular Function"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOMF)")
                    }
                  }
                }
              }
            }
          }
          
          barPlotData <- data.frame("Category" = barPlotTable$termType, "Term" = term_vector, "Statistic" = log(barPlotTable[,termChoice]), "Overlap" = barPlotTable$Overlap)
          if(termChoice == "pvalue" | termChoice == "P.value"){
            ylabel <- "log(P-value)"
          }
          if(termChoice == "p.adjust" | termChoice == "Adjusted.P.value"){
            ylabel <- "log(Adj. P-value)"
          }
          if(reverseOrder){
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          } else {
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, -Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, -Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          }
          output$barPlot <- renderPlot({p}, height = input$barPlotHeight, width = input$barPlotWidth)
          barPlotHeight <- input$barPlotHeight*(2.54/100)
          barPlotWidth <- input$barPlotWidth*(2.54/100)
          output$downloadbarPlot <- downloadHandler(
            filename = function(){paste(input$barPlotFilename, '.png', sep = "")},
            content = function(file){ggsave(file, plot = p, device = "png", height = barPlotHeight, width = barPlotWidth, dpi = 100, units = "cm")}
          )
        } else {
          KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
          GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
          GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
          GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
          termsToInclude <- input$combineGOKEGG
          barPlotTable <- matrix(ncol = 10)
          colnames(barPlotTable) <- column.names
          if("kegg" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, KEGG_table)
          }
          if("gobp" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_BP_table)
          }
          if("gocc" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_CC_table)
          }
          if("gomf" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_MF_table)
          }
          barPlotTable <- barPlotTable[-1,]
          barPlotTable <- barPlotTable[order(barPlotTable[,termChoice], decreasing = F),]
          numberOfTerms <- input$topTermsNumber2
          if(numberOfTerms > nrow(barPlotTable)){
            numberOfTerms <- nrow(barPlotTable)
          }
          barPlotTable <- barPlotTable[1:numberOfTerms,,drop=FALSE]
          
          if(length(grep("enrichR", fileList)) > 0){
            term_vector <- barPlotTable$Term
          } else {
            geneNum <- strsplit(barPlotTable$GeneRatio, split = "/")
            bgNum <- strsplit(barPlotTable$BgRatio, split = "/")
            newOverlap <- c()
            for(i in 1:length(geneNum)){
              newOverlap[i] <- paste(geneNum[[i]][1], "/", bgNum[[i]][1], sep = "")
            }
            barPlotTable$Overlap <- newOverlap
            if(input$rownameChoice2 == "termdesc"){
              term_vector <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
            }
            if(input$rownameChoice2 == "term"){
              term_vector <- barPlotTable$ID
            }
            if(input$rownameChoice2 == "desc"){
              term_vector <- barPlotTable$Description
              if(any(duplicated(term_vector))){
                idx <- which(duplicated(term_vector))
                terms <- term_vector[idx]
                for(i in 1:length(terms)){
                  idx2 <- which(term_vector == terms[i])
                  termTypes <- barPlotTable$termType[idx2]
                  for(j in 1:length(termTypes)){
                    if(termTypes[j] == "KEGG Pathways"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(KEGG)")
                    }
                    if(termTypes[j] == "Biological Processes"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOBP)")
                    }
                    if(termTypes[j] == "Cellular Components"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOCC)")
                    }
                    if(termTypes[j] == "Molecular Function"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOMF)")
                    }
                  }
                }
              }
            }
          }
          
          barPlotData <- data.frame("Category" = barPlotTable$termType, "Term" = term_vector, "Statistic" = log(barPlotTable[,termChoice]), "Overlap" = barPlotTable$Overlap)
          if(termChoice == "pvalue" | termChoice == "P.value"){
            ylabel <- "log(P-value)"
          }
          if(termChoice == "p.adjust" | termChoice == "Adjusted.P.value"){
            ylabel <- "log(Adj. P-value)"
          }
          if(reverseOrder){
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          } else {
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, -Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, -Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          }
          output$barPlot <- renderPlot({p}, height = input$barPlotHeight, width = input$barPlotWidth)
          barPlotHeight <- input$barPlotHeight*(2.54/100)
          barPlotWidth <- input$barPlotWidth*(2.54/100)
          output$downloadbarPlot <- downloadHandler(
            filename = function(){paste(input$barPlotFilename, '.png', sep = "")},
            content = function(file){ggsave(file, plot = p, device = "png", height = barPlotHeight, width = barPlotWidth, dpi = 100, units = "cm")}
          )
        }
      }
    } else if(input$barPlotTermListOption2 == "selectInput"){
      termChoice <- input$termChoice2
      if(length(grep("enrichR", fileList)) > 0){
        if(termChoice == "pvalue"){
          termChoice <- "P.value"
        }
        if(termChoice == "p.adjust"){
          termChoice <- "Adjusted.P.value"
        }
      }
      if(!is.null(fileList)){
        if(UPDOWNcheck == TRUE){
          directionChoice <- as.character(input$chooseDOWNUP2)
          fileList <- fileList[grep(pattern = directionChoice, x = fileList)]
          KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
          GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
          GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
          GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
          termsToInclude <- input$combineGOKEGG
          barPlotTable <- matrix(ncol = 10)
          colnames(barPlotTable) <- column.names
          if("kegg" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, KEGG_table)
          }
          if("gobp" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_BP_table)
          }
          if("gocc" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_CC_table)
          }
          if("gomf" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_MF_table)
          }
          barPlotTable <- barPlotTable[-1,]
          barPlotTable <- barPlotTable[order(barPlotTable[,termChoice], decreasing = F),]
          
          selectedTerms <- input$selectInputTerms2
          
          if(length(grep("enrichR", fileList)) > 0){
            if(!is.null(selectedTerms)){
              barPlotTable <- barPlotTable[which(barPlotTable$Term %in% selectedTerms),]
            } else {
              barPlotTable <- barPlotTable[-c(1:nrow(barPlotTable)),,drop=FALSE]
            }
            term_vector <- barPlotTable$Term
          } else {
            if(!is.null(selectedTerms)){
              barPlotTable <- barPlotTable[which(barPlotTable$ID %in% selectedTerms),]
            } else {
              barPlotTable <- barPlotTable[-c(1:nrow(barPlotTable)),,drop=FALSE]
            }
            geneNum <- strsplit(barPlotTable$GeneRatio, split = "/")
            bgNum <- strsplit(barPlotTable$BgRatio, split = "/")
            newOverlap <- c()
            for(i in 1:length(geneNum)){
              newOverlap[i] <- paste(geneNum[[i]][1], "/", bgNum[[i]][1], sep = "")
            }
            barPlotTable$Overlap <- newOverlap
            if(input$rownameChoice2 == "termdesc"){
              term_vector <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
            }
            if(input$rownameChoice2 == "term"){
              term_vector <- barPlotTable$ID
            }
            if(input$rownameChoice2 == "desc"){
              term_vector <- barPlotTable$Description
              if(any(duplicated(term_vector))){
                idx <- which(duplicated(term_vector))
                terms <- term_vector[idx]
                for(i in 1:length(terms)){
                  idx2 <- which(term_vector == terms[i])
                  termTypes <- barPlotTable$termType[idx2]
                  for(j in 1:length(termTypes)){
                    if(termTypes[j] == "KEGG Pathways"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(KEGG)")
                    }
                    if(termTypes[j] == "Biological Processes"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOBP)")
                    }
                    if(termTypes[j] == "Cellular Components"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOCC)")
                    }
                    if(termTypes[j] == "Molecular Function"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOMF)")
                    }
                  }
                }
              }
            }
          }
          
          barPlotData <- data.frame("Category" = barPlotTable$termType, "Term" = term_vector, "Statistic" = log(barPlotTable[,termChoice]), "Overlap" = barPlotTable$Overlap)
          if(termChoice == "pvalue" | termChoice == "P.value"){
            ylabel <- "log(P-value)"
          }
          if(termChoice == "p.adjust" | termChoice == "Adjusted.P.value"){
            ylabel <- "log(Adj. P-value)"
          }
          if(reverseOrder){
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          } else {
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, -Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, -Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          }
          output$barPlot <- renderPlot({p}, height = input$barPlotHeight, width = input$barPlotWidth)
          barPlotHeight <- input$barPlotHeight*(2.54/100)
          barPlotWidth <- input$barPlotWidth*(2.54/100)
          output$downloadbarPlot <- downloadHandler(
            filename = function(){paste(input$barPlotFilename, '.png', sep = "")},
            content = function(file){ggsave(file, plot = p, device = "png", height = barPlotHeight, width = barPlotWidth, dpi = 100, units = "cm")}
          )
        } else {
          KEGG_table <- read.delim(fileList[grep(pattern = "KEGG", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          KEGG_table$termType <- c(rep("KEGG Pathways", nrow(KEGG_table)))
          GO_BP_table <- read.delim(fileList[grep(pattern = "GO_BP", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_BP_table$termType <- c(rep("Biological Processes", nrow(GO_BP_table)))
          GO_CC_table <- read.delim(fileList[grep(pattern = "GO_CC", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_CC_table$termType <- c(rep("Cellular Components", nrow(GO_CC_table)))
          GO_MF_table <- read.delim(fileList[grep(pattern = "GO_MF", x = fileList)], header = T, sep = "\t", stringsAsFactors = F)
          GO_MF_table$termType <- c(rep("Molecular Function", nrow(GO_MF_table)))
          termsToInclude <- input$combineGOKEGG
          barPlotTable <- matrix(ncol = 10)
          colnames(barPlotTable) <- column.names
          if("kegg" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, KEGG_table)
          }
          if("gobp" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_BP_table)
          }
          if("gocc" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_CC_table)
          }
          if("gomf" %in% termsToInclude){
            barPlotTable <- rbind(barPlotTable, GO_MF_table)
          }
          barPlotTable <- barPlotTable[-1,]
          barPlotTable <- barPlotTable[order(barPlotTable[,termChoice], decreasing = F),]
          
          selectedTerms <- input$selectInputTerms2
          
          if(length(grep("enrichR", fileList)) > 0){
            if(!is.null(selectedTerms)){
              barPlotTable <- barPlotTable[which(barPlotTable$Term %in% selectedTerms),]
            } else {
              barPlotTable <- barPlotTable[-c(1:nrow(barPlotTable)),,drop=FALSE]
            }
            term_vector <- barPlotTable$Term
          } else {
            if(!is.null(selectedTerms)){
              barPlotTable <- barPlotTable[which(barPlotTable$ID %in% selectedTerms),]
            } else {
              barPlotTable <- barPlotTable[-c(1:nrow(barPlotTable)),,drop=FALSE]
            }
            geneNum <- strsplit(barPlotTable$GeneRatio, split = "/")
            bgNum <- strsplit(barPlotTable$BgRatio, split = "/")
            newOverlap <- c()
            for(i in 1:length(geneNum)){
              newOverlap[i] <- paste(geneNum[[i]][1], "/", bgNum[[i]][1], sep = "")
            }
            barPlotTable$Overlap <- newOverlap
            if(input$rownameChoice2 == "termdesc"){
              term_vector <- paste(barPlotTable$ID, barPlotTable$Description, sep = "~")
            }
            if(input$rownameChoice2 == "term"){
              term_vector <- barPlotTable$ID
            }
            if(input$rownameChoice2 == "desc"){
              term_vector <- barPlotTable$Description
              if(any(duplicated(term_vector))){
                idx <- which(duplicated(term_vector))
                terms <- term_vector[idx]
                for(i in 1:length(terms)){
                  idx2 <- which(term_vector == terms[i])
                  termTypes <- barPlotTable$termType[idx2]
                  for(j in 1:length(termTypes)){
                    if(termTypes[j] == "KEGG Pathways"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(KEGG)")
                    }
                    if(termTypes[j] == "Biological Processes"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOBP)")
                    }
                    if(termTypes[j] == "Cellular Components"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOCC)")
                    }
                    if(termTypes[j] == "Molecular Function"){
                      term_vector[idx2[j]] <- paste(term_vector[idx2[j]], "(GOMF)")
                    }
                  }
                }
              }
            }
          }
          
          barPlotData <- data.frame("Category" = barPlotTable$termType, "Term" = term_vector, "Statistic" = log(barPlotTable[,termChoice]), "Overlap" = barPlotTable$Overlap)
          if(termChoice == "pvalue" | termChoice == "P.value"){
            ylabel <- "log(P-value)"
          }
          if(termChoice == "p.adjust" | termChoice == "Adjusted.P.value"){
            ylabel <- "log(Adj. P-value)"
          }
          if(reverseOrder){
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          } else {
            p <- ggplot(data = barPlotData, aes(x = reorder(Term, -Statistic), y = Statistic, fill = Category)) +
              geom_bar(stat = "identity", width = 0.8)+
              theme(panel.background=element_blank(),axis.line.x.bottom=element_line(),axis.line.y.left=element_line(),axis.ticks.y=element_blank(),axis.text=element_text(size=10))+
              ylab(ylabel)+xlab("Term")
            p <- p + coord_flip() + scale_y_continuous(trans="reverse", expand=c(0,0))
            if(addOverlapCount){
              p <- p + geom_text(data = barPlotData, 
                                 aes(x = reorder(Term, -Statistic), y = Statistic, label = Overlap), 
                                 hjust = +1.1, 
                                 size = 4,
                                 inherit.aes = TRUE)
            }
          }
          output$barPlot <- renderPlot({p}, height = input$barPlotHeight, width = input$barPlotWidth)
          barPlotHeight <- input$barPlotHeight*(2.54/100)
          barPlotWidth <- input$barPlotWidth*(2.54/100)
          output$downloadbarPlot <- downloadHandler(
            filename = function(){paste(input$barPlotFilename, '.png', sep = "")},
            content = function(file){ggsave(file, plot = p, device = "png", height = barPlotHeight, width = barPlotWidth, dpi = 100, units = "cm")}
          )
        }
      }
    }
  }
})