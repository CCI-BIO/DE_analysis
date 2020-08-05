## volcano_DEapp.R
## Script containing functions to allow the construction of volcano plots for the DE analysis application.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

output$selectedStatmethod <- renderUI({
  if(input$statSelect == "fdr"){
    numericInput("volcanoFDR", 
                 h4("FDR threshold:"),
                 min = 0, max = 1, step = 0.01,
                 value = 0.05)
  } else if(input$statSelect == "pval"){
    numericInput("volcanoPval", 
                 h4("P-value threshold:"), 
                 min = 0, max = 1, step = 0.01,
                 value = 0.05)
  }
})

output$selectedGeneListmethod <- renderUI({
  if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
    if(input$geneListOption == "topGenes"){
      numericInput("topGenesNumber",
                   h4("Number of top genes to label:"),
                   min = 1, max = max.row, step = 1,
                   value = 10)
    } else if(input$geneListOption == "selectInput"){
      selectizeInput("selectInputGenes", h4("Select genes to label:"), choices = NULL, multiple = TRUE)
    } else if(input$geneListOption == "listInput"){
      fileInput("file4", h4("List of genes (single column):"))
    }
  }
})

observeEvent({input$geneListOption},
             {
               if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                 tryCatch({volcano.data <- read.delim(input$fileChoice2, header = T, stringsAsFactors = F, sep = "\t")
                 updateSelectizeInput(session, "selectInputGenes", choices = rownames(volcano.data), server = TRUE)},
                 error=function(cond){return(NULL)})
               }
             })

# Generate volcano plot
observeEvent({input$refresh2},
             {
               Volcanorefresh <<- TRUE
               statSelect <- input$statSelect
               
               FCCut <- input$volcanoFC
               if(statSelect == "fdr"){
                 statType <- "FDR"
                 statCut <- input$volcanoFDR
               }
               if(statSelect == "pval"){
                 statType <- "PValue"
                 statCut <- input$volcanoPval
               }
               
               volcano.data <- read.delim(input$fileChoice2, header = T, stringsAsFactors = F, sep = "\t")
               volcano.data <- volcano.data[,c(statType, "logFC")]
               volcano.data$genes <- rownames(volcano.data)
               volcano.data$group <- c(rep("Not significant", nrow(volcano.data)))
               for(i in 1:nrow(volcano.data)){
                 if(abs(volcano.data$logFC[i]) >= log2(FCCut) & volcano.data[,c(statType)][i] > statCut){
                   volcano.data$group[i] <- "Sig. FC"
                 }
                 if(abs(volcano.data$logFC[i]) < log2(FCCut) & volcano.data[,c(statType)][i] <= statCut){
                   volcano.data$group[i] <- paste("Sig. ", statType, sep = "")
                 }
                 if(abs(volcano.data$logFC[i]) >= log2(FCCut) & volcano.data[,c(statType)][i] <= statCut){
                   volcano.data$group[i] <- paste("Sig. ", statType, " & Sig. FC", sep = "")
                 }
               }
               
               geneListOption <- input$geneListOption
               if(geneListOption == "topGenes"){
                 number_of_genes <- input$topGenesNumber
                 genes_to_label <- volcano.data$genes[1:number_of_genes]
               }
               if(geneListOption == "selectInput"){
                 genes_to_label <- input$selectInputGenes
               }
               if(geneListOption == "listInput"){
                 inFile <- input$file4
                 if(is.null(inFile)){
                   return(NULL)
                 }
                 conn <- file(inFile$datapath)
                 genes_to_label <- readLines(conn)
                 close(conn)
               }
               
               rownames(volcano.data) <- NULL
               colour.vector <- c("black", "red", "blue", "green")
               names(colour.vector) <- c("Not significant", paste("Sig. ", statType, sep = ""), "Sig. FC", paste("Sig. ", statType, " & Sig. FC", sep = ""))
               
               if(statType == "FDR"){
                 p <- ggplot(volcano.data, aes(logFC, -log10(FDR))) +
                   geom_point(aes(col = group)) +
                   scale_color_manual(values = colour.vector) + 
                   geom_text_repel(data = subset(volcano.data, genes%in%genes_to_label), aes(label = genes), min.segment.length = 0, box.padding = 0.5)
                 
               }
               if(statType == "PValue"){
                 p <- ggplot(volcano.data, aes(logFC, -log10(PValue))) +
                   geom_point(aes(col = group)) +
                   scale_color_manual(values = colour.vector) + 
                   geom_text_repel(data = subset(volcano.data, genes%in%genes_to_label), aes(label = genes), min.segment.length = 0, box.padding = 0.5)
                 
               }
               output$volcanoPlot <- renderPlot({p}, height = input$volcanoPlotHeight, width = input$volcanoPlotWidth)
               volcanoHeight <- input$volcanoPlotHeight*(2.54/100)
               volcanoWidth <- input$volcanoPlotWidth*(2.54/100)
               output$downloadVolcano <- downloadHandler(
                 filename = function(){paste(input$volcanoFilename, '.png', sep = "")},
                 content = function(file){ggsave(file, plot = p, device = "png", height = volcanoHeight, width = volcanoWidth, dpi = 100, units = "cm")}
               )
             }
)

shinyjs::disable("refresh2")

observeEvent(c(input$OutputDir, input$OutputDirNew, input$tabs, outputLocation$output),
             {
               if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                 tryCatch({volcano.data <- read.delim(outFileList()[[1]], header = T, stringsAsFactors = F, sep = "\t")
                 enable("refresh2")},
                 error=function(cond){
                   disable("refresh2")
                 }
                 )
               } else {
                 disable("refresh2")
               }
             }
)

shinyjs::disable("downloadVolcano")

observeEvent({input$volcanoFilename},{
  if(input$volcanoFilename == ""){
    disable("downloadVolcano")
  }
  else{
    if(Volcanorefresh == TRUE){
      enable("downloadVolcano")
    }
  }
})

observeEvent({input$refresh2},{
  if(nchar(input$volcanoFilename) >= 1){
    enable("downloadVolcano")
  }
})