## pathview_DEapp.R
## Script containing functions to allow the construction of KEGG pathview plots for the DE analysis application.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

# Get list of pathways for selection
observeEvent(
  c(input$OutputDir, input$OutputDirNew, input$tabs, input$fileChoice2, outputLocation$output, input$pathviewOrganism),
  {
    data(paths.hsa)
    names(paths.hsa) <- gsub(pattern = "hsa", replacement = "", x = names(paths.hsa))
    pathwayList <- names(paths.hsa)
    names(pathwayList) <- paths.hsa
    output$selectPathviewPath <- renderUI(
      {
        selectInput(
          "selectPathviewPath2",
          label = h4("Select KEGG pathway to plot:"),
          choices = pathwayList
        )
      }
    )
  }
)

# Select statistical threshold value
observeEvent(
  c(input$OutputDir, input$OutputDirNew, input$tabs, input$fileChoice2, outputLocation$output, input$pathviewFilterGeneList),
  {
    if(input$pathviewFilterGeneList == "pvalue"){
      output$pathviewThreshold <- renderUI(
        {
          numericInput(
            "pathviewThreshold2", 
            h4("P-value threshold:"), 
            min = 0, max = 1, step = 0.01,
            value = 0.05
          )
        }
      )
    } else if(input$pathviewFilterGeneList == "fdr"){
      output$pathviewThreshold <- renderUI(
        {
          numericInput(
            "pathviewThreshold2", 
            h4("FDR threshold:"), 
            min = 0, max = 1, step = 0.01,
            value = 0.05
          )
        }
      )
    }
  }
)

# Enable refresh button
observeEvent(
  c(input$OutputDir, input$OutputDirNew, input$tabs, outputLocation$output),
  {
    if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
      tryCatch(
        {
          test.data <- read.delim(outFileList()[[1]], header = T, stringsAsFactors = F, sep = "\t")
          enable("refresh5")
        },
        error = function(cond){
          disable("refresh5")
        }
      )
    } else {
      disable("refresh5")
    }
  }
)

# Generate pathview plot
observeEvent(
  c(input$refresh5),
  {
    if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
      shinyjs::disable("refresh5")
      geneList <- read.delim(input$fileChoice2, header = T, stringsAsFactors = F, sep = "\t")
      statThreshold <- input$pathviewThreshold2
      pathviewOrganism <- input$pathviewOrganism
      pathwayToTest <- input$selectPathviewPath2
      suffix <- input$pathviewSuffix
      if(is.null(suffix)){
        suffix <- "DE"
      }
      if(pathviewOrganism == "mmu"){
        pathwayToTest <- paste("mmu", pathwayToTest, sep = "")
      }
      if(pathviewOrganism == "hsa"){
        pathwayToTest <- paste("hsa", pathwayToTest, sep = "")
      }
      outDir <- outputLocation$output
      if(input$pathviewFilterGeneList == "pvalue"){
        geneList <- geneList[which(geneList$PValue <= statThreshold),]
        inputVector <- geneList$logFC
        names(inputVector) <- geneList$entrez_id
        prevDir <- getwd()
        setwd(outDir)
        pathview(gene.data = inputVector, pathway.id = pathwayToTest, species = pathviewOrganism, out.suffix = suffix)
        setwd(prevDir)
        shinyjs::enable("refresh5")
      }
      if(input$pathviewFilterGeneList == "fdr"){
        geneList <- geneList[which(geneList$FDR <= statThreshold),]
        inputVector <- geneList$logFC
        names(inputVector) <- geneList$entrez_id
        prevDir <- getwd()
        setwd(outDir)
        pathview(gene.data = inputVector, pathway.id = pathwayToTest, species = pathviewOrganism, out.suffix = suffix)
        setwd(prevDir)
        shinyjs::enable("refresh5")
      }
      output$pathviewPlot <- renderImage(
        {
          list(
            src = paste(outDir, "/", pathwayToTest, ".", suffix, ".png", sep = ""),
            contentType = "image/png",
            alt = "Plot"
          )
        }, 
        deleteFile = FALSE
      )
    }
  }
)