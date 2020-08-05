## PCA_DEapp.R
## Script containing functions to allow the construction of PCA plots for the DE analysis application.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

outNames <- reactive({
  inFile <- input[["file1"]]
  if(is.null(inFile)){
    return(NULL)
  }
  test <- read.delim(inFile$datapath, header = TRUE)
  if("length" %in% colnames(test)){
    test <- test[,-which(colnames(test) %in% "length")]
  }
  if("transcript_id.s." %in% colnames(test)){
    test <- test[,-which(colnames(test) %in% "transcript_id.s.")]
  }
  if("gene_id" %in% colnames(test)){
    test <- test[,-which(colnames(test) %in% "gene_id")]
  }
  vars <- rep(list(c()), ncol(test))
  names(vars) <- colnames(test)
  for(i in 1:length(vars)){
    vars[[i]] <- colnames(test)[i]
  }
  return(vars)
})

output$sampleChoice <- renderUI({
  selectInput("sampleChoice2", "Samples:", outNames(), selected = outNames(), multiple = TRUE)
})

outGroups <- reactive({
  inFile <- input[["file3"]]
  if(is.null(inFile)){
    return(NULL)
  }
  test <- read.delim(inFile$datapath, header = TRUE)
  vars <- rep(list(c()), (ncol(test)-1))
  names(vars) <- colnames(test)[2:ncol(test)]
  for(i in 1:length(vars)){
    vars[[i]] <- colnames(test)[i+1]
  }
  return(vars)
})

output$groupChoice <- renderUI({
  selectInput("groupChoice2", "Grouping to colour:", outGroups())
})

outColourBar <- reactive({
  inFile <- input[["file3"]]
  if(is.null(inFile)){
    return(NULL)
  }
  test <- read.delim(inFile$datapath, header = TRUE)
  if("Samples" %in% colnames(test)){
    test <- test[,-which(colnames(test) %in% "Samples"), drop = FALSE]
  }
  groupToCompare <- input$groupChoiceFinal2
  if(!is.null(groupToCompare)){
    if(groupToCompare %in% colnames(test)){
      test <- test[,-which(colnames(test) %in% groupToCompare), drop = FALSE]
    }
    if(ncol(test) == 0){
      return(NULL)
    } else if(ncol(test) > 0){
      vars <- rep(list(c()), ncol(test))
      names(vars) <- colnames(test)
      for(i in 1:length(vars)){
        vars[[i]] <- colnames(test)[i]
      }
      return(vars)
    }
  }
})

output$colourBar <- renderUI({
  checkboxGroupInput("colourBar2", h4("Additional colour bars:"), choices = outColourBar(), selected = outColourBar())
})

observeEvent(c(input$refresh1),
             {
               PCArefresh <<- TRUE
               samplesToKeep <- input$sampleChoice2
               inFile <- input[["file1"]]
               if(is.null(inFile)){
                 return(NULL)
               }
               test <- read.delim(inFile$datapath, header = TRUE)
               if("length" %in% colnames(test)){
                 test <- test[,-which(colnames(test) %in% "length")]
               }
               if("transcript_id.s." %in% colnames(test)){
                 test <- test[,-which(colnames(test) %in% "transcript_id.s.")]
               }
               if("gene_id" %in% colnames(test)){
                 test <- test[,-which(colnames(test) %in% "gene_id")]
               }
               test <- test[,which(colnames(test) %in% samplesToKeep)]
               test <- as.matrix(test)
               mode(test) <- "numeric"
               test[test == 0] <- 0.0001
               test <- log(test)
               test <- t(test)
               test <- test[,apply(test, 2, var) != 0]
               test.pca <- prcomp(test, center = TRUE, scale. = TRUE)
               inFile3 <- input[["file3"]]
               if(is.null(inFile3)){
                 return(NULL)
               }
               meta <- read.delim(inFile3$datapath, header = TRUE)
               rownames(meta) <- meta$Samples
               meta$Samples <- gsub(pattern = "-", replacement = ".", x = meta$Samples)
               meta <- meta[which(meta$Samples %in% samplesToKeep),]
               colourGroup <- input$groupChoice2
               output$pcaPlot <- renderPlot({autoplot(test.pca, label = TRUE, size = 0, data = meta, colour = colourGroup, label.show.legend = FALSE)}, height = 700)
               output$downloadPCA <- downloadHandler(
                 filename = function(){paste(input$PCAfilename, '.png', sep = "")},
                 content = function(file){ggsave(file, plot = autoplot(test.pca, label = TRUE, size = 0, data = meta, colour = colourGroup, label.show.legend = FALSE), device = "png", height = 18.5, width = 21.2, dpi = 100, units = "cm")}
               )
             })

shinyjs::disable("downloadPCA")

observeEvent({input$PCAfilename},{
  if(input$PCAfilename == ""){
    disable("downloadPCA")
  }
  else{
    if(PCArefresh == TRUE){
      inFile <- input[["file1"]]
      if(!is.null(inFile)){
        enable("downloadPCA")
      }
    }
  }
})

observeEvent({input$refresh1},{
  if(nchar(input$PCAfilename) >= 1){
    enable("downloadPCA")
  }
})