## heatmap_DEapp.R
## Script containing functions to allow the construction of heatmaps for the DE analysis application.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

shinyjs::disable("heatmapGeneSamplesOption")
shinyjs::disable("heatmapGeneListOption")

# List files
heatmapOutFileList <- eventReactive(c(input$OutputDir, input$OutputDirNew, input$tabs),
                                    {
                                      if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                                        outDir <- outputLocation$output
                                        fileList <- list.files(path = outDir, pattern = ".txt", full.names = T)
                                        fileList <- fileList[grep(pattern = "all_sample_tpm_DEanalysis.txt", x = fileList)]
                                        tryCatch({heatmap.data <- read.delim(fileList, header = T, stringsAsFactors = F, sep = "\t")
                                        updateSelectizeInput(session, "selectInputSamples", choices = colnames(heatmap.data), server = TRUE)
                                        enable("heatmapGeneSamplesOption")
                                        enable("heatmapGeneListOption")
                                        return(fileList)},
                                        error=function(cond){
                                          updateRadioButtons(session, "heatmapGeneSamplesOption", selected = "all")
                                          updateRadioButtons(session, "heatmapGeneListOption", selected = "topGenes")
                                          disable("heatmapGeneSamplesOption")
                                          disable("heatmapGeneListOption")
                                          return(NULL)})
                                      }
                                    }
)

output$heatmapSampleList <- renderUI({
  if(input$heatmapGeneSamplesOption == "select"){
    selectizeInput("selectInputSamples", h4("Select Samples to plot:"), choices = NULL, multiple = TRUE)
  }
})

observeEvent({input$heatmapGeneSamplesOption},
             {
               if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                 tryCatch({heatmap.data <- read.delim(heatmapOutFileList(), header = T, stringsAsFactors = F, sep = "\t")
                 updateSelectizeInput(session, "selectInputSamples", choices = colnames(heatmap.data), selected = colnames(heatmap.data), server = TRUE)},
                 error=function(cond){return(NULL)})
               }
             })

output$heatmapGeneList <- renderUI({
  if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
    if(input$heatmapGeneListOption == "topGenes"){
      numericInput("topGenesNumber2",
                   h4("Number of top genes to plot:"),
                   min = 1, max = max.row, step = 1,
                   value = 10)
    } else if(input$heatmapGeneListOption == "selectInput"){
      selectizeInput("selectInputGenes2", h4("Select genes to plot:"), choices = NULL, multiple = TRUE)
    } else if(input$heatmapGeneListOption == "listInput"){
      fileInput("file5", h4("List of genes (single column):"))
    }
  }
})

observeEvent({input$heatmapGeneListOption},
             {
               if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                 tryCatch({heatmap.data <- read.delim(heatmapOutFileList(), header = T, stringsAsFactors = F, sep = "\t")
                 updateSelectizeInput(session, "selectInputGenes2", choices = rownames(heatmap.data), server = TRUE)},
                 error=function(cond){return(NULL)})
               }
             })



observeEvent({input$file6},
             {
               differentialExpressionTest <- read.delim(heatmapOutFileList(), header = T, sep = "\t", stringsAsFactors = F)
               inFile6 <- input[["file6"]]
               if(is.null(inFile6)){
                 return(NULL)
               }
               if(!is.null(inFile6)){
                 test <- read.delim(inFile6$datapath, header = TRUE)
                 sample.names1 <- as.character(test$Samples)
                 sample.names2 <- as.character(colnames(differentialExpressionTest))
                 sample.names2 <- gsub(pattern = "\\.", replacement = "-", sample.names2)
                 if(all(sample.names2 %in% sample.names1)){
                   output$sampleMatch2 <- renderUI(HTML(paste("")))
                   heatmapGroupCheck <<- TRUE
                 } else {
                   output$sampleMatch2 <- renderUI(HTML(paste(tags$span(style="color:red", "Missing samples in group data!"))))
                   heatmapGroupCheck <<- FALSE
                 }
               }
             })

observeEvent({input$file6},
             {
               inFile6 <- input[["file6"]]
               if(is.null(inFile6)){
                 return(NULL)
               }
               if(heatmapGroupCheck == TRUE){
                 output$heatmapGroupOption <- renderUI({
                   radioButtons("heatmapGroupOption2",
                                h4("Goupings to label:"),
                                choices = list("No groups" = "none",
                                               "All groups" = "all",
                                               "Selected groups" = "selectGroups"),
                                selected = "none")
                 })
                 output$heatmapAverageOption <- renderUI({
                   radioButtons("heatmapAverageOption2",
                                h4("Average grouping:"),
                                choices = list("Yes" = "yes",
                                               "No" = "no"),
                                selected = "no")  
                 })
               } else {
                 output$heatmapGroupOption <- renderUI({})
                 output$heatmapAverageOption <- renderUI({})
               }
             })

observeEvent({input$heatmapGroupOption2}, {
  if(heatmapGroupCheck == TRUE){
    if(input$heatmapGroupOption2 == "none"){
      # heatmapGroupLastSelected <<- NULL
      enable("heatmapAverageOption2")
    }
    if(input$heatmapGroupOption2 == "all"){
      # heatmapGroupLastSelected <<- NULL
      updateRadioButtons(session, "heatmapAverageOption2", selected = "no")
      disable("heatmapAverageOption2")
    }
    if(input$heatmapGroupOption2 == "selectGroups"){
      inFile6 <- input[["file6"]]
      if(is.null(inFile6)){
        return(NULL)
      }
      test <- read.delim(inFile6$datapath, header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
      groupings <- colnames(test)
      groupings <- groupings[-which(groupings == "Samples")]
      names(groupings) <- groupings
      output$groupSelection <- renderUI({checkboxGroupInput("groupSelection2", h4("Select groupings to label:"), choices = as.list(groupings), selected = as.list(groupings))})
      updateRadioButtons(session, "heatmapAverageOption2", selected = "no")
      disable("heatmapAverageOption2")
    } else {
      output$groupSelection <- renderUI({})
    }
  } else {
    output$groupSelection <- renderUI({})
  }
})

observe({
  if(length(input$groupSelection2) == 1){
    heatmapGroupLastSelected <<- input$groupSelection2
  }
  if(length(input$groupSelection2) < 1){
    updateCheckboxGroupInput(session, "groupSelection2", selected = heatmapGroupLastSelected)
  }
})

observeEvent({input$heatmapAverageOption2}, {
  if(heatmapGroupCheck == TRUE){
    if(input$heatmapAverageOption2 == "yes"){
      inFile6 <- input[["file6"]]
      if(is.null(inFile6)){
        return(NULL)
      }
      test <- read.delim(inFile6$datapath, header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
      groupings <- colnames(test)
      groupings <- groupings[-which(groupings == "Samples")]
      names(groupings) <- groupings
      output$groupSelectionAverage <- renderUI({radioButtons("groupSelectionAverage2", h4("Select grouping to average:"), choices = as.list(groupings))})
    } else {
      output$groupSelectionAverage <- renderUI({})
    }
  } else {
    output$groupSelectionAverage <- renderUI({})
  }
})

# Generate heatmap
observeEvent({input$refresh3},
             {
               Heatmaprefresh <<- TRUE
               
               heatmapDistance <- input$distanceMetricCustom
               heatmapAgglomeration <- input$agglomerationMethodCustom
               
               heatmapScaling <- input$heatmapScaling
               heatmapRowNames <- input$heatmapRowNames
               heatmapColumnNames <- input$heatmapColumnNames
               
               if(heatmapRowNames == "yes"){
                 heatmapRowNames <- TRUE
               }
               if(heatmapRowNames == "no"){
                 heatmapRowNames <- FALSE
               }
               
               if(heatmapColumnNames == "yes"){
                 heatmapColumnNames <- TRUE
               }
               if(heatmapColumnNames == "no"){
                 heatmapColumnNames <- FALSE
               }
               
               heatmapTPM <- read.delim(heatmapOutFileList(), header = T, sep = "\t", stringsAsFactors = F)
               heatmapTPM <- as.matrix(heatmapTPM)
               mode(heatmapTPM) <- "numeric"
               heatmapTPM[heatmapTPM == 0] <- 0.0001
               
               differentialExpression <- read.delim(input$fileChoice2, header = T, sep = "\t", stringsAsFactors = F)
               
               heatmapGeneListOption <- input$heatmapGeneListOption
               if(heatmapGeneListOption == "topGenes"){
                 number_of_genes <- input$topGenesNumber2
                 genes_to_label <- rownames(differentialExpression)[1:number_of_genes]
               }
               if(heatmapGeneListOption == "selectInput"){
                 genes_to_label <- input$selectInputGenes2
               }
               if(heatmapGeneListOption == "listInput"){
                 inFile <- input$file5
                 if(is.null(inFile)){
                   return(NULL)
                 }
                 conn <- file(inFile$datapath)
                 genes_to_label <- readLines(conn)
                 close(conn)
               }
               
               if(is.null(genes_to_label)){
                 return(NULL)
               }
               else{
                 heatmapTPM <- heatmapTPM[which(rownames(heatmapTPM) %in% genes_to_label),]
               }
               
               if(input$heatmapGeneSamplesOption == "select"){
                 heatmapTPM <- heatmapTPM[,which(colnames(heatmapTPM) %in% input$selectInputSamples),drop = FALSE]
               }
               
               p <- Heatmap(matrix=heatmapTPM,
                            column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                            show_heatmap_legend = T,
                            clustering_distance_rows = heatmapDistance,
                            clustering_method_rows = heatmapAgglomeration,
                            clustering_distance_columns = heatmapDistance,
                            clustering_method_columns = heatmapAgglomeration,
                            column_names_gp = gpar(fontsize = 10),
                            row_names_gp = gpar(fontsize = 10),
                            show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
               
               if(heatmapGroupCheck == TRUE){
                 inFile6 <- input[["file6"]]
                 if(is.null(inFile6)){
                   return(NULL)
                 }
                 meta <- read.delim(inFile6$datapath, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
                 meta$Samples <- gsub(pattern = "-", replacement = ".", meta$Samples)
                 meta <- meta[which(meta$Samples %in% colnames(heatmapTPM)),]
                 
                 if(input$heatmapGroupOption2 == "none"){
                   if(input$heatmapAverageOption2 == "yes"){
                     groupToAverage <- input$groupSelectionAverage2
                     groupToAverageVector <- meta[,which(colnames(meta) == groupToAverage)]
                     groupToAverageVectorUnique <- unique(groupToAverageVector)
                     new.heatmapTPM <- matrix(nrow = nrow(heatmapTPM), ncol = length(groupToAverageVectorUnique))
                     rownames(new.heatmapTPM) <- rownames(heatmapTPM)
                     colnames(new.heatmapTPM) <- groupToAverageVectorUnique
                     for(i in 1:length(groupToAverageVectorUnique)){
                       groupToAverageVectorUnique.test <- groupToAverageVectorUnique[i]
                       heatmap.subset <- heatmapTPM[,which(groupToAverageVector == groupToAverageVectorUnique.test)]
                       heatmap.subset.average <- rowMeans(heatmap.subset)
                       new.heatmapTPM[,i] <- heatmap.subset.average
                     }
                     heatmapTPM <- new.heatmapTPM
                     
                     if(heatmapScaling == "row"){
                       heatmapTPM <- t(scale(t(log2(heatmapTPM))))
                     }
                     else{
                       heatmapTPM <- log2(heatmapTPM)
                     }
                     if(input$heatmapTranspose == "yes"){
                       heatmapTPM <- t(heatmapTPM)
                     }
                     p <- Heatmap(matrix=heatmapTPM,
                                  column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                                  show_heatmap_legend = T,
                                  clustering_distance_rows = heatmapDistance,
                                  clustering_method_rows = heatmapAgglomeration,
                                  clustering_distance_columns = heatmapDistance,
                                  clustering_method_columns = heatmapAgglomeration,
                                  column_names_gp = gpar(fontsize = 10),
                                  row_names_gp = gpar(fontsize = 10),
                                  show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
                     
                   } else {
                     if(heatmapScaling == "row"){
                       heatmapTPM <- t(scale(t(log2(heatmapTPM))))
                     }
                     else{
                       heatmapTPM <- log2(heatmapTPM)
                     }
                     if(input$heatmapTranspose == "yes"){
                       heatmapTPM <- t(heatmapTPM)
                     }
                     p <- Heatmap(matrix=heatmapTPM,
                                  column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                                  show_heatmap_legend = T,
                                  clustering_distance_rows = heatmapDistance,
                                  clustering_method_rows = heatmapAgglomeration,
                                  clustering_distance_columns = heatmapDistance,
                                  clustering_method_columns = heatmapAgglomeration,
                                  column_names_gp = gpar(fontsize = 10),
                                  row_names_gp = gpar(fontsize = 10),
                                  show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
                   }
                 }
                 
                 if(input$heatmapGroupOption2 == "all"){
                   # Set colours for ComplexHeatmap
                   complex_all_groups <- meta[,-which(colnames(meta) == "Samples"), drop = FALSE]
                   complex_groupList <- unique(as.vector(as.matrix(complex_all_groups)))
                   complex_ngroup <- length(complex_groupList)
                   complex_cols <- rainbow(complex_ngroup)
                   complex_cols <- as.data.frame(cbind(group = unique(complex_groupList), colour = complex_cols))
                   complex_all_colours <- list()
                   for(i in 1:ncol(complex_all_groups)){
                     complex_groupings <- as.character(unique(complex_all_groups[,i])) 
                     complex_cols_subset <- complex_cols[which(complex_cols$group %in% complex_groupings),]
                     complex_colGrp <- complex_cols_subset$colour
                     names(complex_colGrp) <- complex_cols_subset$group
                     complex_all_colours[[i]] <- complex_colGrp
                     groupName <- colnames(complex_all_groups)[i]
                     names(complex_all_colours)[i] <- groupName
                   }
                   
                   if(heatmapScaling == "row"){
                     heatmapTPM <- t(scale(t(log2(heatmapTPM))))
                   }
                   else{
                     heatmapTPM <- log2(heatmapTPM)
                   }
                   if(input$heatmapTranspose == "yes"){
                     heatmapTPM <- t(heatmapTPM)
                     ha1 <- HeatmapAnnotation(df = complex_all_groups, col = complex_all_colours, which = "row")
                     p <- ha1 + Heatmap(matrix=heatmapTPM,
                                        column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                                        show_heatmap_legend = T,
                                        clustering_distance_rows = heatmapDistance,
                                        clustering_method_rows = heatmapAgglomeration,
                                        clustering_distance_columns = heatmapDistance,
                                        clustering_method_columns = heatmapAgglomeration,
                                        column_names_gp = gpar(fontsize = 10),
                                        row_names_gp = gpar(fontsize = 10),
                                        show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
                   } else {
                     ha1 <- HeatmapAnnotation(df = complex_all_groups, col = complex_all_colours)
                     p <- Heatmap(matrix=heatmapTPM, top_annotation = ha1,
                                  column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                                  show_heatmap_legend = T,
                                  clustering_distance_rows = heatmapDistance,
                                  clustering_method_rows = heatmapAgglomeration,
                                  clustering_distance_columns = heatmapDistance,
                                  clustering_method_columns = heatmapAgglomeration,
                                  column_names_gp = gpar(fontsize = 10),
                                  row_names_gp = gpar(fontsize = 10),
                                  show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
                   }
                 }
                 
                 if(input$heatmapGroupOption2 == "selectGroups"){
                   # Set colours for ComplexHeatmap
                   complex_all_groups <- meta[,-which(colnames(meta) == "Samples"), drop = FALSE]
                   complex_all_groups <- complex_all_groups[,which(colnames(complex_all_groups) %in% input$groupSelection2), drop = FALSE]
                   complex_groupList <- unique(as.vector(as.matrix(complex_all_groups)))
                   complex_ngroup <- length(complex_groupList)
                   complex_cols <- rainbow(complex_ngroup)
                   complex_cols <- as.data.frame(cbind(group = unique(complex_groupList), colour = complex_cols))
                   complex_all_colours <- list()
                   for(i in 1:ncol(complex_all_groups)){
                     complex_groupings <- as.character(unique(complex_all_groups[,i])) 
                     complex_cols_subset <- complex_cols[which(complex_cols$group %in% complex_groupings),]
                     complex_colGrp <- complex_cols_subset$colour
                     names(complex_colGrp) <- complex_cols_subset$group
                     complex_all_colours[[i]] <- complex_colGrp
                     groupName <- colnames(complex_all_groups)[i]
                     names(complex_all_colours)[i] <- groupName
                   }
                   
                   if(heatmapScaling == "row"){
                     heatmapTPM <- t(scale(t(log2(heatmapTPM))))
                   }
                   else{
                     heatmapTPM <- log2(heatmapTPM)
                   }
                   if(input$heatmapTranspose == "yes"){
                     heatmapTPM <- t(heatmapTPM)
                     ha1 <- HeatmapAnnotation(df = complex_all_groups, col = complex_all_colours, which = "row")
                     p <- ha1 + Heatmap(matrix=heatmapTPM,
                                        column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                                        show_heatmap_legend = T,
                                        clustering_distance_rows = heatmapDistance,
                                        clustering_method_rows = heatmapAgglomeration,
                                        clustering_distance_columns = heatmapDistance,
                                        clustering_method_columns = heatmapAgglomeration,
                                        column_names_gp = gpar(fontsize = 10),
                                        row_names_gp = gpar(fontsize = 10),
                                        show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
                   } else {
                     ha1 <- HeatmapAnnotation(df = complex_all_groups, col = complex_all_colours)
                     p <- Heatmap(matrix=heatmapTPM, top_annotation = ha1,
                                  column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                                  show_heatmap_legend = T,
                                  clustering_distance_rows = heatmapDistance,
                                  clustering_method_rows = heatmapAgglomeration,
                                  clustering_distance_columns = heatmapDistance,
                                  clustering_method_columns = heatmapAgglomeration,
                                  column_names_gp = gpar(fontsize = 10),
                                  row_names_gp = gpar(fontsize = 10),
                                  show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
                   }
                 }
                 
               } else {
                 if(heatmapScaling == "row"){
                   heatmapTPM <- t(scale(t(log2(heatmapTPM))))
                 }
                 else{
                   heatmapTPM <- log2(heatmapTPM)
                 }
                 if(input$heatmapTranspose == "yes"){
                   heatmapTPM <- t(heatmapTPM)
                 }
                 p <- Heatmap(matrix=heatmapTPM,
                              column_title_gp = gpar(fontsize = 10, fontface = "bold"), name = "Scale",
                              show_heatmap_legend = T,
                              clustering_distance_rows = heatmapDistance,
                              clustering_method_rows = heatmapAgglomeration,
                              clustering_distance_columns = heatmapDistance,
                              clustering_method_columns = heatmapAgglomeration,
                              column_names_gp = gpar(fontsize = 10),
                              row_names_gp = gpar(fontsize = 10),
                              show_row_names = heatmapRowNames, show_column_names = heatmapColumnNames)
               }
               
               output$heatmapPlot <- renderPlot({p}, height = input$heatmapPlotHeight, width = input$heatmapPlotWidth)
               output$downloadHeatmap <- downloadHandler(
                 filename = function(){paste(input$heatmapFilename, '.png', sep = "")},
                 content = function(file){
                   png(file, width = input$heatmapPlotWidth, height = input$heatmapPlotHeight, res = 100)
                   print(p)
                   dev.off()
                 }
               )
             }
)

shinyjs::disable("refresh3")
shinyjs::disable("file6")

observeEvent(c(input$OutputDir, input$OutputDirNew, input$tabs, outputLocation$output),
             {
               if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                 tryCatch({heatmapTPM <- read.delim(heatmapOutFileList(), header = T, sep = "\t", stringsAsFactors = F)
                 enable("file6")
                 enable("refresh3")},
                 error=function(cond){
                   disable("file6")
                   disable("refresh3")
                 }
                 )
               } else {
                 disable("file6")
                 disable("refresh3")
               }
             }
)

shinyjs::disable("downloadHeatmap")

observeEvent({input$heatmapFilename},{
  if(input$heatmapFilename == ""){
    disable("downloadHeatmap")
  }
  else{
    if(Heatmaprefresh == TRUE){
      enable("downloadHeatmap")
    }
  }
})

observeEvent({input$refresh3},{
  if(nchar(input$heatmapFilename) >= 1){
    enable("downloadHeatmap")
  }
})