## main_DEapp.R
## Script containing functions for the edgeR DE analysis.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

## UI commands
output$groupChoiceFinal <- renderUI({
  selectInput("groupChoiceFinal2", h4("Groupings to compare:"), outGroups())
})

outComb <- reactive({
  if(is.null(outGroups())){
    return(NULL)
  }
  samplesToKeep <- input$sampleChoice2
  groupToTest <- input$groupChoiceFinal2
  inFile3 <- input[["file3"]]
  if(is.null(inFile3)){
    return(NULL)
  }
  meta <- read.delim(inFile3$datapath, header = TRUE)
  
  rownames(meta) <- meta$Samples
  meta$Samples <- gsub(pattern = "-", replacement = ".", x = meta$Samples)
  meta <- meta[which(meta$Samples %in% samplesToKeep),]
  group <- meta[,groupToTest]
  
  tryCatch(expr = {
    design <- model.matrix(~ -1+factor(group))
    colnames(design) <- unique(group)
    
    design.pairs <- function(levels){
      n <- length(levels)
      design <- matrix(0, n, choose(n, 2))
      rownames(design) <- levels
      colnames(design) <- 1:choose(n,2)
      k <- 0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          k <- k+1
          design[i,k] <- 1
          design[j,k] <- -1
          colnames(design)[k] <- paste(levels[i], "-", levels[j], sep = "")
        }
      }
      design
    }
    contrast.matrix <- design.pairs(unique(group))
    vars <- rep(list(c()), ncol(contrast.matrix))
    names(vars) <- colnames(contrast.matrix)
    for(i in 1:length(vars)){
      vars[[i]] <- colnames(contrast.matrix)[i]
    }
    return(vars)},
    error = function(e){return(NULL)},
    warning = function(w){return(NULL)})
})

output$combChoice1 <- renderUI({
  checkboxGroupInput(inputId = "combChoice2", label = h4("Comparisons:"), choices = outComb(), selected = outComb())
})

observe({
  if(length(input$combChoice2) == 1){
    combinationLastSelected <<- input$combChoice2
  }
  if(length(input$combChoice2) < 1){
    updateCheckboxGroupInput(session, "combChoice2", selected = combinationLastSelected)
  }
})

observeEvent(c(input$offlineMode, input$tabs),{
  if(input$offlineMode == "david"){
    output$davidEmail <- renderUI({textInput("davidEmail2", label = h4("DAVID registered email:"))}) 
    output$internalKEGG <- renderUI({})
    output$statisticToUse <- renderUI({radioButtons("statisticToUse2",
                                                    h4("Statistic for filtering input for GO/KEGG enrichment:"),
                                                    choices = list("P-value" = "PValue",
                                                                   "FDR" = "FDR"))})
    output$separateDirections <- renderUI({radioButtons("separateDirections2",
                                                        h4("Split UP and DOWN regulated for GO/KEGG enrichment:"),
                                                        choices = list("Yes" = "yes",
                                                                       "No" = "no"),
                                                        selected = "yes")})
  }
  if(input$offlineMode == "clusterprofiler"){
    output$davidEmail <- renderUI({})
    output$internalKEGG <- renderUI({radioButtons("internalKEGG2", label = h4("Use KEGG offline database (not recommended):"), 
                                                  choices = list("Yes" = "yes",
                                                                 "No" = "no"),
                                                  selected = "no")})
    output$statisticToUse <- renderUI({radioButtons("statisticToUse2",
                                                    h4("Statistic for filtering input for GO/KEGG enrichment:"),
                                                    choices = list("P-value" = "PValue",
                                                                   "FDR" = "FDR"))})
    output$separateDirections <- renderUI({radioButtons("separateDirections2",
                                                        h4("Split UP and DOWN regulated for GO/KEGG enrichment:"),
                                                        choices = list("Yes" = "yes",
                                                                       "No" = "no"),
                                                        selected = "yes")})
  }
  if(input$offlineMode == "enrichr"){
    output$davidEmail <- renderUI({})
    output$internalKEGG <- renderUI({})
    output$statisticToUse <- renderUI({radioButtons("statisticToUse2",
                                                    h4("Statistic for filtering input for GO/KEGG enrichment:"),
                                                    choices = list("P-value" = "PValue",
                                                                   "FDR" = "FDR"))})
    output$separateDirections <- renderUI({radioButtons("separateDirections2",
                                                        h4("Split UP and DOWN regulated for GO/KEGG enrichment:"),
                                                        choices = list("Yes" = "yes",
                                                                       "No" = "no"),
                                                        selected = "yes")})
  }
  if(input$offlineMode == "none"){
    output$davidEmail <- renderUI({})
    output$internalKEGG <- renderUI({})
    output$statisticToUse <- renderUI({})
    output$separateDirections <- renderUI({})
  }
})

observeEvent(c(input$offlineMode, input$davidEmail2), {
  if(input$offlineMode == "david"){
    email <- input$davidEmail2
    valid_url <- paste("https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=", email, sep = "")
    webpage <- GET(valid_url)
    webpage <- httr::content(webpage, as = "text")
    if(length(grep(pattern = "true", webpage)) == 1){
      output$emailWarn1 <- renderText("Valid email!")
      output$emailWarn2 <- NULL
      emailCheck <<- TRUE
    } else {
      output$emailWarn1 <- NULL
      output$emailWarn2 <- renderText("Unregistered email, please register with DAVID!")
      emailCheck <<- FALSE
    }
  } else{
    output$emailWarn1 <- NULL
    output$emailWarn2 <- NULL
    emailCheck <<- FALSE
  }
})

shinyjs::disable("submit1")

## Main algorithm

observeEvent(c(input$davidEmail2, input$combChoice2, input$offlineMode, input$OutputDir), {
  if(emailCheck == TRUE){
    if(!is.null(input$combChoice2)){
      OutputDir <- input$OutputDir
      if(class(OutputDir[1]) == "list"){
        outDir <- paste(OutputDir$root, paste(unlist(OutputDir$path[-1]), collapse = "/"))
        if(nchar(outDir) > 1){
          shinyjs::enable("submit1")
        } else {
          shinyjs::disable("submit1")
        }
      } else {
        shinyjs::disable("submit1")
      }
    } else {
      shinyjs::disable("submit1")
    }
  } 
  if(emailCheck == FALSE){
    if(input$offlineMode != "david"){
      if(!is.null(input$combChoice2)){
        OutputDir <- input$OutputDir
        if(class(OutputDir[1]) == "list"){
          outDir <- paste(OutputDir$root, paste(unlist(OutputDir$path[-1]), collapse = "/"))
          if(nchar(outDir) > 1){
            shinyjs::enable("submit1")
          } else {
            shinyjs::disable("submit1")
          }
        } else {
          shinyjs::disable("submit1")
        }
      } else {
        shinyjs::disable("submit1")
      }
    } else {
      shinyjs::disable("submit1")
    }
  }
})

observeEvent(input$submit1,
             {withProgress(message = "Running Differential Expression", value = 0, 
                           { incProgress(amount = 0, message = "Preparing data")
                             
                             # Set working directory
                             OutputDir <- input$OutputDir
                             setwd(paste(OutputDir$root, paste(unlist(OutputDir$path[-1]), collapse = "/"), sep = ""))
                             
                             # Get parameters for analysis
                             cpmCut <- as.numeric(input$num1)
                             percentCut <- (as.numeric(input$slider1)/100)
                             FCCut <- as.numeric(input$FCgroup)
                             PvalCut <- input$PvalInput
                             FDRCut <- input$FDRInput
                             
                             clusterMethod <- input$agglomerationMethod
                             distanceMethod <- input$distanceMetric
                             rowNameCheck <- input$plotRowNames
                             zScoreScaling <- input$zScoreScaling
                             
                             colourBars <- input$colourBar2
                             
                             genome <- input$Organism
                             if(genome == "human"){
                               genomeUse <- org.Hs.egSYMBOL
                               orgEnrichment <- "hsa"
                               genomeFull <- "Homo sapiens"
                             }
                             if(genome == "mouse"){
                               genomeUse <- org.Mm.egSYMBOL
                               orgEnrichment <- "mmu"
                               genomeFull <- "Mus musculus"
                             }
                             
                             # Get expression and meta data
                             samplesToKeep <- input$sampleChoice2
                             
                             counts <- read.delim(isolate(input$file1$datapath), stringsAsFactors = FALSE, header = TRUE, row.names = 1)
                             if("length" %in% colnames(counts)){
                               counts <- counts[,-which(colnames(counts) %in% "length")]
                             }
                             counts <- counts[, which(colnames(counts) %in% samplesToKeep)]
                             counts <- as.matrix(counts)
                             mode(counts) <- "numeric"
                             countsToPCA <- counts
                             tpm <- read.delim(isolate(input$file2$datapath), stringsAsFactors = FALSE, header = TRUE, row.names = 1)
                             if("transcript_id.s." %in% colnames(tpm)){
                               tpm <- tpm[,-which(colnames(tpm) %in% "transcript_id.s.")]
                             }
                             tpm <- tpm[, which(colnames(tpm) %in% samplesToKeep)]
                             tpm <- as.matrix(tpm)
                             mode(tpm) <- "numeric"
                             meta <- read.delim(isolate(input$file3$datapath), stringsAsFactors = FALSE, header = TRUE)
                             meta$Samples <- gsub(pattern = "-", replacement = ".", x = meta$Samples)
                             meta <- meta[which(meta$Samples %in% samplesToKeep),]
                             
                             # Get groups
                             groupToTest <- input$groupChoiceFinal2
                             group <- meta[,groupToTest]
                             
                             # Set colours for ComplexHeatmap
                             complex_all_groups <- meta[,-which(colnames(meta) == "Samples"), drop = FALSE]
                             if(is.null(colourBars)){
                               complex_all_groups <- complex_all_groups[,which(colnames(complex_all_groups) == groupToTest), drop = FALSE]
                             }
                             if(!is.null(colourBars)){
                               complex_all_groups <- complex_all_groups[,which(colnames(complex_all_groups) %in% c(groupToTest, colourBars)), drop = FALSE]
                             }
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
                             ha1 <- HeatmapAnnotation(df = complex_all_groups, col = complex_all_colours)
                             
                             complex_test_group <- unique(complex_all_groups[,groupToTest, drop = FALSE])
                             complex_test_cols_subset <- complex_cols[which(complex_cols$group %in% complex_test_group[,1]),]
                             complex_test_colours <- complex_test_cols_subset$colour
                             names(complex_test_colours) <- complex_test_cols_subset$group
                             complex_test_colours <- list(complex_test_colours)
                             names(complex_test_colours) <- groupToTest
                             ha2 <- HeatmapAnnotation(df = complex_test_group, col = complex_test_colours)
                             
                             # Filter genes
                             keep <- rowSums(cpm(counts) > cpmCut) >= (percentCut*ncol(counts)+1)
                             counts <- counts[keep,]
                             tpm <- tpm[keep,]
                             
                             write.table(tpm, "all_sample_tpm_DEanalysis.txt", quote = F, sep = "\t")
                             
                             # Run DGE
                             y <- DGEList(counts = counts, group = group)
                             y <- calcNormFactors(y)
                             
                             incProgress(amount = 1/3, message = "Creating QC plots")
                             
                             # Plot PCA from Step 2
                             test <- countsToPCA
                             test[test == 0] <- 0.0001
                             test <- log(test)
                             test <- t(test)
                             test <- test[,apply(test, 2, var) != 0]
                             test.pca <- prcomp(test, center = TRUE, scale. = TRUE)
                             metaToPCA <- meta
                             rownames(metaToPCA) <- metaToPCA$Samples
                             ggsave(filename = paste(groupToTest, "_PCA_plot.png", sep = ""), plot = autoplot(test.pca, label = TRUE, size = 0, data = metaToPCA, colour = groupToTest, label.show.legend = FALSE), device = "png", height = 18.5, width = 21.2, dpi = 100, units = "cm")
                             
                             # Plot QC metrics
                             pdf(file = "QC_metrics_plots.pdf")
                             par(mar = c(5.1,4.1,4.1,10))
                             plotMDS(y, main = "MDS plot of RNA-seq data")
                             plotMD(y, column = 1)
                             y <- estimateCommonDisp(y)
                             y <- estimateGLMTrendedDisp(y)
                             y <- estimateTagwiseDisp(y)
                             plotBCV(y)
                             dev.off()
                             
                             design <- model.matrix(~ -1+factor(group))
                             colnames(design) <- unique(group)
                             
                             design.pairs <- function(levels){
                               n <- length(levels)
                               design <- matrix(0, n, choose(n, 2))
                               rownames(design) <- levels
                               colnames(design) <- 1:choose(n,2)
                               k <- 0
                               for(i in 1:(n-1)){
                                 for(j in (i+1):n){
                                   k <- k+1
                                   design[i,k] <- 1
                                   design[j,k] <- -1
                                   colnames(design)[k] <- paste(levels[i], "-", levels[j], sep = "")
                                 }
                               }
                               design
                             }
                             
                             contrast.matrix <- design.pairs(unique(group))
                             
                             combToKeep <- input$combChoice2
                             
                             contrast.matrix <- contrast.matrix[,which(colnames(contrast.matrix) %in% combToKeep), drop = FALSE]
                             
                             fit <- glmFit(y, design)
                             egSymbol <- toTable(genomeUse)
                             
                             incProgress(amount = 1/3, message = "Running DGE comparisons")
                             
                             withProgress(message = "Comparing groups", value = 0,
                                          {for(i in 1:ncol(contrast.matrix)){
                                            comparisonName <- colnames(contrast.matrix)[i]
                                            incProgress(1/ncol(contrast.matrix), detail = paste("Running comparison: ", comparisonName, sep = ""))
                                            deg_lrt <- glmLRT(fit, contrast = contrast.matrix[,i])
                                            deg_all <- topTags(deg_lrt, n = nrow(deg_lrt), adjust.method = "BH", sort.by = "PValue")
                                            deg_all_table <- as.data.frame(deg_all)
                                            FC <- logratio2foldchange(deg_all_table$logFC, base = 2)
                                            deg_expTab <- cbind.data.frame(deg_all_table, FC)
                                            m <- match(rownames(deg_expTab), egSymbol$symbol)
                                            deg_expTab$entrez_id <- egSymbol$gene_id[m]
                                            fileName <- paste(comparisonName, "_DEanalysis.txt", sep = "")
                                            write.table(deg_expTab, fileName, sep = "\t", quote = F)
                                            temp <- merge(deg_expTab, tpm, "row.names")
                                            reOrdCol <- c(1, 8, 2:7, 9:ncol(temp))
                                            deg_expTab_tpm <- temp[order(temp$PValue), reOrdCol]
                                            fileName <- paste(comparisonName, "_tpm_DEanalysis.txt", sep = "")
                                            write.table(deg_expTab_tpm, fileName, sep = "\t", row.names = F, quote = F)
                                            
                                            # GSEA input files 
                                            gsea_expression_table <- merge(deg_expTab, tpm, "row.names")
                                            gsea_expression_table <- gsea_expression_table[order(gsea_expression_table$PValue),]
                                            compare1 <- strsplit(comparisonName, "-")[[1]][1]
                                            compare2 <- strsplit(comparisonName, "-")[[1]][2]
                                            comparison_idx <- c(which(design[,compare1] == 1), which(design[,compare2] == 1))
                                            comparison_idx <- comparison_idx + 8
                                            desc_column <- c(rep("na", nrow(gsea_expression_table)))
                                            gsea_expression_table_final <- cbind("NAME" = gsea_expression_table[,1], "DESCRIPTION" = desc_column, gsea_expression_table[,c(comparison_idx)])
                                            fileName <- paste(comparisonName, "_GSEA_expression_input.txt", sep = "")
                                            write.table(gsea_expression_table_final, fileName, sep = "\t", row.names = F, quote = F)
                                            
                                            num.samples <- length(comparison_idx)
                                            num.classes <- 2
                                            line1 <- paste(num.samples, num.classes, "1", sep = "\t")
                                            line2 <- paste("#", compare1, compare2, sep = "\t")
                                            line3 <- paste(paste(c(rep("0", length(which(design[,compare1] == 1)))), collapse = "\t"), paste(c(rep("1", length(which(design[,compare2] == 1)))), collapse = "\t"), sep = "\t")
                                            writeLines(c(line1, line2, line3), con = paste(comparisonName, "_GSEA_phenotype_input.cls", sep = ""))
                                            
                                            gsea_expression_table$fcsign <- sign(as.numeric(gsea_expression_table$logFC))
                                            gsea_expression_table$logP <- -log10(gsea_expression_table$PValue)
                                            gsea_expression_table$metric <- gsea_expression_table$logP/gsea_expression_table$fcsign
                                            gsea_expression_table <- gsea_expression_table[order(gsea_expression_table$metric, decreasing = T),]
                                            fileName <- paste(comparisonName, "_GSEA_rank_input.rnk", sep = "")
                                            write.table(gsea_expression_table[,c("Row.names", "metric")], fileName, sep = "\t", row.names = F, quote = F)
                                            
                                            # KEGG/GO enrichment
                                            offlineMode <- input$offlineMode
                                            
                                            ####################################
                                            # CLUSTERPROFILER analysis
                                            ####################################
                                            
                                            if(offlineMode == "clusterprofiler"){
                                              # Run clusterProfiler GO and KEGG analysis and GSEA
                                              statisticToUse <- input$statisticToUse2
                                              if(statisticToUse == "PValue"){
                                                thresholdToUse <- PvalCut
                                              } else {
                                                thresholdToUse <- FDRCut
                                              }
                                              internalKEGG <- input$internalKEGG2
                                              if(internalKEGG == "yes"){
                                                internalKEGG <- TRUE
                                              } else {
                                                internalKEGG <- FALSE
                                              }
                                              separateDirections <- input$separateDirections2
                                              if(separateDirections == "yes"){
                                                # DOWN REGULATED
                                                down_filtered_list <- deg_expTab$entrez_id[which(deg_expTab[,statisticToUse] < thresholdToUse & deg_expTab[,"FC"] <= -FCCut)]
                                                down_filtered_list <- na.exclude(down_filtered_list)
                                                
                                                if(length(down_filtered_list) > 1){
                                                  #KEGG
                                                  DOWN_KEGG_enrichment <- enrichKEGG(gene = down_filtered_list, organism = orgEnrichment, use_internal_data = internalKEGG)
                                                  if(is.null(DOWN_KEGG_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    DOWN_KEGG_enrichment <- data.frame(DOWN_KEGG_enrichment)
                                                    ID_vec <- DOWN_KEGG_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_KEGG_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_KEGG_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  #GO
                                                  if(orgEnrichment == "hsa"){
                                                    DOWN_BP_GO_enrichment <- enrichGO(gene = down_filtered_list, OrgDb = org.Hs.eg.db, ont = "BP")
                                                    DOWN_CC_GO_enrichment <- enrichGO(gene = down_filtered_list, OrgDb = org.Hs.eg.db, ont = "CC")
                                                    DOWN_MF_GO_enrichment <- enrichGO(gene = down_filtered_list, OrgDb = org.Hs.eg.db, ont = "MF")
                                                  }
                                                  if(orgEnrichment == "mmu"){
                                                    DOWN_BP_GO_enrichment <- enrichGO(gene = down_filtered_list, OrgDb = org.Mm.eg.db, ont = "BP")
                                                    DOWN_CC_GO_enrichment <- enrichGO(gene = down_filtered_list, OrgDb = org.Mm.eg.db, ont = "CC")
                                                    DOWN_MF_GO_enrichment <- enrichGO(gene = down_filtered_list, OrgDb = org.Mm.eg.db, ont = "MF")
                                                  }
                                                  
                                                  if(is.null(DOWN_BP_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    DOWN_BP_GO_enrichment <- data.frame(DOWN_BP_GO_enrichment)
                                                    ID_vec <- DOWN_BP_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_BP_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_BP_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  if(is.null(DOWN_CC_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    DOWN_CC_GO_enrichment <- data.frame(DOWN_CC_GO_enrichment)
                                                    ID_vec <- DOWN_CC_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_CC_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_CC_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  if(is.null(DOWN_MF_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    DOWN_MF_GO_enrichment <- data.frame(DOWN_MF_GO_enrichment)
                                                    ID_vec <- DOWN_MF_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_MF_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_MF_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                }
                                                
                                                # UP REGULATED
                                                up_filtered_list <- deg_expTab$entrez_id[which(deg_expTab[,statisticToUse] < thresholdToUse & deg_expTab[,"FC"] >= FCCut)]
                                                up_filtered_list <- na.exclude(up_filtered_list)
                                                
                                                if(length(up_filtered_list) > 1){
                                                  #KEGG
                                                  UP_KEGG_enrichment <- enrichKEGG(gene = up_filtered_list, organism = orgEnrichment, use_internal_data = internalKEGG)
                                                  if(is.null(UP_KEGG_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    UP_KEGG_enrichment <- data.frame(UP_KEGG_enrichment)
                                                    ID_vec <- UP_KEGG_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_KEGG_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_KEGG_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  #GO
                                                  if(orgEnrichment == "hsa"){
                                                    UP_BP_GO_enrichment <- enrichGO(gene = up_filtered_list, OrgDb = org.Hs.eg.db, ont = "BP")
                                                    UP_CC_GO_enrichment <- enrichGO(gene = up_filtered_list, OrgDb = org.Hs.eg.db, ont = "CC")
                                                    UP_MF_GO_enrichment <- enrichGO(gene = up_filtered_list, OrgDb = org.Hs.eg.db, ont = "MF")
                                                  }
                                                  if(orgEnrichment == "mmu"){
                                                    UP_BP_GO_enrichment <- enrichGO(gene = up_filtered_list, OrgDb = org.Mm.eg.db, ont = "BP")
                                                    UP_CC_GO_enrichment <- enrichGO(gene = up_filtered_list, OrgDb = org.Mm.eg.db, ont = "CC")
                                                    UP_MF_GO_enrichment <- enrichGO(gene = up_filtered_list, OrgDb = org.Mm.eg.db, ont = "MF")
                                                  }
                                                  
                                                  if(is.null(UP_BP_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    UP_BP_GO_enrichment <- data.frame(UP_BP_GO_enrichment)
                                                    ID_vec <- UP_BP_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_BP_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_BP_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  if(is.null(UP_CC_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    UP_CC_GO_enrichment <- data.frame(UP_CC_GO_enrichment)
                                                    ID_vec <- UP_CC_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_CC_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_CC_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  if(is.null(UP_MF_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    UP_MF_GO_enrichment <- data.frame(UP_MF_GO_enrichment)
                                                    ID_vec <- UP_MF_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_MF_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_MF_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                }
                                              } else {
                                                # COMBINED
                                                filtered_list <- deg_expTab$entrez_id[which(deg_expTab[,statisticToUse] < thresholdToUse & abs(deg_expTab[,"FC"]) >= FCCut)]
                                                filtered_list <- na.exclude(filtered_list)
                                                
                                                if(length(filtered_list) > 1){
                                                  #KEGG
                                                  KEGG_enrichment <- enrichKEGG(gene = filtered_list, organism = orgEnrichment, use_internal_data = internalKEGG)
                                                  if(is.null(KEGG_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    KEGG_enrichment <- data.frame(KEGG_enrichment)
                                                    ID_vec <- KEGG_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    KEGG_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(KEGG_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  #GO
                                                  if(orgEnrichment == "hsa"){
                                                    BP_GO_enrichment <- enrichGO(gene = filtered_list, OrgDb = org.Hs.eg.db, ont = "BP")
                                                    CC_GO_enrichment <- enrichGO(gene = filtered_list, OrgDb = org.Hs.eg.db, ont = "CC")
                                                    MF_GO_enrichment <- enrichGO(gene = filtered_list, OrgDb = org.Hs.eg.db, ont = "MF")
                                                  }
                                                  if(orgEnrichment == "mmu"){
                                                    BP_GO_enrichment <- enrichGO(gene = filtered_list, OrgDb = org.Mm.eg.db, ont = "BP")
                                                    CC_GO_enrichment <- enrichGO(gene = filtered_list, OrgDb = org.Mm.eg.db, ont = "CC")
                                                    MF_GO_enrichment <- enrichGO(gene = filtered_list, OrgDb = org.Mm.eg.db, ont = "MF")
                                                  }
                                                  
                                                  if(is.null(BP_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    BP_GO_enrichment <- data.frame(BP_GO_enrichment)
                                                    ID_vec <- BP_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    BP_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(BP_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                    
                                                  }
                                                  
                                                  if(is.null(CC_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    CC_GO_enrichment <- data.frame(CC_GO_enrichment)
                                                    ID_vec <- CC_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    CC_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(CC_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                  
                                                  if(is.null(MF_GO_enrichment)){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  } else {
                                                    MF_GO_enrichment <- data.frame(MF_GO_enrichment)
                                                    ID_vec <- MF_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    MF_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(MF_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  }
                                                }
                                              }
                                            }
                                            
                                            ####################################
                                            # DAVID analysis
                                            ####################################
                                            
                                            if(offlineMode == "david"){
                                              statisticToUse <- input$statisticToUse2
                                              if(statisticToUse == "PValue"){
                                                thresholdToUse <- PvalCut
                                              } else {
                                                thresholdToUse <- FDRCut
                                              }
                                              separateDirections <- input$separateDirections2
                                              if(separateDirections == "yes"){
                                                # DOWN REGULATED
                                                down_filtered_list <- deg_expTab$entrez_id[which(deg_expTab[,statisticToUse] < thresholdToUse & deg_expTab[,"FC"] <= -FCCut)]
                                                down_filtered_list <- na.exclude(down_filtered_list)
                                                if(length(down_filtered_list) > 1){
                                                  #KEGG
                                                  tryCatch(expr = {
                                                    DOWN_KEGG_enrichment <- enrichDAVID(gene = down_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    DOWN_KEGG_enrichment <- data.frame(DOWN_KEGG_enrichment)
                                                    ID_vec <- DOWN_KEGG_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_KEGG_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_KEGG_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  #GO
                                                  tryCatch(expr = {
                                                    DOWN_BP_GO_enrichment <- enrichDAVID(gene = down_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    DOWN_BP_GO_enrichment <- data.frame(DOWN_BP_GO_enrichment)
                                                    ID_vec <- DOWN_BP_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_BP_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_BP_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  tryCatch(expr = {
                                                    DOWN_CC_GO_enrichment <- enrichDAVID(gene = down_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    DOWN_CC_GO_enrichment <- data.frame(DOWN_CC_GO_enrichment)
                                                    ID_vec <- DOWN_CC_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_CC_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_CC_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  tryCatch(expr = {
                                                    DOWN_MF_GO_enrichment <- enrichDAVID(gene = down_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    DOWN_MF_GO_enrichment <- data.frame(DOWN_MF_GO_enrichment)
                                                    ID_vec <- DOWN_MF_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    DOWN_MF_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(DOWN_MF_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DOWN_DAVID_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                }
                                                
                                                # UP REGULATED
                                                up_filtered_list <- deg_expTab$entrez_id[which(deg_expTab[,statisticToUse] < thresholdToUse & deg_expTab[,"FC"] >= FCCut)]
                                                up_filtered_list <- na.exclude(up_filtered_list)
                                                
                                                if(length(up_filtered_list) > 1){
                                                  #KEGG
                                                  tryCatch(expr = {
                                                    UP_KEGG_enrichment <- enrichDAVID(gene = up_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    UP_KEGG_enrichment <- data.frame(UP_KEGG_enrichment)
                                                    ID_vec <- UP_KEGG_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_KEGG_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_DAVID_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_KEGG_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_DAVID_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  #GO
                                                  tryCatch(expr = {
                                                    UP_BP_GO_enrichment <- enrichDAVID(gene = up_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    UP_BP_GO_enrichment <- data.frame(UP_BP_GO_enrichment)
                                                    ID_vec <- UP_BP_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_BP_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_DAVID_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_BP_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_DAVID_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  tryCatch(expr = {
                                                    UP_CC_GO_enrichment <- enrichDAVID(gene = up_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    UP_CC_GO_enrichment <- data.frame(UP_CC_GO_enrichment)
                                                    ID_vec <- UP_CC_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_CC_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_DAVID_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_CC_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_DAVID_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  tryCatch(expr = {
                                                    UP_MF_GO_enrichment <- enrichDAVID(gene = up_filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    UP_MF_GO_enrichment <- data.frame(UP_MF_GO_enrichment)
                                                    ID_vec <- UP_MF_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    UP_MF_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_UP_DAVID_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(UP_MF_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_UP_DAVID_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                }
                                              } else {
                                                # COMBINED
                                                filtered_list <- deg_expTab$entrez_id[which(deg_expTab[,statisticToUse] < thresholdToUse & abs(deg_expTab[,"FC"]) >= FCCut)]
                                                filtered_list <- na.exclude(filtered_list)
                                                
                                                if(length(filtered_list) > 1){
                                                  #KEGG
                                                  tryCatch(expr = {
                                                    KEGG_enrichment <- enrichDAVID(gene = filtered_list, idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    KEGG_enrichment <- data.frame(KEGG_enrichment)
                                                    ID_vec <- KEGG_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    KEGG_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DAVID_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(KEGG_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DAVID_KEGG_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  #GO
                                                  tryCatch(expr = {
                                                    BP_GO_enrichment <- enrichDAVID(gene = filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    BP_GO_enrichment <- data.frame(BP_GO_enrichment)
                                                    ID_vec <- BP_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    BP_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DAVID_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(BP_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DAVID_GO_BP_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  tryCatch(expr = {
                                                    CC_GO_enrichment <- enrichDAVID(gene = filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    CC_GO_enrichment <- data.frame(CC_GO_enrichment)
                                                    ID_vec <- CC_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    CC_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DAVID_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(CC_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DAVID_GO_CC_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                  tryCatch(expr = {
                                                    MF_GO_enrichment <- enrichDAVID(gene = filtered_list, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_DIRECT", species = genomeFull, david.user = input$davidEmail2, minGSSize = 0, maxGSSize = NA)
                                                    MF_GO_enrichment <- data.frame(MF_GO_enrichment)
                                                    ID_vec <- MF_GO_enrichment$geneID
                                                    if(length(ID_vec) > 0){
                                                      for(j in 1:length(ID_vec)){
                                                        gene_split <- unlist(strsplit(ID_vec[j], "/"))
                                                        for(k in 1:length(gene_split)){
                                                          gene_split[k] <- rownames(deg_expTab)[which(deg_expTab$entrez_id == gene_split[k])]
                                                        }
                                                        ID_vec[j] <- paste(gene_split, collapse = "/")
                                                      }
                                                    }
                                                    MF_GO_enrichment$geneID <- ID_vec
                                                    fileName <- paste(comparisonName, "_DAVID_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    write.table(MF_GO_enrichment, fileName, sep = "\t", quote = F, row.names = F)
                                                  },
                                                  warning = function(w){
                                                    headers <- paste(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"), collapse = "\t")
                                                    fileName <- paste(comparisonName, "_DAVID_GO_MF_enrichment_DEanalysis.txt", sep = "")
                                                    conn <- file(fileName)
                                                    writeLines(headers, conn)
                                                    close(conn)
                                                  })
                                                }
                                              }
                                            }
                                            
                                            ####################################
                                            # ENRICHR analysis
                                            ####################################
                                            
                                            if(offlineMode == "enrichr"){
                                              statisticToUse <- input$statisticToUse2
                                              if(statisticToUse == "PValue"){
                                                thresholdToUse <- PvalCut
                                              } else {
                                                thresholdToUse <- FDRCut
                                              }
                                              if(orgEnrichment == "hsa"){
                                                databases <- c("KEGG_2019_Human", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
                                              }
                                              if(orgEnrichment == "mmu"){
                                                databases <- c("KEGG_2019_Mouse", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
                                              }
                                              separateDirections <- input$separateDirections2
                                              if(separateDirections == "yes"){
                                                # DOWN REGULATED
                                                down_filtered_list <- rownames(deg_expTab)[which(deg_expTab[,statisticToUse] < thresholdToUse & deg_expTab[,"FC"] <= -FCCut)]
                                                if(length(down_filtered_list) > 1){
                                                  down_enriched_terms <- enrichr(down_filtered_list, databases = databases)
                                                  if(orgEnrichment == "hsa"){
                                                    write.table(down_enriched_terms$KEGG_2019_Human, paste(comparisonName, "_DOWN_enrichR_KEGG_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  }
                                                  if(orgEnrichment == "mmu"){
                                                    write.table(down_enriched_terms$KEGG_2019_Mouse, paste(comparisonName, "_DOWN_enrichR_KEGG_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  }
                                                  write.table(down_enriched_terms$GO_Biological_Process_2018, paste(comparisonName, "_DOWN_enrichR_GO_BP_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  write.table(down_enriched_terms$GO_Cellular_Component_2018, paste(comparisonName, "_DOWN_enrichR_GO_CC_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  write.table(down_enriched_terms$GO_Molecular_Function_2018, paste(comparisonName, "_DOWN_enrichR_GO_MF_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                }
                                                # UP REGULATED
                                                up_filtered_list <- rownames(deg_expTab)[which(deg_expTab[,statisticToUse] < thresholdToUse & deg_expTab[,"FC"] >= FCCut)]
                                                if(length(up_filtered_list) > 1){
                                                  up_enriched_terms <- enrichr(up_filtered_list, databases = databases)
                                                  if(orgEnrichment == "hsa"){
                                                    write.table(up_enriched_terms$KEGG_2019_Human, paste(comparisonName, "_UP_enrichR_KEGG_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  }
                                                  if(orgEnrichment == "mmu"){
                                                    write.table(up_enriched_terms$KEGG_2019_Mouse, paste(comparisonName, "_UP_enrichR_KEGG_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  }
                                                  write.table(up_enriched_terms$GO_Biological_Process_2018, paste(comparisonName, "_UP_enrichR_GO_BP_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  write.table(up_enriched_terms$GO_Cellular_Component_2018, paste(comparisonName, "_UP_enrichR_GO_CC_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  write.table(up_enriched_terms$GO_Molecular_Function_2018, paste(comparisonName, "_UP_enrichR_GO_MF_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                }
                                              } else {
                                                filtered_list <- rownames(deg_expTab)[which(deg_expTab[,statisticToUse] < thresholdToUse & abs(deg_expTab[,"FC"]) >= FCCut)]
                                                if(length(filtered_list) > 1){
                                                  enriched_terms <- enrichr(filtered_list, databases = databases)
                                                  if(orgEnrichment == "hsa"){
                                                    write.table(enriched_terms$KEGG_2019_Human, paste(comparisonName, "_enrichR_KEGG_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  }
                                                  if(orgEnrichment == "mmu"){
                                                    write.table(enriched_terms$KEGG_2019_Mouse, paste(comparisonName, "_enrichR_KEGG_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  }
                                                  write.table(enriched_terms$GO_Biological_Process_2018, paste(comparisonName, "_enrichR_GO_BP_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  write.table(enriched_terms$GO_Cellular_Component_2018, paste(comparisonName, "_enrichR_GO_CC_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                  write.table(enriched_terms$GO_Molecular_Function_2018, paste(comparisonName, "_enrichR_GO_MF_enrichment_DEanalysis.txt", sep = ""), sep = "\t", quote = F, row.names = F)
                                                }
                                              }
                                            }
                                            
                                            # Volcano Plot with points FC>=|2| AND FDR<0.05 highlighted
                                            fileName <- paste(comparisonName, "_volcanoPlot.png", sep = "")
                                            png(fileName)
                                            title <- paste(comparisonName, ": FC vs FDR", sep = "")
                                            with(deg_expTab,plot(deg_expTab$logFC,-log10(deg_expTab$FDR),pch=20,main=title, xlab="log2FoldChange",ylab="-log10(FDR)"))
                                            with(subset(deg_expTab,FDR<FDRCut),points(logFC,-log10(FDR),pch=20,col="red"))
                                            with(subset(deg_expTab,abs(logFC)>log(FCCut)),points(logFC,-log10(FDR),pch=20,col="blue"))
                                            with(subset(deg_expTab,FDR<FDRCut&abs(logFC)>log(FCCut)),points(logFC,-log10(FDR),pch=20,col="green"))
                                            dev.off()
                                            
                                            # Filter biologically significant results FC
                                            deg_FC2 <- deg_expTab[which(abs(deg_expTab$FC) >= FCCut),]
                                            deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
                                            if(nrow(deg_FC2) > 1){
                                              gene_list <- rownames(deg_FC2)
                                              mat <- tpm
                                              mat[mat==0] <- 0.001
                                              gn <- rownames(mat)
                                              index <- which(gn %in% gene_list)
                                              tpm_FC2 <- as.matrix(mat[index,])
                                              # Z-score scaling on log-transformed data
                                              if(zScoreScaling == "row"){
                                                tpm_FC2 <- t(scale(t(log2(tpm_FC2))))
                                              }else{
                                                tpm_FC2 <- log2(tpm_FC2)
                                              }
                                              
                                              # Display row name parameter
                                              if("FConly" %in% rowNameCheck){
                                                displayRowName <- TRUE
                                              }else{
                                                displayRowName <- FALSE
                                              }
                                              
                                              fileName <- paste(comparisonName, "_SigGenes_FC", FCCut, ".png", sep = "")
                                              title <- paste(comparisonName, "\n FC>=|", FCCut, "|", sep = "")
                                              png(fileName, height = 700, width = 800, res = 150)
                                              print(Heatmap(matrix=tpm_FC2, column_title=title, 
                                                            name="scale", top_annotation = ha1,
                                                            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                            show_heatmap_legend = T,  
                                                            clustering_distance_rows = distanceMethod,
                                                            clustering_method_rows = clusterMethod,
                                                            clustering_distance_columns = distanceMethod,
                                                            clustering_method_columns = clusterMethod, 
                                                            column_names_gp = gpar(fontsize = 5),
                                                            row_names_gp = gpar(fontsize = 5),
                                                            show_row_names = displayRowName))
                                              dev.off()
                                            }
                                            # Filter biologically significant results FC and statistically significant results pval
                                            deg_FC2 <- deg_expTab[which((abs(deg_expTab$FC) >= FCCut) & (deg_expTab$FDR < PvalCut)),]
                                            deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
                                            if(nrow(deg_FC2) > 1){
                                              gene_list <- rownames(deg_FC2)
                                              mat <- tpm
                                              mat[mat==0] <- 0.001
                                              gn <- rownames(mat)
                                              index <- which(gn %in% gene_list)
                                              tpm_FC2 <- as.matrix(mat[index,])
                                              # Z-score scaling on log-transformed data
                                              if(zScoreScaling == "row"){
                                                tpm_FC2 <- t(scale(t(log2(tpm_FC2))))
                                              }else{
                                                tpm_FC2 <- log2(tpm_FC2)
                                              }
                                              
                                              # Display row name parameter
                                              if("FCsigPval" %in% rowNameCheck){
                                                displayRowName <- TRUE
                                              }else{
                                                displayRowName <- FALSE
                                              }
                                              
                                              fileName<-paste(comparisonName,"_SigGenes_FC", FCCut, "_sigPval.png", sep = "")
                                              title<-paste(comparisonName,"\n FC>=|", FCCut, "| & Pval<", PvalCut, sep = "")
                                              png(fileName,width=800,height=700,res=150)
                                              print(Heatmap(matrix=tpm_FC2, column_title=title, 
                                                            name="scale", top_annotation = ha1, 
                                                            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                            show_heatmap_legend = T, 
                                                            clustering_distance_rows = distanceMethod,
                                                            clustering_method_rows = clusterMethod,
                                                            clustering_distance_columns = distanceMethod,
                                                            clustering_method_columns = clusterMethod,
                                                            column_names_gp = gpar(fontsize = 5),
                                                            row_names_gp = gpar(fontsize = 5),
                                                            show_row_names = displayRowName))
                                              dev.off()
                                            }
                                            # Filter biologically significant results FC and statistically significant results fdr
                                            deg_FC2 <- deg_expTab[which((abs(deg_expTab$FC) >= FCCut) & (deg_expTab$FDR < FDRCut)),]
                                            deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
                                            if(nrow(deg_FC2) > 1){
                                              gene_list <- rownames(deg_FC2)
                                              mat <- tpm
                                              mat[mat==0] <- 0.001
                                              gn <- rownames(mat)
                                              index <- which(gn %in% gene_list)
                                              tpm_FC2 <- as.matrix(mat[index,])
                                              # Z-score scaling on log-transformed data
                                              if(zScoreScaling == "row"){
                                                tpm_FC2 <- t(scale(t(log2(tpm_FC2))))
                                              }else{
                                                tpm_FC2 <- log2(tpm_FC2)
                                              }
                                              
                                              # Display row name parameter
                                              if("FCsigFDR" %in% rowNameCheck){
                                                displayRowName <- TRUE
                                              }else{
                                                displayRowName <- FALSE
                                              }
                                              
                                              fileName<-paste(comparisonName,"_SigGenes_FC", FCCut, "_sigFDR.png", sep = "")
                                              title<-paste(comparisonName,"\n FC>=|", FCCut, "| & FDR<", FDRCut, sep = "")
                                              png(fileName,width=800,height=700,res=150)
                                              print(Heatmap(matrix=tpm_FC2, column_title=title, 
                                                            name="scale", top_annotation = ha1, 
                                                            column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                                                            show_heatmap_legend = T, 
                                                            clustering_distance_rows = distanceMethod,
                                                            clustering_method_rows = clusterMethod,
                                                            clustering_distance_columns = distanceMethod,
                                                            clustering_method_columns = clusterMethod,
                                                            column_names_gp = gpar(fontsize = 5),
                                                            row_names_gp = gpar(fontsize = 5),
                                                            show_row_names = displayRowName))
                                              dev.off()
                                            }
                                            # Filter biologically significant results FC and statistically significant results pval and average groups
                                            deg_FC2 <- deg_expTab[which((abs(deg_expTab$FC) >= FCCut) & (deg_expTab$FDR < PvalCut)),]
                                            deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
                                            if(nrow(deg_FC2) > 1){
                                              gene_list <- rownames(deg_FC2)
                                              mat <- tpm
                                              mat[mat==0] <- 0.001
                                              gn <- rownames(mat)
                                              index <- which(gn %in% gene_list)
                                              tpm_FC2 <- as.matrix(mat[index,])
                                              for(j in 1:length(unique(group))){
                                                groupToAverage <- unique(group)[j]
                                                groupIdx <- which(group == groupToAverage)
                                                if(j == 1){
                                                  groupData <- rowMeans(tpm_FC2[,groupIdx, drop = FALSE])
                                                  deg_FC2_group_average <- groupData
                                                }
                                                if(j > 1){
                                                  groupData <- rowMeans(tpm_FC2[,groupIdx])
                                                  deg_FC2_group_average <- cbind(deg_FC2_group_average, groupData)
                                                }
                                              }
                                              colnames(deg_FC2_group_average) <- unique(group)
                                              # Z-score scaling on log-transformed data
                                              if(zScoreScaling == "row"){
                                                deg_FC2_group_average <- t(scale(t(log2(deg_FC2_group_average))))
                                              }else{
                                                deg_FC2_group_average <- log2(deg_FC2_group_average)
                                              }
                                              
                                              # Display row name parameter
                                              if("FCsigPvalave" %in% rowNameCheck){
                                                displayRowName <- TRUE
                                              }else{
                                                displayRowName <- FALSE
                                              }
                                              
                                              fileName<-paste(comparisonName, "_SigGenes_FC", FCCut, "_sigPval_averaged_group.png", sep = "")
                                              title<-paste(comparisonName,"\n FC>=|", FCCut, "| & Pval<", PvalCut, " averaged groups", sep = "")
                                              png(fileName,width=800,height=700,res=150)
                                              print(Heatmap(matrix=deg_FC2_group_average, column_title=title, 
                                                            name="scale", top_annotation = ha2, 
                                                            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                            show_heatmap_legend = T, 
                                                            clustering_distance_rows = distanceMethod,
                                                            clustering_method_rows = clusterMethod,
                                                            clustering_distance_columns = distanceMethod,
                                                            clustering_method_columns = clusterMethod,
                                                            column_names_gp = gpar(fontsize = 5),
                                                            row_names_gp = gpar(fontsize = 5),
                                                            show_row_names = displayRowName))
                                              dev.off()
                                            }
                                            # Filter biologically significant results FC and statistically significant results fdr and average groups
                                            deg_FC2 <- deg_expTab[which((abs(deg_expTab$FC) >= FCCut) & (deg_expTab$FDR < FDRCut)),]
                                            deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
                                            if(nrow(deg_FC2) > 1){
                                              gene_list <- rownames(deg_FC2)
                                              mat <- tpm
                                              mat[mat==0] <- 0.001
                                              gn <- rownames(mat)
                                              index <- which(gn %in% gene_list)
                                              tpm_FC2 <- as.matrix(mat[index,])
                                              for(j in 1:length(unique(group))){
                                                groupToAverage <- unique(group)[j]
                                                groupIdx <- which(group == groupToAverage)
                                                if(j == 1){
                                                  groupData <- rowMeans(tpm_FC2[,groupIdx, drop = FALSE])
                                                  deg_FC2_group_average <- groupData
                                                }
                                                if(j > 1){
                                                  groupData <- rowMeans(tpm_FC2[,groupIdx])
                                                  deg_FC2_group_average <- cbind(deg_FC2_group_average, groupData)
                                                }
                                              }
                                              colnames(deg_FC2_group_average) <- unique(group)
                                              # Z-score scaling on log-transformed data
                                              if(zScoreScaling == "row"){
                                                deg_FC2_group_average <- t(scale(t(log2(deg_FC2_group_average))))
                                              }else{
                                                deg_FC2_group_average <- log2(deg_FC2_group_average)
                                              }
                                              
                                              # Display row name parameter
                                              if("FCsigFDRave" %in% rowNameCheck){
                                                displayRowName <- TRUE
                                              }else{
                                                displayRowName <- FALSE
                                              }
                                              
                                              fileName<-paste(comparisonName, "_SigGenes_FC", FCCut, "_sigFDR_averaged_group.png", sep = "")
                                              title<-paste(comparisonName,"\n FC>=|", FCCut, "| & FDR<", FDRCut, " averaged groups", sep = "")
                                              png(fileName,width=800,height=700,res=150)
                                              print(Heatmap(matrix=deg_FC2_group_average, column_title=title, 
                                                            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                            name="scale", top_annotation = ha2, 
                                                            show_heatmap_legend = T, 
                                                            clustering_distance_rows = distanceMethod,
                                                            clustering_method_rows = clusterMethod,
                                                            clustering_distance_columns = distanceMethod,
                                                            clustering_method_columns = clusterMethod,
                                                            column_names_gp = gpar(fontsize = 5),
                                                            row_names_gp = gpar(fontsize = 5),
                                                            show_row_names = displayRowName))
                                              dev.off()
                                            }
                                          }
                                          })
                             incProgress(amount = 1/3, message = NULL)
                           })
             })