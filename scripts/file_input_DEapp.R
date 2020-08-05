## file_input_DEapp.R
## Script containing functions to allow file input.
##
## Created by Nisitha Jayatilleke
## Date: 19/08/2019
## Last updated: 19/08/2019

shinyDirChoose(input = input, id = "OutputDir", roots = c("R:/" = "R:/", "C:/" = "C:/"))
observeEvent(input$OutputDir, {output$dirLoc <- renderText({parseDirPath(roots = c("R:/" = "R:/", "C:/" = "C:/"), input$OutputDir)})})

output$table_preview <- renderDataTable({
  inFile <- input[[input$tableSelect]]
  if(is.null(inFile)){
    return(NULL)
  }
  read.delim(inFile$datapath, header = TRUE)
})

observeEvent(input$file1,
             {
               inFile <- input[["file1"]]
               if(is.null(inFile)){
                 return(NULL)
               }
               test <- read.delim(inFile$datapath, header = TRUE)
               if(!("gene_id" %in% colnames(test))){
                 output$columnMissing1 <- renderText({"Missing gene_id column"})
               }
               if("gene_id" %in% colnames(test)){
                 output$columnMissing1 <- renderText({""})
               }
             })

observeEvent(input$file2,
             {
               inFile <- input[["file2"]]
               if(is.null(inFile)){
                 return(NULL)
               }
               test <- read.delim(inFile$datapath, header = TRUE)
               if(!("gene_id" %in% colnames(test))){
                 output$columnMissing2 <- renderText({"Missing gene_id column"})
               }
               if("gene_id" %in% colnames(test)){
                 output$columnMissing2 <- renderText({""})
               }
             })

observeEvent({input$file1
  input$file2
  input$file3},
  {
    inFile1 <- input[["file1"]]
    inFile2 <- input[["file2"]]
    if(is.null(inFile1) & is.null(inFile2)){
      return(NULL)
    }
    inFile3 <- input[["file3"]]
    if(is.null(inFile3)){
      return(NULL)
    }
    if(!is.null(inFile1)){
      test <- read.delim(inFile1$datapath, header = TRUE)
      if("length" %in% colnames(test)){
        test <- test[,-which(colnames(test) %in% "length")]
      }
      if("transcript_id.s." %in% colnames(test)){
        test <- test[,-which(colnames(test) %in% "transcript_id.s.")]
      }
      if("gene_id" %in% colnames(test)){
        test <- test[,-which(colnames(test) %in% "gene_id")]
      }
      test2 <- read.delim(inFile3$datapath, header = TRUE)
      if(nrow(test2) != ncol(test)){
        output$sampleMatch <- renderText({"Number of rows do not match number of samples in count data"})
      }
      if(nrow(test2) == ncol(test)){
        output$sampleMatch <- renderText({""})
      }
    }
    if(!is.null(inFile2)){
      test <- read.delim(inFile2$datapath, header = TRUE)
      if("length" %in% colnames(test)){
        test <- test[,-which(colnames(test) %in% "length")]
      }
      if("transcript_id.s." %in% colnames(test)){
        test <- test[,-which(colnames(test) %in% "transcript_id.s.")]
      }
      if("gene_id" %in% colnames(test)){
        test <- test[,-which(colnames(test) %in% "gene_id")]
      }
      test2 <- read.delim(inFile3$datapath, header = TRUE)
      if(nrow(test2) != ncol(test)){
        output$sampleMatch <- renderText({"Number of rows do not match number of samples in count data"})
      }
      if(nrow(test2) == ncol(test)){
        output$sampleMatch <- renderText({""})
      }
    }
  })