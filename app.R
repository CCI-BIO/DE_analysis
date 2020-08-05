# DEanalysis Rshiny app
# Shiny app to streamline DGE script written by Chelsea Mayoh
# 
# Created by Nisitha Jayatilleke
# Date: 19/03/2019
# Last Updated: 13/05/2019 

##### BIOCONDUCTOR INSTALL
Bioconductor_list <- c("edgeR", "ComplexHeatmap", "clusterProfiler", "RDAVIDWebService", "AnnotationDbi", "GO.db", "KEGG.db", "org.Hs.eg.db", "org.Mm.eg.db", "pathview")
for(i in 1:length(Bioconductor_list)){
  if(!requireNamespace(Bioconductor_list[i], quietly = TRUE)){
    BiocManager::install(Bioconductor_list[i], update = FALSE)
  }
}
# Import relevant libraries
library(shiny) # CRAN
library(shinyjs) # CRAN
library(shinyFiles) # CRAN
library(edgeR) # Bioconductor
library(gtools) # CRAN
library(gplots) # CRAN
library(ggplot2) # CRAN
library(ggrepel) # CRAN
library(ComplexHeatmap) # Bioconductor
library(ggfortify) # CRAN
library(org.Mm.eg.db) # Bioconductor
library(org.Hs.eg.db) # Bioconductor
library(knitr) # CRAN
library(clusterProfiler) # Bioconductor
library(RDAVIDWebService) # Bioconductor
library(enrichR) # CRAN
library(httr) # CRAN
library(V8) # CRAN
library(pathview) #Bioconductor

# Increase max file size
options(shiny.maxRequestSize = 500*1024^2)

# Knit R markdown files
rmdfiles <- c("instructions.Rmd", "about.Rmd")
sapply(rmdfiles, knit, quiet = T)

# Button checks
PCArefresh <- FALSE
combinationLastSelected <- NULL
Volcanorefresh <- FALSE
Heatmaprefresh <- FALSE
heatmapGroupCheck <- FALSE
heatmapGroupLastSelected <- NULL
emailCheck <- FALSE
UPDOWNcheck <- FALSE
barPlotGroupLastSelected <- NULL
barPlotrefresh <- FALSE

# Other variables
max.row <- 0

# JS code
jscode <- "shinyjs.toTop = function() {window.scrollTo(0, 0)}"

# Define UI for application
ui <- fluidPage(
  useShinyjs(),
  # Initalise navigation bar
  navbarPage(id = "tabs", title = "RNA-seq DE analysis",
             # Instructions tab
             tabPanel("Instructions",
                      fluidRow(
                        includeMarkdown("instructions.md")
                      ),
                      extendShinyjs(text = jscode),
                      actionButton("toTop", "Jump to top"),
                      br(), br(), br()
             ),
             # First step tab
             tabPanel(value = "step1", "Step 1: Data Input", 
                      # Page title
                      titlePanel("Differential Expression Analysis using EdgeR"),
                      # Input widget for raw expression counts
                      fluidRow(
                        column(4,
                               h3("Input Raw Expression data"),
                               fileInput(inputId = "file1", HTML(paste(h5("File input: Raw counts file. File containing raw expression counts output from
                                                               the RNA-seq pipeline."), '\n', h5("(eg. GeneExpression_rawCounts.txt)"))))),
                        column(4, br(), br(), br(), br(), br(), br(),
                               textOutput("columnMissing1")),
                        tags$head(tags$style("#columnMissing1{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"
                        )
                        )
                      ),
                      # Input widget for TPM expression counts
                      fluidRow(
                        column(4,
                               h3("Input TPM Expression data"),
                               fileInput("file2", HTML(paste(h5("File input: TPM counts file. File containing normalised/TPM counts output from
                                                      the RNA-seq pipeline."), '\n', h5("(eg. GeneExpression_TPM_Counts.txt)"))))),
                        column(4,br(), br(), br(), br(), br(), br(),
                               textOutput("columnMissing2")),
                        tags$head(tags$style("#columnMissing2{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"
                        )
                        )
                      ),
                      # Input widget for group information
                      fluidRow(
                        column(4,
                               h3("Input Group data"),
                               fileInput("file3", HTML(paste(h5("File input: Tab-delimited table with minimum of 2 columns with 'Samples' and 'Group' headers (see instructions)."), 
                                                             '\n', 
                                                             h5("'Samples' column should correspond with sample names in raw/TPM counts file."),
                                                             '\n',
                                                             h5("'Group' column(s) can contain any character strings for given groups a sample belongs to."),
                                                             '\n',
                                                             h5("Multiple group columns can exist with different header names.")))
                               )
                        ),
                        column(4,br(), br(), br(), br(), br(), br(), br(), br(),
                               textOutput("sampleMatch")),
                        tags$head(tags$style("#sampleMatch{color: red;
                                  font-size: 20px;
                                  font-style: bold;
                                  }"
                                  )
                                  )
                      ),
                      # Select widget for tables
                      fluidRow(br(),
                               column(4,
                                      selectInput("tableSelect", 
                                                  h3("View table:"),
                                                  choices = list("Raw expression data" = "file1",
                                                                 "TPM data" = "file2",
                                                                 "Group data"= "file3")
                                      ))
                      ),
                      # Visualise selected table
                      fluidRow(
                        mainPanel(
                          dataTableOutput("table_preview")
                        )
                      )
             ),
             
             tabPanel(value = "step2", "Step 2: PCA Sample Filtering",
                      sidebarLayout(
                        fluid = TRUE, position = "right",
                        sidebarPanel(width = 5,
                                     # Widget to select samples
                                     uiOutput("sampleChoice"),
                                     uiOutput("groupChoice"),
                                     actionButton("refresh1", "Generate/Refresh Plot"), h1(),
                                     textInput("PCAfilename", "PCA filename to save:"),
                                     downloadButton("downloadPCA")
                        ),
                        mainPanel(width = 7,
                          # Page title
                          h2("Choose samples to keep"),
                          h5("Select samples in box and press delete to remove. 
                                         Generate PCA plot by clicking Generate/Refresh Plot. 
                                         To add samples back to the list, click on the white space inside the box to bring up a selection list. 
                                         Simply click on the sample names to add them back to the list."),
                          h5("The PCA plot is generated using log2(TPM counts) and coloured by the specifed groupings."),
                          # Visualise resulting PCA plot
                          h3("PCA plot:"),
                          plotOutput("pcaPlot")
                          )
                        )
             ),
             
             tabPanel(value = "step3", "Step 3: DE Parameters",
                      # Page title
                      titlePanel("Change DE analysis parameters"),
                      fluidRow(
                        mainPanel("Alter Differential Expression analysis parameters based on your own requirements."),
                        br(), br()
                      ),
                      # Widget to select CPM and percentage samples
                      fluidRow(
                        column(4,
                               numericInput("num1", 
                                            h4("Counts per million (cpm) threshold:"), 
                                            min = 1, value = 1)
                        ),
                        column(4,
                               sliderInput("slider1", 
                                           h4("Percentage of samples to meet cpm threshold:"),
                                           min = 0, max = 100, value = 50)
                        )
                      ),
                      # Widget to select FC, Pval and FDR thresholds
                      fluidRow(
                        column(4,
                               numericInput("FCgroup", 
                                            h4("Fold change threshold:"), 
                                            min = 1.1, max = 3, step = 0.1,
                                            value = 2)
                        ),
                        column(4,
                               numericInput("PvalInput", 
                                            h4("P-value threshold:"), 
                                            min = 0, max = 1, step = 0.01,
                                            value = 0.05)
                        ),
                        column(4,
                               numericInput("FDRInput", 
                                            h4("FDR threshold:"),
                                            min = 0, max = 1, step = 0.01,
                                            value = 0.05)
                        )
                      ),
                      # Widget to select organism
                      fluidRow(
                        column(4,
                               radioButtons("Organism",
                                            h4("Organism:"),
                                            choices = list("Homo sapiens (Human)" = "human",
                                                           "Mus musculus (Mouse)" = "mouse"))
                        )
                      ),
                      # Title separator for GO/KEGG enrichment section
                      fluidRow(
                        h3("GO/KEGG enrichment parameters")
                      ),
                      # Widget to choose parameters for GO/KEGG enrichment
                      fluidRow(
                        column(4,
                               radioButtons("offlineMode",
                                            h4("KEGG/GO enrichment method:"),
                                            choices = list("Enrichr (recommended)" = "enrichr",
                                                           "DAVID web service" = "david",
                                                           "ClusterProfiler" = "clusterprofiler",
                                                           "None" = "none"),
                                            selected = "none")
                        ),
                        column(4,
                               uiOutput("davidEmail"),
                               uiOutput("internalKEGG")),
                        column(4, br(), br(),
                               textOutput("emailWarn1"),
                               textOutput("emailWarn2"),
                               tags$head(tags$style("#emailWarn1{color: green;font-size: 17px;font-style: bold}")),
                               tags$head(tags$style("#emailWarn2{color: red;font-size: 17px;font-style: bold}"))
                        )
                      ),
                      fluidRow(
                        column(4,
                               uiOutput("statisticToUse")
                        ),
                        column(4,
                               uiOutput("separateDirections")
                        )
                      ),
                      # Widget to choose groupings to compare
                      fluidRow(
                        column(4,
                               uiOutput("groupChoiceFinal")
                        )
                      ),
                      # Widget to select comparisons to make
                      fluidRow(
                        column(4,
                               uiOutput("combChoice1")
                        )
                      ),
                      # Title separator for heatmap section
                      fluidRow(
                        h3("Heatmap parameters")
                      ),
                      # Widget to select Heatmap clustering methods 
                      fluidRow(
                        column(4,
                               selectInput("distanceMetric",
                                           h4("Distance metric:"),
                                           choices = list("Euclidean" = "euclidean",
                                                          "Maximum" = "maximum",
                                                          "Manhattan" = "manhattan",
                                                          "Canberra" = "canberra",
                                                          "Binary" = "binary",
                                                          "Minkowski" = "minkowski"),
                                           selected = "euclidean")
                        ),
                        column(4,
                               selectInput("agglomerationMethod",
                                           h4("Agglomeration method:"),
                                           choices = list("Ward.D" = "ward.D",
                                                          "Ward.D2" = "ward.D2",
                                                          "Single linkage" = "single",
                                                          "Complete linkage" = "complete",
                                                          "Average linkage" = "average",
                                                          "Mcquitty" = "mcquitty",
                                                          "Median linkage" = "median",
                                                          "Centroid" = "centroid"),
                                           selected = "average")
                        )
                      ),
                      fluidRow(
                        column(4,
                               uiOutput("colourBar")
                        ),
                        column(4,
                               checkboxGroupInput("plotRowNames",
                                                  h4("Display row names on output plot:"),
                                                  choices = list("Fold-change" = "FConly",
                                                                 "Fold-change & significant FDR" = "FCsigFDR",
                                                                 "Fold-change & significant FDR (averaged)" = "FCsigFDRave",
                                                                 "Fold-change & significant P-value" = "FCsigPval",
                                                                 "Fold-change & significant P-value (averaged)" = "FCsigPvalave"), 
                                                  selected = c("FConly", "FCsigFDR", "FCsigFDRave", "FCsigPval", "FCsigPvalave"))
                        ),
                        column(4,
                               radioButtons("zScoreScaling", 
                                            h4("Z-score scaling:"),
                                            choices = list("Yes" = "row",
                                                           "No" = "none"), 
                                            selected = "row")
                        )
                      ),
                      # Title separator for heatmap section
                      fluidRow(
                        h3("Save location")
                      ),
                      # Widget to select output directory
                      fluidRow(
                        column(4,
                               h4("Output Directory:"),
                               shinyDirButton(id = "OutputDir",
                                              label = "Choose directory", 
                                              title = "Select")
                        )
                      ),
                      # Print selected output directory
                      fluidRow(br(),
                               column(8,
                                      verbatimTextOutput("dirLoc", placeholder = TRUE))
                      ),
                      # Button to submit and run pipeline
                      fluidRow(br(),
                               column(4,
                                      actionButton(inputId = "submit1",
                                                   label = "Submit")),
                               br(),br(),br()
                      )
             ),
             
             tabPanel(value = "step4", "Step 4: Generate Custom Plots",
                      # Section header and description
                      fluidRow(column(12, h3("Data location"))),
                      fluidRow(column(12, h5("Location containing Differential Expression output files. The location will default to the selected output directory in step 3.
                                  Plot generation requires all output files from this program to be present in the directory."))),
                      # Widget to select output file directory
                      fluidRow(
                        column(4,
                               h4("Output file directory:"),
                               shinyDirButton(id = "OutputDirNew",
                                              label = "Choose directory",
                                              title = "Select")
                        )
                      ),
                      # Print selected output directory
                      fluidRow(br(),
                               column(8,
                                      verbatimTextOutput("dirLocNew", placeholder = TRUE))
                      ),
                      # List all potential files
                      fluidRow(br(),
                               column(8,
                                      uiOutput("fileChoice")
                              )
                      ),
                      fluidRow(br(),
                               tabsetPanel(
                                 tabPanel(title = "Volcano Plot", value = "tab1", 
                                          sidebarLayout(
                                            fluid = TRUE, position = "right",
                                            sidebarPanel(
                                              radioButtons("statSelect",
                                                           h4("Select threshold statistic to use:"), 
                                                           choices = list("P-value" = "pval",
                                                                          "FDR" = "fdr"),
                                                           selected = "fdr"),
                                              numericInput("volcanoFC", 
                                                           h4("Fold change threshold:"), 
                                                           min = 1.1, max = 3, step = 0.1,
                                                           value = 2),
                                              uiOutput("selectedStatmethod"),
                                              radioButtons("geneListOption",
                                                           h4("Gene list upload method:"),
                                                           choices = list("Top genes (specified amount)" = "topGenes",
                                                                          "Select input" = "selectInput",
                                                                          "Upload list" = "listInput"),
                                                           selected = "topGenes"),
                                              uiOutput("selectedGeneListmethod"),
                                              fluidRow(column(6, numericInput("volcanoPlotWidth", h4("Plot width (px):"), value = 700, min = 0, max = 2000)),
                                                       column(6, numericInput("volcanoPlotHeight", h4("Plot height (px):"), value = 700, min = 0, max = 2000))),
                                              actionButton("refresh2", "Generate/Refresh Plot"), h3(),
                                              textInput("volcanoFilename", h4("Volcano plot filename to save:")),
                                              downloadButton("downloadVolcano")
                                            ),
                                            mainPanel(
                                              plotOutput("volcanoPlot")
                                            )
                                          )
                                 ),
                                 tabPanel(title = "Heatmap Plot", value = "tab2",
                                        sidebarLayout(
                                          fluid = TRUE, position = "right",
                                          sidebarPanel(
                                            selectInput("distanceMetricCustom",
                                                        h4("Distance metric:"),
                                                        choices = list("Euclidean" = "euclidean",
                                                                       "Maximum" = "maximum",
                                                                       "Manhattan" = "manhattan",
                                                                       "Canberra" = "canberra",
                                                                       "Binary" = "binary",
                                                                       "Minkowski" = "minkowski"),
                                                        selected = "euclidean"),
                                            selectInput("agglomerationMethodCustom",
                                                        h4("Agglomeration method:"),
                                                        choices = list("Ward.D" = "ward.D",
                                                                       "Ward.D2" = "ward.D2",
                                                                       "Single linkage" = "single",
                                                                       "Complete linkage" = "complete",
                                                                       "Average linkage" = "average",
                                                                       "Mcquitty" = "mcquitty",
                                                                       "Median linkage" = "median",
                                                                       "Centroid" = "centroid"),
                                                        selected = "average"),
                                            radioButtons("heatmapGeneSamplesOption",
                                                         h4("Select sample criteria:"), 
                                                         choices = list("All samples" = "all",
                                                                        "Selected samples" = "select"),
                                                         selected = "all"),
                                            uiOutput("heatmapSampleList"),
                                            radioButtons("heatmapGeneListOption",
                                                         h4("Gene list upload method:"),
                                                         choices = list("Top genes (specified amount)" = "topGenes",
                                                                        "Select input" = "selectInput",
                                                                        "Upload list" = "listInput"),
                                                         selected = "topGenes"),
                                            uiOutput("heatmapGeneList"),
                                            fileInput("file6", h4("File input: Group data")),
                                            uiOutput("sampleMatch2"), h3(),
                                            uiOutput("heatmapGroupOption"),
                                            uiOutput("groupSelection"),
                                            uiOutput("heatmapAverageOption"),
                                            uiOutput("groupSelectionAverage"),
                                            radioButtons("heatmapScaling", h4("Z-score scaling:"),
                                                         choices = list("Yes" = "row",
                                                                        "No" = "none"),
                                                         selected = "row"),
                                            radioButtons("heatmapRowNames", h4("Display row names:"),
                                                         choices = list("Yes" = "yes",
                                                                        "No" = "no"),
                                                         selected = "yes"),
                                            radioButtons("heatmapColumnNames", h4("Display column names:"),
                                                         choices = list("Yes" = "yes",
                                                                        "No" = "no"),
                                                         selected = "yes"),
                                            radioButtons("heatmapTranspose", h4("Transpose heatmap"),
                                                         choices = list("Yes" = "yes",
                                                                        "No" = "no"),
                                                         selected = "no"),
                                            fluidRow(column(6, numericInput("heatmapPlotWidth", h4("Plot width (px):"), value = 700, min = 0, max = 2000)),
                                                     column(6, numericInput("heatmapPlotHeight", h4("Plot height (px):"), value = 700, min = 0, max = 2000))),
                                            actionButton("refresh3", "Generate/Refresh Plot"), h3(),
                                            textInput("heatmapFilename", h4("Heatmap plot filename to save:")),
                                            downloadButton("downloadHeatmap")
                                          ),
                                          mainPanel(
                                            plotOutput("heatmapPlot")
                                          )
                                        )         
                                 ),
                                 tabPanel(title = "Enrichment Bar Plot", value = "tab3",
                                          sidebarLayout(
                                            fluid = TRUE, position = "right",
                                            sidebarPanel(
                                              uiOutput("chooseDOWNUP"),
                                              checkboxGroupInput("combineGOKEGG", h4("Enrichment terms to plot:"),
                                                                 choices = list("KEGG pathways" = "kegg",
                                                                                "GO Biological Process" = "gobp",
                                                                                "GO Cellular Component" = "gocc",
                                                                                "GO Molecular Function" = "gomf"), 
                                                                 selected = c("kegg", "gobp", "gocc", "gomf")),
                                              radioButtons("plotMethod", h4("Method to select terms to plot:"),
                                                           choices = list("Filter by statistical threshold" = "filter",
                                                                          "Specify top or specific terms" = "select")),
                                              uiOutput("barPlotTermListOption"),
                                              uiOutput("barPlotTermList"),
                                              uiOutput("termChoice"),
                                              uiOutput("termThreshold"),
                                              uiOutput("rownameChoice"),
                                              radioButtons("addOverlapCount", h4("Add gene overlap fractions:"),
                                                           choices = list("Yes" = T,
                                                                           "No" = F),
                                                           selected = F),
                                              radioButtons("reverseOrder", h4("Reverse ordering of terms:"),
                                                           choices = list("Yes" = T,
                                                                          "No" = F),
                                                           selected = F),
                                              fluidRow(column(6, numericInput("barPlotWidth", h4("Plot width (px):"), value = 700, min = 0, max = 2000)),
                                                       column(6, numericInput("barPlotHeight", h4("Plot height (px):"), value = 700, min = 0, max = 2000))),
                                              actionButton("refresh4", "Generate/Refresh Plot"), h3(),
                                              textInput("barPlotFilename", h4("Bar plot filename to save:")),
                                              downloadButton("downloadbarPlot")
                                            ),
                                            mainPanel(
                                              plotOutput("barPlot")
                                            )
                                          )
                                 ),
                                 tabPanel(title = "KEGG Pathview", value = "tab3",
                                          sidebarLayout(
                                            fluid = TRUE, position = "right",
                                            sidebarPanel(
                                              radioButtons("pathviewOrganism", h4("Select organism:"),
                                                           choices = list("Human" = "hsa",
                                                                          "Mouse" = "mmu")),
                                              uiOutput("selectPathviewPath"),
                                              radioButtons("pathviewFilterGeneList", h4("Method for filtering DE gene list:"),
                                                           choices = list("Filter by P-value" = "pvalue",
                                                                          "Filter by FDR" = "fdr")),
                                              uiOutput("pathviewThreshold"),
                                              textInput("pathviewSuffix", label = h4("Suffix for filename:"), value = ""),
                                              actionButton("refresh5", "Generate/Refresh Plot")
                                            ),
                                            mainPanel(
                                              tags$div(
                                                class = "plotBox",
                                                style = "overflow-y:visible; overflow-x:scroll; height:800px",
                                                imageOutput("pathviewPlot")
                                              )
                                            )
                                          )
                                 )
                               )
                      )
             ),
             
             # More section in navigation bar
             navbarMenu("More",
                        # About information
                        tabPanel("About",
                                 fluidRow(
                                   includeMarkdown("about.md")
                                 )
                        )
             )
  )
)


# Server logic for DEanalysis 
server <- function(input, output, session) {
  
  observeEvent(input$toTop, {
    js$toTop();
  })
  
  # Import file input scripts (step 1)
  source("scripts/file_input_DEapp.R", local = T)
  
  # Import PCA plot scripts (step 2)
  source("scripts/PCA_DEapp.R", local = T)
  
  # Import main DE analysis function (step 3)
  source("scripts/main_DEapp.R", local = T)
  
  # Select file location (and update where neccessary) for step 4
  outputLocation <- reactiveValues(output = NULL)
  observeEvent(input$OutputDir,
               {
                 output$dirLocNew <- renderText({parseDirPath(roots = c("R:/" = "R:/", "C:/" = "C:/"), input$OutputDir)})
                 outputLocation$output <- parseDirPath(roots = c("R:/" = "R:/", "C:/" = "C:/"), input$OutputDir)
               }
              )

  shinyDirChoose(input = input, id = "OutputDirNew", roots = c("R:/" = "R:/", "C:/" = "C:/"))
  observeEvent(input$OutputDirNew,
               {
                 output$dirLocNew <- renderText({parseDirPath(roots = c("R:/" = "R:/", "C:/" = "C:/"), input$OutputDirNew)})
                 outputLocation$output <- parseDirPath(roots = c("R:/" = "R:/", "C:/" = "C:/"), input$OutputDirNew)
               }
  )

  shinyjs::disable("geneListOption")
  
  # List files 
  outFileList <- eventReactive(c(input$OutputDir, input$OutputDirNew, input$tabs, outputLocation$output),
                               {
                                 if(!is.null(outputLocation$output) & length(outputLocation$output) >= 1){
                                   outDir <- outputLocation$output
                                   fileList <- list.files(path = outDir, pattern = ".txt", full.names = T)
                                   fileList <- fileList[-grep(pattern = "_tpm_DEanalysis.txt", x = fileList)]
                                   fileList <- fileList[-grep(pattern = "_enrichment_DEanalysis.txt", x = fileList)]
                                   fileList <- fileList[-grep(pattern = "_GSEA_expression_input.txt", x = fileList)]
                                   fileList <- fileList[grep(pattern = "_DEanalysis.txt", x = fileList)]
                                   fileComparisonNames <- gsub(pattern = "_DEanalysis.txt", replacement = "", x = fileList)
                                   fileComparisonNames <- gsub(pattern = paste(outDir, "/", sep = ""), replacement = "", x = fileComparisonNames)
                                   names(fileList) <- fileComparisonNames
                                   tryCatch({volcano.data <- read.delim(fileList[1], header = T, stringsAsFactors = F, sep = "\t")
                                   max.row <<- nrow(volcano.data)
                                   fileList <- as.list(fileList)
                                   enable("geneListOption")
                                   return(fileList)},
                                   error=function(cond){
                                     updateRadioButtons(session, "geneListOption", selected = "topGenes")
                                     disable("geneListOption")
                                     return(NULL)})
                                   
                                 }
                               }
  )
  
  output$fileChoice <- renderUI({
    selectInput("fileChoice2", h4("Comparison to plot:"), choices = NULL)
  })
  
  observe({
    if(!is.null(outFileList())){
      updateSelectInput(session, "fileChoice2", choices = outFileList(), selected = outFileList()[[1]])
    } else {
      updateSelectInput(session, "fileChoice2", choices = list())
    }
  })
  
  # Import heatmap plot scripts (step 4.1)
  source("scripts/volcano_DEapp.R", local = T)
  
  # Import heatmap plot scripts (step 4.2)
  source("scripts/heatmap_DEapp.R", local = T)
   
  # Import bar plot scripts (step 4.3)
  source("scripts/barplot_DEapp.R", local = T)
  
  # Import pathview scripts (step 4.4)
  source("scripts/pathview_DEapp.R", local = T)
  
  # Stop application command required for Inno deployment
  if(!interactive()){
    session$onSessionEnded(function(){
      stopApp()
      q("no")
    })
  }
}

# Run the application 
shinyApp(ui = ui, server = server)

