library(shiny)
library(dplyr)
library(readr)
library(tools)
library(clustifyr)
library(rsconnect)
library(ExperimentHub)
library(Seurat)
library(shinydashboard)


options(shiny.maxRequestSize = 1500*1024^2)
options(repos = BiocManager::repositories())
options(shiny.reactlog = TRUE)

eh <- ExperimentHub()
refs <- query(eh, "clustifyrdatahub")

# Define UI for data upload app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("Clustifyr RShiny App"),
    # #dashboardPage(
    #     dashboardHeader(),
    #     dashboardSidebar(),
    #     dashboardBody()
    # ),
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Select a file ----
            fileInput("file1", "Choose Matrix File",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv",
                                 '.xlsx',
                                 ".tsv",
                                 ".rds",
                                 ".rda")
            ),
            actionButton("matrixPopup", "Display UMI Matrix in popup window"),
            fileInput("file2", "Choose Metadata File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv",
                                 '.xlsx',
                                 ".tsv",
                                 ".rds",
                                 ".rda")
                      ),
            
            actionButton("metadataPopup", "Display Metadata table in popup"),
            tags$hr(),
            
            # Input: Checkbox if file has header ----
            checkboxInput("header", "Header", TRUE),
            
            
            # Horizontal line ----
            tags$hr(),
            
            # Input: Select separator ----
            radioButtons("sep", "Separator",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t"),
                         selected = ","),
            
            # Input: Select number of rows to display ----
            radioButtons("disp", "Display",
                         choices = c(Head = "head",
                                     All = "all"),
                         selected = "head"),
            
            selectInput("metadataCellType", "Cell Type Metadata Column:",
                        choice=list("")
                        ),
            
            
            helpText("Choose cell type metadata column for average_clusters function"),
            hr(),
            selectInput("dataHubReference", "ClustifyrDataHub Reference:", 
                        choices=list("ref_MCA","ref_tabula_muris_drop","ref_tabula_muris_facs",
                        "ref_mouse.rnaseq", "ref_moca_main", "ref_immgen", "ref_hema_microarray",
                        "ref_cortex_dev", "ref_pan_indrop", "ref_pan_smartseq2", 
                        "ref_mouse_atlas")),
            helpText("Choose reference cell atlas for clustify function"),
            hr(),
            helpText("Choose cell reference for clustify function"),
            downloadButton("downloadReference", "Download reference matrix"),
            downloadButton("downloadClustify", "Download clustify matrix")
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            #tableOutput("contents1"), #UMI Count Matrix
            #tags$hr(),
            #tableOutput("contents2"), #Metadata table
            #tags$hr(),
            tableOutput("reference"), #Reference Matrix
            tags$hr(),
            tableOutput("clustify") #Clustify Matrix
        )
        
    )
)

# Define server logic to read selected file ----
server <- function(input, output, session) {
    
    data1Display <- reactive({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        file <- input$file1
        fileTypeFile1 <- tools::file_ext(file$datapath)
        req(file)
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        if (fileTypeFile1 == "csv")
        {
            # df1 <- read_csv(file$datapath,
            #                   header = input$header,
            #                   sep = input$sep)
            df1 <- read.csv(file$datapath,
                            header = input$header,
                            sep = input$sep)
            rownames(df1) <- df1[, 1]
            #df1[, 1] <- NULL
        }
        else if (fileTypeFile1 == "tsv")
        {
            df1 <- read_tsv(file$datapath,
                            header = input$header)
            rownames(df1) <- df1[, 1]
            #df1[, 1] <- NULL
        }
        else
        {
            df1 <- load(file$datapath)
        }
        df1
    })
    
    data1 <- reactive({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        file <- input$file1
        fileTypeFile1 <- tools::file_ext(file$datapath)
        req(file)
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        if (fileTypeFile1 == "csv")
        {
            # df1 <- read_csv(file$datapath,
            #                   header = input$header,
            #                   sep = input$sep)
            df1 <- read.csv(file$datapath,
                            header = input$header,
                            sep = input$sep)
            rownames(df1) <- df1[, 1]
            df1[, 1] <- NULL
        }
        else if (fileTypeFile1 == "tsv")
        {
            df1 <- read_tsv(file$datapath,
                            header = input$header)
            rownames(df1) <- df1[, 1]
            df1[, 1] <- NULL
        }
        else
        {
            df1 <- load(file$datapath)
        }
        df1
    })
    
    data2 <- reactive({
        file <- input$file2    
        fileTypeFile2 <- tools::file_ext(file$datapath)
        req(file)
        if (fileTypeFile2 == "csv")
        {
            df2 <- read.csv(file$datapath,
                            header = input$header,
                            sep = input$sep)
            # df2 <- read_csv(file$datapath,
            #                   header = input$header,
            #                   sep = input$sep)
            
        }
        else if (fileTypeFile2 == "tsv")
        {
            df2 <- read_tsv(file$datapath,
                            header = input$header) 
        }
        else
        {
            df2 <- load(file$datapath)
        }
        updateSelectInput(session, "metadataCellType",
                          choices = c(NA, colnames(df2)),
                          selected = NA
        )
        df2
    })
   
    output$contents1 <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        # file <- input$file1
        # fileTypeFile1 <- tools::file_ext(file$datapath)
        # req(file)
        # # when reading semicolon separated files,
        # # having a comma separator causes `read.csv` to error
        #     if (fileTypeFile1 == "csv")
        #     {
        #         df1 <- read_csv(file$datapath,
        #                         header = input$header,
        #                         sep = input$sep,
        #                         quote = input$quote)
        #         rownames(df1) <- df1[, 1]
        #         df1[, 1] <- NULL
        #     }
        #     else if (fileTypeFile1 == "tsv")
        #     {
        #         df1 <- read_tsv(file$datapath,
        #                         header = input$header,
        #                         quote = input$quote)
        #         rownames(df1) <- df1[, 1]
        #         df1[, 1] <- NULL
        #     }
        #     else
        #     {
        #         df1 <- load(file$datapath)
        #     }
            df1 <- data1Display()   
            #file 1
            if(input$disp == "head") {
                return(head(df1, cols = 5))
            }
            else {
                return(df1)
            }
        })
    
    output$contents2 <- renderTable({
        # file <- input$file2    
        # fileTypeFile2 <- tools::file_ext(file$datapath)
        # req(file)
        # if (fileTypeFile2 == "csv")
        # {
        #     df2 <- read_delim(file$datapath, 
        #                      ",", 
        #                     header = input$header,
        #                     sep = input$sep,
        #                     quote = input$quote)
        # }
        # else if (fileTypeFile2 == "tsv")
        # {
        #     df2 <- read_tsv(file$datapath,
        #                     header = input$header,
        #                     quote = input$quote)   
        # }
        # else
        # {
        #     df2 <- load(file$datapath)
        # }
        df2 <- data2()
        #file 2
        if(input$disp == "head") {
            return(head(df2))
        }
        else {
            return(df2)
        }
    })
    
    observeEvent(input$matrixPopup, {
        showModal(modalDialog(
            tags$caption("UMI Count Matrix"),
            DT::renderDataTable({
                matrixRender <- head(data1Display())
                DT::datatable(matrixRender, escape = FALSE)
            }),
            easyClose = TRUE
        ))
    })
    
    observeEvent(input$metadataPopup, {
            showModal(modalDialog(
                tags$caption("Metadata table"),
                DT::renderDataTable({
                    matrixRender <- head(data2())
                    DT::datatable(matrixRender, escape = FALSE)
                }),
                easyClose = TRUE
            ))
        })
    
    # dataRef <- reactive({
    #     reference_matrix <- average_clusters(mat = data1(), metadata = data2()[[input$metadataCellType]], if_log = FALSE)
    #     reference_matrix
    # })
    
    #dataClustify <- reactive({
    #    eh <- ExperimentHub()
    # query
    #    refs <- query(eh, "clustifyrdatahub")
    #refs <- listResources(eh, "clustifyrdatahub")
    #    benchmarkRef <- loadResources(eh, "clustifyrdatahub", input$dataHubReference)[[1]]
    
    #    UMIMatrix <- data1()
    #    matrixSeuratObject <- CreateSeuratObject(counts = UMIMatrix, project = "Seurat object matrix", min.cells = 0, min.features = 0) 
    #   matrixSeuratObject <- NormalizeData(matrixSeuratObject)
    #   matrixSeuratObject <- FindVariableFeatures(matrixSeuratObject, selection.method = "vst", nfeatures = 2000)
    
    #  metadataCol <- data2()[[input$metadataCellType]]
    #   # use for classification of cell types
    # res <- clustify(
    #   input = matrixSeuratObject@assays$RNA@data, 
    #    metadata = metadataCol,
    #  ref_mat = benchmarkRef,
    # query_genes = VariableFeatures(matrixSeuratObject)
    #)
    #res
    #})
    
  dataRef <- reactive({
        reference_matrix <- average_clusters(mat = data1(), metadata = data2()[[input$metadataCellType]], if_log = FALSE)
        reference_matrix

  })
  
  dataClustify <- reactive({
        benchmarkRef <- loadResources(eh, "clustifyrdatahub", input$dataHubReference)[[1]]
        
        UMIMatrix <- data1()
        matrixSeuratObject <- CreateSeuratObject(counts = UMIMatrix, project = "Seurat object matrix", min.cells = 0, min.features = 0)
        matrixSeuratObject <- FindVariableFeatures(matrixSeuratObject, selection.method = "vst", nfeatures = 2000)
        
        metadataCol <- data2()[[input$metadataCellType]]
        # use for classification of cell types
        res <- clustify(
            input = matrixSeuratObject@assays$RNA@data, 
            metadata = metadataCol,
            ref_mat = benchmarkRef,
            query_genes = VariableFeatures(matrixSeuratObject)
        )
        res
    })
    

    output$reference <- renderTable({
        reference_matrix <- dataRef()
        return(head(reference_matrix))
    })
    
    output$clustify <- renderTable({
        res <- dataClustify()
        return(head(res))
    })
    #Make plots such as heat maps to compare benchmarking with clustify with actual cell types
    
    referenceDownload <- reactive({
        referenceMatrix <- dataRef()
    })
    
    clustifyDownload <- reactive({
        clustifyMatrix <- dataClustify()
    })
    
    output$downloadReference <- downloadHandler(
        filename = function() {
            paste("reference-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(referenceDownload(), file)
            #mat %>% as_tibble(rownames = "rowname") %>% write_csv("mat.csv")
        }
    )
    output$downloadClustify <- downloadHandler(
        filename = function() {
            paste("clustify-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(clustifyDownload(), file)
            #mat %>% as_tibble(rownames = "rowname") %>% write_csv("mat.csv")
        }
    )
}

# Create Shiny app ----
shinyApp(ui, server)