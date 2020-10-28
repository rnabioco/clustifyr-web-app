library(shiny)
library(dplyr)
library(readr)
library(tools)
library(clustifyr)
library(rsconnect)
#library(ExperimentHub)
options(shiny.maxRequestSize = 1500*1024^2)
options(repos = BiocManager::repositories())
options(shiny.reactlog = TRUE)

# Define UI for data upload app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("Clustifyr Reference Matrix Generation"),
    
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
                                 ".rda"),
            ),
            fileInput("file2", "Choose Metadata File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv",
                                 '.xlsx',
                                 ".tsv",
                                 ".rds",
                                 ".rda"),
                      ),
            
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
            
            textInput("metadataCellType", label = h3("Metadata Cell Type"), placeholder = "Enter metadata cell type column"),
            
            selectInput("dataHubReference", "ClustifyrDataHub Reference:", 
                        choices=list("ref_MCA","ref_tabula_muris_drop","ref_tabula_muris_facs",
                        "ref_mouse.rnaseq", "ref_moca_main", "ref_immgen", "ref_hema_microarray",
                        "ref_cortex_dev", "ref_pan_indrop", "ref_pan_smartseq2", 
                        "ref_mouse_atlas")),
            hr(),
            helpText("Choose cell reference for clustify function")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            tableOutput("contents1"),
            tags$hr(),
            tableOutput("contents2"),
            tags$hr(),
            tableOutput("reference"),
            tableOutput("clustify")
        )
        
    )
)

# Define server logic to read selected file ----
server <- function(input, output) {
    
    output$contents1 <- renderTable({
        
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
                df1 <- read.csv(file$datapath,
                                header = input$header,
                                sep = input$sep,
                                quote = input$quote)
                rownames(df1) <- df1[, 1]
                df1[, 1] <- NULL
            }
            else if (fileTypeFile1 == "tsv")
            {
                df1 <- read_tsv(file$datapath,
                                header = input$header,
                                quote = input$quote)
                rownames(df1) <- df1[, 1]
                df1[, 1] <- NULL
            }
            else
            {
                df1 <- load(file$datapath)
            }
                
            #file 1
            if(input$disp == "head") {
                return(head(df1, col = 5))
            }
            else {
                return(df1)
            }
        }
    )
    
    output$contents2 <- renderTable({
        file <- input$file2    
        fileTypeFile2 <- tools::file_ext(file$datapath)
        req(file)
        if (fileTypeFile2 == "csv")
        {
            df2 <- read.csv(file$datapath,
                            header = input$header,
                            sep = input$sep,
                            quote = input$quote)
        }
        else if (fileTypeFile2 == "tsv")
        {
            df2 <- read_tsv(file$datapath,
                            header = input$header,
                            quote = input$quote)   
        }
        else
        {
            df2 <- load(file$datapath)
        }
        
        #file 2
        if(input$disp == "head") {
            return(head(df2))
        }
        else {
            return(df2)
        }
    }
)
    
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
            df1 <- read.csv(file$datapath,
                            header = input$header,
                            sep = input$sep,
                            quote = input$quote)
            rownames(df1) <- df1[, 1]
            df1[, 1] <- NULL
        }
        else if (fileTypeFile1 == "tsv")
        {
            df1 <- read_tsv(file$datapath,
                            header = input$header,
                            quote = input$quote)
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
                            sep = input$sep,
                            quote = input$quote)
            
        }
        else if (fileTypeFile2 == "tsv")
        {
            df2 <- read_tsv(file$datapath,
                            header = input$header,
                            quote = input$quote) 
        }
        else
        {
            df2 <- load(file$datapath)
        }
        df2
    })
        
    output$reference <- renderTable({
        #req(input$file1)
        #req(input$file2)
        #fileTypeFile1 <- file_ext(input$file1)
        #fileTypeFile2 <- file_ext(input$file2)
        #df1 <- read_csv(input$file1$datapath,
         #               header = input$header,
          #              sep = input$sep,
           #             quote = input$quote)
        #df2 <- read_csv(input$file2$datapath,
          #              header = input$header,
         #               sep = input$sep,
           #             quote = input$quote)
        reference_matrix <- average_clusters(mat = data1(), metadata = data2()[[input$metadataCellType]], if_log = TRUE)
        head(reference_matrix)
    })
    
    output$clustify <- renderTable({
        #eh <- ExperimentHub()
        ## query
        #refs <- query(eh, "clustifyrdatahub")
        #refs
        #refs <- listResources(eh, "clustifyrdatahub")
        #benchmarkRef <- loadResources(eh, "clustifyrdatahub", input$dataHubReference)[[1]]
        
        ## use for classification of cell types
        #res <- clustify(
        #    input = data1(),
        #    metadata = data2()[[input$metadataCellType]],
        #    ref_mat = benchmarkRef,
        #    query_genes = pbmc_vargenes
        #)
    })
    
}

# Create Shiny app ----
shinyApp(ui, server)