library(shiny)
library(dplyr)
library(readr)
library(tools)
library(clustifyr)
library(rsconnect)
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
                         choices = c(Tab = "\t"),
                         selected = ","),
            
            # Input: Select number of rows to display ----
            radioButtons("disp", "Display",
                         choices = c(Head = "head",
                                     All = "all"),
                         selected = "head")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            tableOutput("contents")
            
        )
        
    )
)

# Define server logic to read selected file ----
server <- function(input, output) {
    
    output$contents <- renderTable({
        
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
            }
            else if (fileTypeFile1 == "tsv")
            {
                df1 <- read_tsv(file$datapath,
                                header = input$header,
                                quote = input$quote)
            }
            else
            {
                df1 <- load(file$datapath)
            }
                
            #file 1
            if(input$disp == "head") {
                return(head(df1))
            }
            else {
                return(df1)
            }
        }
    )
    
    output$file2Contents <- renderTable({
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
        
    output$reference <- renderTable({
        req(input$file1)
        req(input$file2)
        fileTypeFile1 <- file_ext(input$file1)
        fileTypeFile2 <- file_ext(input$file2)
        df1 <- read.csv(input$file1$datapath,
                        header = input$header,
                        sep = input$sep,
                        quote = input$quote)
        df2 <- read.csv(input$file2$datapath,
                        header = input$header,
                        sep = input$sep,
                        quote = input$quote)
        reference_matrix <- average_clusters(mat = df1, metadata = df2$cellCol, if_log = TRUE)
        head(reference_matrix)
    })
    
}

# Create Shiny app ----
shinyApp(ui, server)