# Define UI for data upload app ----
ui <- dashboardPage(
  dashboardHeader(title = "Clustifyr RShiny App"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Load Matrix", tabName = "matrixLoad", icon = icon("th")),
      menuItem("Load Metadata", tabName = "metadataLoad", icon = icon("th")),
      menuItem("Choose cluster and ref column", tabName = "clusterRefCol", icon = icon("th"))
    )
    
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "dashboard",
        # js stuff ----
        useShinyjs(),
        tags$head(
          tags$script(HTML(js2))
        ),
        tags$head(tags$style(HTML('
            .skin-blue .sidebar .inactiveLink {
                color: black;
                opacity : 25%;
            }'
        ))),
        tags$head(tags$style(".inactiveLink {
                           pointer-events: none;
                           cursor: not-allowed;
                           }")),
        tags$head(tags$style(HTML('
            .skin-blue .sidebar .doneLink {
                color: green;
            }'
        ))),
        tags$head(tags$style(HTML('
            .skin-blue .sidebar .doneLink.active > a {
                color: green;
                border-left-color: green;
            }'
        ))),
        tags$head(tags$style(HTML('
            .skin-blue .sidebar .doneLink:hover {
                color: green;
                border-left-color: green;
            }'
        ))),
        
        # waiter stuff ----
        use_waiter(),
        
        # load example data ----
        actionButton("example",
                     "load example data",
                     icon = icon("space-shuttle")
        ),
        
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        
        # Horizontal line ----
        tags$hr(),
        
        # Input: Select separator ----
        radioButtons("sepMat", "Separator - Matrix",
                     choices = c(
                       Comma = ",",
                       Semicolon = ";",
                       Tab = "\t"
                     ),
                     selected = ","
        ),
        
        radioButtons("sepMeta", "Separator - Metadata",
                     choices = c(
                       Comma = ",",
                       Semicolon = ";",
                       Tab = "\t"
                     ),
                     selected = ","
        ),
        
        # Input: Select number of rows to display ----
        radioButtons("dispMat", "Display - Matrix",
                     choices = c(
                       Head = "head",
                       All = "all"
                     ),
                     selected = "head"
        ),
        radioButtons("dispMeta", "Display - Metadata",
                     choices = c(
                       Head = "head",
                       All = "all"
                     ),
                     selected = "head"
        )
      ),
      tabItem(
        tabName = "matrixLoad",
        h2("Load UMI Counts Matrix"),
        # Input: Select a file ----
        fileInput("file1", "Choose Matrix File",
                  multiple = TRUE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".xlsx",
                    ".tsv",
                    ".rds",
                    ".rda"
                  )
        ),
        
        # GEO id load ----
        actionButton("geo1", 
                     "from GEO id",
                     icon = icon("search")),
        
        actionButton("matrixPopup", "Display UMI Matrix in popup"),
        DTOutput("contents1"), # UMI Count Matrix
        tags$hr()
      ),
      tabItem(
        tabName = "metadataLoad",
        h2("Load Metadata table"),
        fileInput("file2", "Choose Metadata File",
                  multiple = FALSE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".xlsx",
                    ".tsv",
                    ".rds",
                    ".rda"
                  )
        ),
        
        # GEO id load ----
        actionButton("geo2", 
                     "from GEO id",
                     icon = icon("search")),
        
        actionButton("metadataPopup", "Display Metadata table in popup"),
        fluidRow(column(12, DTOutput('contents2'))),
        #DT::dataTableOutput("contents2"), # Metadata table
        tags$hr(),
        textOutput("colclicked"),
        
        h2("Choose cluster and reference column (cell types)"),
        selectInput("metadataCellType", "Cell Type Metadata Column:",
                    choice = list("")
        ),
        
        helpText("Choose cell type metadata column for average_clusters function"),
        hr(),
        selectInput("dataHubReference", "ClustifyrDataHub Reference:",
                    choices = list(
                      "ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs",
                      "ref_mouse.rnaseq", "ref_moca_main", "ref_immgen", "ref_hema_microarray",
                      "ref_cortex_dev", "ref_pan_indrop", "ref_pan_smartseq2",
                      "ref_mouse_atlas"
                    )
        ),
        helpText("Choose reference cell atlas for clustify function"),
        hr(),
        helpText("Choose cell reference for clustify function"),
      ),
      tabItem(
        tabName = "clusterRefCol",
        box(id = "box_clustifym",
            collapsible = TRUE,
            collapsed = TRUE,
            solidHeader = TRUE,
            status = "info",
            title = "clustifyr messages",
            htmlOutput("clustifym")),
        downloadButton("downloadReference", "Download reference matrix"),
        downloadButton("downloadClustify", "Download clustify matrix"),
        actionButton("uploadClustify", "Upload reference matrix"),
        
        DT::dataTableOutput("reference"), # Reference Matrix
        tags$hr(),
        DT::dataTableOutput("clustify"), # Clustify Matrix
        tags$hr(),
        plotOutput("hmap", height = "600px")
      )
    )
  )
)