# Define server logic to read selected file ----
server <- function(input, output, session) {
  # reactive file location to make interactivity easier
  rv <- reactiveValues()
  rv$matrixloc <- NULL
  rv$metaloc <- NULL
  rv$step <- 0
  rv$clustifym <- "clustifyr not yet run"
  rv$lastgeo <- "GSE151974"#GSE113049"
  rv$ref <- NULL
  rv$ref_visited <- 0
  rv$ref_link <- NULL
  rv$res_visited <- 0
  rv$obj <- NULL

  # waiter checkpoints
  w1 <- Waiter$new(
    id = "contents1",
    html = tagList(
      spin_flower(),
      h4("Matrix loading..."),
      h4("")
    )
  )

  w2 <- Waiter$new(
    id = "contents2",
    html = tagList(
      spin_flower(),
      h4("Metadata loading..."),
      h4("")
    )
  )

  w3 <- Waiter$new(
    id = "reference",
    html = tagList(
      spin_flower(),
      h4("Reference building..."),
      h4("")
    )
  )

  w4 <- Waiter$new(
    id = "clustify",
    html = tagList(
      spin_flower(),
      h4("Clustifyr running..."),
      h4("")
    )
  )

  w5 <- Waiter$new(
    id = "hmap",
    html = tagList(
      spin_flower(),
      h4("Heatmap drawing..."),
      h4("")
    )
  )

  w6 <- Waiter$new(id = "modalgeo",
                   html = tagList(
                     spin_flower(),
                     h4("Info fetching..."),
                     h4("")
                   ))

  w7 <- Waiter$new(id = "modalfiles",
                   html = tagList(
                     spin_flower(),
                     h4("File previewing..."),
                     h4("")
                   ))

  w8 <- Waiter$new(
    id = "contents3",
    html = tagList(
      spin_flower(),
      h4("Reference loading..."),
      h4("")
    )
  )

  data1 <- reactive({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$file1) | !is.null(rv$matrixloc)) {
      if (!is.null(input$file1)) {
        rv$matrixloc <- input$file1
      }
      file <- rv$matrixloc

    if (!is.null(file)) {
      w1$show()
      print(file)
    }

    fileTypeFile1 <- tools::file_ext(file$datapath)
    req(file)
    if (str_to_lower(fileTypeFile1) == "rds") {
      df1 <- readRDS(file$datapath)
      if (any(class(df1) %in% c("SingleCellExperiment", "Seurat"))) {
        rv$obj <- df1
        df1 <- object_data(rv$obj, "data")
      }
    } else if (str_to_lower(fileTypeFile1) == "rdata") {
      df1 <- load_rdata(file$datapath)
      if (any(class(df1) %in% c("SingleCellExperiment", "Seurat"))) {
        rv$obj <- df1
        df1 <- object_data(df1, "data")
      }
    } else {
      df1 <- fread(file$datapath)
    }
    } else if (!is.null(rv$obj)) {
      df1 <- object_data(rv$obj, "data")
    } else {
      return(NULL)
    }

    df1 <- df1 %>% as.data.frame()
    if (!has_rownames(df1)) {
        rownames(df1) <- df1[, 1]
        df1[, 1] <- NULL
    }

    w1$hide()
    df1
  })

  data2 <- reactive({
    if (!is.null(input$file2) | !is.null(rv$metaloc)) {
      if (!is.null(input$file2)) {
        rv$metaloc <- input$file2
      }
      file <- rv$metaloc

    if (!is.null(file)) {
      w2$show()
      print(file)
    }

    fileTypeFile2 <- tools::file_ext(file$datapath)
    req(file)
    if (str_to_lower(fileTypeFile2) == "rds") {
      df2 <- readRDS(file$datapath)
      if (any(class(df2) %in% c("SingleCellExperiment", "Seurat"))) {
        rv$obj <- df2
        df2 <- object_data(df2, "meta.data")
      }
    } else if (str_to_lower(fileTypeFile2) == "rdata") {
      df2 <- load_rdata(file$datapath)
      if (any(class(df2) %in% c("SingleCellExperiment", "Seurat"))) {
        rv$obj <- df2
        df2 <- object_data(rv$obj, "meta.data")
      }
    } else {
      df2 <- fread(file$datapath)
    }
    } else if (!is.null(rv$obj)) {
      df2 <- object_data(rv$obj, "meta.data")
    } else {
      return(NULL)
    }

    df2 <- df2 %>% as.data.frame()
    if (!has_rownames(df2)) {
      rownames(df2) <- df2[, 1]
      df2[, 1] <- NULL
    }

    w2$hide()
    df2
  })

  data3a <- reactive({
    if (!is.null(input$file3)) {
      rv$ref <- input$file3
    }
    if (rv$ref == "built-in") {
      return(NULL)
    }
    file <- rv$ref

    if (!is.null(file)) {
      w8$show()
      print(file)
    }

    fileTypeFile3 <- tools::file_ext(file$datapath)
    req(file)

    if (str_to_lower(fileTypeFile3) == "rds") {
      df3 <- readRDS(file$datapath) %>% as.data.frame()
    } else if (str_to_lower(fileTypeFile3) == "rdata") {
      df3 <- load_rdata(file$datapath) %>% as.data.frame()
    } else {
      df3 <- fread(file$datapath) %>% # , header = input$header, sep = input$sepMat) %>%
        as.data.frame()
    }

    if (!has_rownames(df3)) {
      rownames(df3) <- df3[, 1]
      df3[, 1] <- NULL
    }

    w8$hide()
    df3
  })

  output$contents1 <- DT::renderDataTable({
    if (is.null(rv$matrixloc) & is.null(rv$obj)) {
      return(df1 <- data.frame(`nodata` = rep("", 6)))
    } else {
      df1 <- data1()
    }

    # file 1
    if (input$dispMat == "head") {
      cols <- ncol(df1)
      df1 <- df1[, 1:min(cols, 5)]
      return(head(df1))
    }
    else {
      return(df1)
    }
  })

  output$contents2 <- DT::renderDataTable({
    if (is.null(rv$metaloc) & is.null(rv$obj)) {
      return(df2 <- data.frame(`nodata` = rep("", 6)))
    } else {
      df2 <- data2()
    }

    updateSelectInput(session, "metadataCellType",
                      choices = c("", colnames(df2)),
                      selected = ""
    )

    # file 2
    if (input$dispMeta == "head") {
      return(head(df2))
    }
    else {
      return(df2)
    }

  },
  callback = DT::JS(js),
  selection = list(target = 'column', mode = "single"))

  output$colclicked <- renderUI({
    if (is.null(input[["column_clicked"]])) {
      "please select cluster column in drop-down menu, or click in the table"
    } else {
      input$metadataCellType
    }
  })

  observeEvent(input[["column_clicked"]], {
    updateSelectInput(session, "metadataCellType",
      selected = input[["column_clicked"]]
    )
  })

  output$ref_summary <- renderUI({
    HTML(paste0("cell types: ", ncol(data3()),
           "<br>",
           "genes: ", nrow(data3())))
  })

  data3b <- reactive({
    w8$show()
    rv$ref <- "built-in"
    ref <- refs[[ref_dict[input$dataHubReference]]]
    rv$ref_link <- refs_meta[ref_dict[input$dataHubReference], ] %>% pull(sourceurl)
    w8$show()

    ref
  })

  data3 <- reactive({
    b <- data3b()
    a <- data3a()
    if (is.null(rv$ref)) {
      return(data.frame(`nodata` = rep("", 6)))
    } else if (rv$ref == "built-in") {
      df3 <- b
    } else {
      df3 <- a
    }
    df3
  })

  output$contents3 <- DT::renderDataTable({
    df3 <- data3()

    # file 3
    if (input$dispMat == "head") {
      return(head(df3))
    }
    else {
      return(df3)
    }
  })

  observeEvent(input$matrixPopup, {
    showModal(modalDialog(
      tags$caption("Matrix table"),
      DT::renderDataTable({
        matrixRender <- head(data1())
        DT::datatable(matrixRender)
      }),
      easyClose = TRUE
    ))
  })

  observeEvent(input$metadataPopup, {
    showModal(modalDialog(
      tags$caption("Metadata table"),
      DT::renderDataTable({
        matrixRender <- head(data2())
        DT::datatable(matrixRender)
      }),
      easyClose = TRUE
    ))
  })

  data_avg <- reactive({
    if (input$metadataCellType == "") {
      return(NULL)
    }
    w3$show()
    reference_matrix <- average_clusters(mat = data1(), metadata = data2()[[input$metadataCellType]], if_log = FALSE)
    w3$hide()
    reference_matrix
  })

  dataClustify <- reactive({
    if (input$metadataCellType == "") {
      return(NULL)
    }
    w4$show()
    benchmarkRef <- data3()

    if (!is.null(rv$obj)) {
      message("Single cell object detected")
      matrixSeuratObject <- rv$obj
      if (any(class(matrixSeuratObject) == "SingleCellExperiment")) {
        matrixSeuratObject <- as.Seurat(matrixSeuratObject)
      }
    } else {
      UMIMatrix <- data1()
      matrixSeuratObject <- CreateSeuratObject(counts = UMIMatrix, project = "Seurat object matrix", min.cells = 0, min.features = 0)
    }
    if (VariableFeatures(matrixSeuratObject) %>% length() == 0) {
      matrixSeuratObject <- FindVariableFeatures(matrixSeuratObject, selection.method = "vst", nfeatures = 2000)
    } else {
      message("Using variable genes in object")
    }

    metadataCol <- data2()[[input$metadataCellType]]
    # use for classification of cell types
    messages <<- capture.output(
      res <- clustify(
        input = object_data(matrixSeuratObject, "data"),
        metadata = metadataCol,
        ref_mat = benchmarkRef,
        query_genes = VariableFeatures(matrixSeuratObject),
        verbose = TRUE
        ),
      type = "message"
    )
    rv$clustifym <<- messages

    w4$hide()
    res
  })

  output$reference <- DT::renderDataTable({
    if (rv$res_visited == 1) {
      message("empty")
      return(df1 <- data.frame(`nodata` = rep("", 6)))
    }
    reference_matrix <- data_avg()
    rownames_to_column(as.data.frame(reference_matrix), input$metadataCellType)
  })

  output$clustify <- DT::renderDataTable({
    if (rv$res_visited == 1) {
      message("empty2")
      return(df1 <- data.frame(`nodata` = rep("", 6)))
    }
    res <- dataClustify()
    rownames_to_column(as.data.frame(res), input$metadataCellType)
  })

  # Make plots such as heat maps to compare benchmarking with clustify with actual cell types

  output$hmap <- renderPlot({
    if (input$metadataCellType == "") {
      return(NULL)
    }

    # could expose as an option
    cutoff_to_display <- 0.5
    tmp_mat <<- dataClustify()

    if (!is.null(tmp_mat)) {
      w5$show()
    }
    tmp_mat <- tmp_mat[, colSums(tmp_mat > 0.5) > 1]
    plot_hmap(tmp_mat)
  })

  referenceDownload <- reactive({
    avgMatrix <- data_avg()
  })

  clustifyDownload <- reactive({
    clustifyMatrix <- dataClustify()
  })

  output$downloadReference <- downloadHandler(
    filename = function() {
      paste("reference-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(referenceDownload(), file, quote = FALSE)
    }
  )
  output$downloadClustify <- downloadHandler(
    filename = function() {
      cat("clustify-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(clustifyDownload(), file, quote = FALSE)
    }
  )

  # load example data
  observeEvent(
    input$example,
    {
      message("loading prepackaged data")
      rv$matrixloc <- list(datapath = "../data/example-input/matrix.csv")
      rv$metaloc <- list(datapath = "../data/example-input/meta-data.csv")
      updateTabItems(session, "tabs", "metadataLoad")
    }
  )

  output$clustifym <- renderUI(
    HTML(paste0(c(rv$clustifym, ""), collapse = "<br/><br/>"))
  )

  # modal for GEO id
  observeEvent(
    input$geo1 | input$geo2,
    showModal(modalDialog(
      div(id = "modalgeo",
          textInput("geoid", "query GEO id", value = rv$lastgeo),
          actionButton("geogo", "Fetch file info", icon = icon("eye"))
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = NULL
    )),
    ignoreInit = T
  )

  observeEvent(
    input$geogo,
    {
      w6$show()
      rv$lastgeo <- input$geoid
      rv$links <- list_geo(rv$lastgeo)
      print(rv$links)
      links2 <- cbind(rv$links %>% mutate(size = map(link, get_file_size)) %>% select(-link),
                      button = sapply(1:nrow(rv$links), make_button("tbl1")),
                      stringsAsFactors = FALSE) %>%
        data.table::data.table()
      links2 <- links2 %>%
        DT::datatable(options = list(
          dom = "ftp",
          searchHighlight = TRUE,
          paging = TRUE,
          pageLength = 5,
          scrollY = FALSE),
          escape = ncol(links2) - 1, fillContainer = TRUE)
      url <- str_c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$geoid)
      w6$hide()
      showModal(modalDialog(
        size = "l",
        div(id = "modalfiles",
            DT::renderDataTable(links2)
        ),
        easyClose = TRUE,
        fade = FALSE,
        footer = tagList(
          actionButton("geopage", label = "Go to GEO page",
                       onclick = paste0('window.open("',
                                        url,
                                        '", "_blank")'),
                       icon = icon("link")),
          actionButton("email", label = "Email author for missing data",
                       onclick = paste0('location.href="',
                                        prep_email(rv$lastgeo),
                                        '"'),
                       icon = icon("envelope-open-text")))
      ))
    }
  )

  observeEvent(input[["button"]], {
    print("button")
    w7$show()
    splitID <- strsplit(input[["button"]], "_")[[1]]
    tbl <- splitID[2]
    row <- splitID[3]
    rv$loadinglink <<- rv$links$link[as.numeric(row)]
    print(rv$links)

    # if tar, read a file list
    if (str_detect(rv$links$link[as.numeric(row)], "/GSE[0-9]+_RAW.tar")) {
      fullb <- F
      previewdata <- preview_link(get_tar(rv$links$link[as.numeric(row)]))
    } else {
      fullb <- T
      previewdata <- preview_link(rv$links$link[as.numeric(row)])
    }

    if (is.null(previewdata)) {
      fullb <- F
      previewdata <- data.frame(unreadable = rep("", 4))
    } else {
      cols <- ncol(previewdata)
      previewdata <-previewdata[, 1:min(cols, 5)]
    }

    w7$hide()
    showModal(modalDialog(
      size = "l",
      div(id = "modalback",
          title = "preview",
          DT::renderDataTable(previewdata),
          if (fullb) {
            actionButton("full", "Start full loading", icon = icon("running"))
          } else {
            disabled(actionButton("full", "Start full loading", icon = icon("running")))
          },
          actionButton("back", "Back to file list", icon = icon("step-backward"))
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = NULL
    ))
  })

  observeEvent(input$full, {
    print(rv$loadinglink)
    if (input[["activeTab"]] == "matrixLoad") {
      rv$matrixloc <- list(datapath = rv$loadinglink)
    } else {
      rv$metaloc <- list(datapath = rv$loadinglink)
    }
    print(rv$loadinglink)
    removeModal()
  })

  observeEvent(input$back, {
    links2 <- cbind(rv$links %>% mutate(size = map(link, get_file_size)) %>% select(-link),
                    button = sapply(1:nrow(rv$links), make_button("tbl1")),
                    stringsAsFactors = FALSE) %>%
      data.table::data.table()
    links2 <- links2 %>%
      DT::datatable(options = list(
        dom = "ftp",
        searchHighlight = TRUE,
        paging = TRUE,
        pageLength = 5,
        scrollY = FALSE),
        escape = ncol(links2)-1, fillContainer = TRUE)
    url <- str_c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$geoid)
    showModal(modalDialog(
      size = "l",
      div(id = "modalfiles",
          DT::renderDataTable(links2)
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = tagList(
        actionButton("geopage", label = "Go to GEO page",
                     onclick = paste0('window.open("',
                                      url,
                                      '", "_blank")'),
                     icon = icon("link")),
        actionButton("email", label = "Email author for missing data",
                     onclick = paste0('location.href="',
                                      prep_email(rv$lastgeo),
                                      '"'),
                     icon = icon("envelope-open-text")))
    ))
  })

  # disable menu at load
  addCssClass(selector = "a[data-value='clustifyres']", class = "inactiveLink")
  addCssClass(selector = "ul li:eq(4)", class = "inactiveLink")

  # check if data is loaded
  observeEvent((!is.null(data1())) + (!is.null(data2())) + (!is.null(data3())) +
                 (!is.null(input$metadataCellType)) +
                 (input$metadataCellType != ""), {
    if ((!is.null(data1())) + (!is.null(data2())) + (!is.null(data3())) +
        (!is.null(input$metadataCellType)) +
        (input$metadataCellType != "") == 5) {
      removeCssClass(selector = "a[data-value='clustifyres']", class = "inactiveLink")
      removeClass(selector = "ul li:eq(4)", class = "inactiveLink")
    }
  })

  observeEvent(data1(), {
    if (!is.null(data1())) {
      addCssClass(selector = "a[data-value='matrixLoad']", class = "doneLink")
      addClass(selector = "ul li:eq(1)", class = "doneLink")
    }
  })

  observeEvent(input$metadataCellType, {
    if (input$metadataCellType != "") {
      addCssClass(selector = "a[data-value='metadataLoad']", class = "doneLink")
      addClass(selector = "ul li:eq(2)", class = "doneLink")
    }
  })

  observeEvent(input[["activeTab"]], {
    if (input[["activeTab"]] == "clusterRef") {
      rv$ref_visited <<- 1
    } else if (input[["activeTab"]] == "clustifyres") {
      if (rv$res_visited == 0)
      rv$res_visited <<- 1
    }
  })

  observeEvent(rv$ref_visited, {
    if (rv$ref_visited == 1 & !is.null(data3())) {
      addCssClass(selector = "a[data-value='clusterRef']", class = "doneLink")
      addClass(selector = "ul li:eq(3)", class = "doneLink")
    }
  })

  observeEvent(rv$res_visited, {
    if (rv$res_visited == 1) {
      print("change")
      rv$res_visited <- 2
    }
  }, ignoreInit = FALSE)

  observeEvent(rv$ref_link, {
    runjs(paste0("document.getElementById('ref_linkgo').onclick = function() {
           window.open('", rv$ref_link, "', '_blank');};"))
  })
}
