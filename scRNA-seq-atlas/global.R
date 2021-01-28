library(shiny)
library(shinyjs)
library(waiter)
library(dplyr)
library(readr)
library(tools)
library(clustifyr)
library(rsconnect)
library(ExperimentHub)
library(Seurat)
library(shinydashboard)
library(tidyverse)
library(data.table)
library(R.utils)
library(DT)
# library(ComplexHeatmap)

options(shiny.maxRequestSize = 1500 * 1024^2)
options(repos = BiocManager::repositories())
options(shiny.reactlog = TRUE)
options(DT.options = list(
  dom = "tp", 
  paging = TRUE,
  pageLength = 6,
  scrollX = TRUE
)
)

eh <- ExperimentHub()
refs <- query(eh, "clustifyrdatahub")
ref_dict <- refs$ah_id %>% setNames(refs$title)

js <- c(
  "table.on('click', 'td', function(){",
  "  var cell = table.cell(this);",
  "  var colindex = cell.index().column;",
  "  var colname = table.column(colindex).header().innerText;",
  "  Shiny.setInputValue('column_clicked', colname);",
  "});"
)