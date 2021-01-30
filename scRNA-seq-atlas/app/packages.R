# require(rsconnect)
# options(repos = BiocManager::repositories())
# pack <- rsconnect::appDependencies()
# write.csv(pack, "data/dependencies.csv")
install_packages_app <- function() {
  options(repos=structure(c(CRAN="http://cloud.r-project.org/")))
  pack <- read.csv("data/dependencies.csv")
  inst <- installed.packages()[,"Package"]
  pack_need <- !unlist(vapply(pack$package, function(x) {x %in% inst}, logical(1)))
  pack2 <- pack[pack_need, ]
  if (!("BiocManager" %in% inst)) {
    install.packages("BiocManager")
  }
  if (!("remotes" %in% inst)) {
    install.packages("remotes")
  }
  while (nrow(pack2) != 0) {
    message("Installing ", nrow(pack2), " package(s)...")
    if (pack2$source[1] == "CRAN") {
      install.packages(pack2$package[1], dependencies = TRUE, ask = FALSE)
    } else if (pack2$source[1] == "Bioconductor") {
      BiocManager::install(pack2$package[1], ask = FALSE, update = FALSE)
    } else {
      message("Installing latest GitHub devel version of clustifyr")
      remotes::install_github("rnabioco/clustifyr", dependencies = TRUE, upgrade = FALSE)
    }
    inst <- installed.packages()[,"Package"]
    pack_need <- !unlist(vapply(pack2$package, function(x) {x %in% inst}, logical(1)))
    pack2 <- pack2[pack_need, ]
  }
  message("All dependent packages installed")
}
