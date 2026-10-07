if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
if (!requireNamespace("shiny", quietly = TRUE)) stop("Run source('setup.R') first.")
shiny::runApp(".", host = "127.0.0.1", port = 7319, launch.browser = TRUE)
