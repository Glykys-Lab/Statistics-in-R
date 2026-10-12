# Start the local Shiny app from the project directory after setup.R has completed.
# Use the project's package library first; bind to localhost and open port 7319
# in the browser. The R process must stay running while the app is in use.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
if (!requireNamespace("shiny", quietly = TRUE)) stop("Run source('setup.R') first.")
shiny::runApp(".", host = "127.0.0.1", port = 7319, launch.browser = TRUE)
