# Run once from the project directory. Packages stay inside this project.
dir.create(".R-library", showWarnings = FALSE)
.libPaths(c(normalizePath(".R-library"), .Library), include.site = FALSE)
packages <- c("shiny", "ggplot2", "readxl", "survival", "nlme", "jsonlite", "nortest")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  # Repair incomplete transitive dependencies as well as top-level packages.
  available <- available.packages(repos = "https://cloud.r-project.org")
  dependencies <- tools::package_dependencies(packages, db = available,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
  needed <- setdiff(unique(c(packages, unlist(dependencies))), rownames(installed.packages()))
  install.packages(unique(c(missing, needed)), lib = .libPaths()[1], repos = "https://cloud.r-project.org")
}
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
message("Setup complete. Run source('launch.R') to open LabStats.")
