# Run once from the project directory. Packages stay inside this project.
# Install missing app and nested-model dependencies from CRAN into .R-library.
# Excluding the site library makes this check reflect what the project itself
# can load. This installs compatible packages; it does not lock package versions.
dir.create(".R-library", showWarnings = FALSE)
.libPaths(c(normalizePath(".R-library"), .Library), include.site = FALSE)
packages <- c("shiny", "ggplot2", "readxl", "survival", "nlme", "jsonlite", "nortest",
  "lme4", "lmerTest", "emmeans", "gridExtra", "tidyr", "dplyr", "commonmark", "svglite", "ggbeeswarm", "zip", "rstudioapi")
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
