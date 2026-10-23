# Data boundary for LabStats: examples, file/paste import, worksheet editing and
# conversion to analysis rows. Keep column names, missing values and table order.
# reproduction_script() embeds the data and analysis code in a standalone R file;
# it is an export containing measurements, not just a link back to this project.
demo_data <- function(name = "groups") {
  switch(name,
    groups = data.frame(subject = rep(sprintf("S%02d", 1:8), 3), group = rep(c("Control", "Drug A", "Drug B"), each = 8),
      response = c(9.1, 10.4, 8.8, 11.2, 9.8, 10.1, 8.9, 11.7, 12.1, 14.2, 11.9, 15.8, 13.7, 14.9, 12.5, 15.1, 15.4, 17.2, 14.8, 18.6, 16.9, 18.2, 15.6, 19.1)),
    paired = data.frame(subject = rep(sprintf("S%02d", 1:8), 2), group = rep(c("Before", "After"), each = 8),
      response = c(9.1, 10.4, 8.8, 11.2, 9.8, 10.1, 8.9, 11.7, 12.1, 14.2, 11.9, 15.8, 13.7, 14.9, 12.5, 15.1)),
    xy = data.frame(concentration = rep(1:8, each = 2), response = c(3.1, 4.2, 5.1, 6.7, 7.9, 6.8, 10.1, 9.2, 11.3, 13.1, 14.2, 12.9, 16.1, 17.4, 18.2, 17.3)),
    dose = { x <- rep(10^seq(-2, 2, length.out = 9), each = 3); data.frame(concentration = x, response = 5 + 95/(1 + (x/1.8)^1.2) + rep(c(-2, .3, 1.7), 9)) },
    contingency = data.frame(group = rep(c("Control", "Treatment"), each = 2), outcome = rep(c("Response", "No response"), 2), count = c(12, 28, 25, 15)),
    survival = data.frame(time = c(2, 3, 4, 5, 6, 7, 8, 9, 4, 6, 7, 8, 9, 10, 11, 12), event = c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0), group = rep(c("Control", "Treatment"), each = 8)),
    factorial = data.frame(group = rep(c("Control", "Treatment"), each = 8), sex = rep(rep(c("Female", "Male"), each = 4), 2), response = c(8, 10, 9, 11, 10, 13, 11, 12, 13, 15, 14, 17, 17, 18, 16, 20)),
    stop("Unknown example."))
}

validate_import <- function(data) {
  assert(ncol(data) > 0 && nrow(data) > 0, "The imported table is empty.")
  assert(!anyDuplicated(names(data)) && all(nzchar(trimws(names(data)))), "Use unique, nonempty column headers.")
  as.data.frame(data, check.names = FALSE)
}

read_uploaded_table <- function(file, separator = ",", decimal = ".", sheet = NULL) {
  # Read only the selected upload and Excel sheet; blank numeric fields become NA.
  # validate_import checks the resulting table before it enters the workspace.
  assert(is.data.frame(file) && nrow(file) == 1, "Select one uploaded file to analyze.")
  extension <- tolower(tools::file_ext(file$name))
  if (extension %in% c("xlsx", "xls")) {
    return(validate_import(readxl::read_excel(file$datapath,
      sheet = if (is.null(sheet)) 1 else sheet, .name_repair = "minimal")))
  }
  assert(extension %in% c("csv", "tsv", "txt"), "Choose a CSV, TSV, or Excel file.")
  separator <- if (extension %in% c("tsv", "txt")) "\t" else separator
  validate_import(read.table(file$datapath, header = TRUE, sep = separator, dec = decimal,
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"),
    comment.char = "", quote = "\""))
}

column_demo <- function(name = "groups") {
  d <- demo_data(if (name == "missing") "groups" else name)
  groups <- unique(d$group); subjects <- unique(d$subject)
  wide <- data.frame(Subject = subjects, check.names = FALSE)
  for (g in groups) wide[[g]] <- d$response[d$group == g][match(subjects, d$subject[d$group == g])]
  if (name == "missing") { wide[2, 2] <- NA_real_; wide[5, 3] <- NA_real_ }
  wide
}

paste_table <- function(text, decimal = ".") {
  assert(length(text) == 1 && nzchar(trimws(text)), "Paste a table with a header row first.")
  first <- strsplit(text, "\n", fixed = TRUE)[[1]][1]
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(";", first, fixed = TRUE)) ";" else ","
  validate_import(read.table(text = text, header = TRUE, sep = sep, dec = decimal, check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"), comment.char = "", quote = "\""))
}

edit_cell <- function(data, row, column, value, decimal = ".") {
  # Preserve the existing column type. Clearing a numeric cell records missing
  # data, while a nonnumeric edit must fail rather than alter the whole column.
  assert(length(row) == 1 && row %in% seq_len(nrow(data)) && length(column) == 1 && column %in% seq_len(ncol(data)), "Select a valid table cell.")
  text <- trimws(value)
  if (is.numeric(data[[column]])) {
    value <- if (text %in% c("", "NA", "NaN")) NA_real_ else suppressWarnings(as.numeric(if (decimal == ",") sub(",", ".", text, fixed = TRUE) else text))
    assert(text %in% c("", "NA", "NaN") || (length(value) == 1 && is.finite(value)), "Measurement cells need a number or a blank for missing data.")
  } else value <- if (text == "") NA_character_ else text
  data[row, column] <- value
  data
}

wide_to_long <- function(data, columns, id = "") {
  assert(length(columns) > 0 && all(columns %in% names(data)), "Select at least one measurement column.")
  assert(!id %in% columns, "The subject ID column cannot also be a measurement.")
  # The worksheet, rather than the order of selector clicks, defines the sequence.
  columns <- names(data)[names(data) %in% columns]
  for (name in columns) numeric_column(data, name)
  subject <- if (nzchar(id)) as.character(data[[id]]) else as.character(seq_len(nrow(data)))
  long <- do.call(rbind, lapply(columns, function(name) data.frame(subject = subject, group = name, response = data[[name]], stringsAsFactors = FALSE)))
  long$group <- factor(long$group, levels = columns)
  long
}

reproduction_script <- function(data, config, plot_options = list(), graph_kind = "raw", worksheet = NULL) {
  # Embed helpers and optional nested references so the export can run from a
  # different directory. Worksheet metadata is limited to what result sheets need.
  if (!is.null(worksheet)) worksheet <- worksheet[intersect(c("name", "layout", "columns", "data"), names(worksheet))]
  c("# LabStats reproducible analysis: data and functions are embedded below.",
    "# This file contains your data. Share it only with intended recipients.",
    "# Required packages: ggplot2; survival for survival; nlme for mixed effects; nortest for residual normality.",
    if (config$method %in% nested_methods) "# Linear mixed analyses additionally need lme4, lmerTest, emmeans, gridExtra, tidyr, dplyr, commonmark, svglite, ggbeeswarm and zip.",
    if (config$method %in% nested_methods) "# install.packages(c('ggplot2','nlme','lme4','lmerTest','emmeans','gridExtra','tidyr','dplyr','commonmark','svglite','ggbeeswarm','zip'))" else "# install.packages(c('ggplot2', 'survival', 'nlme', 'nortest'))", "",
    readLines(file.path("R", "nested.R"), warn = FALSE, encoding = "UTF-8"), "",
    if (config$method %in% nested_methods) c("embedded_nested_references <-", capture.output(dput(setNames(lapply(c("between-groups", "repeated"), function(n) readLines(file.path("R", paste0("nested-reference-", n, ".inc")), warn = FALSE, encoding = "UTF-8")), c("between-groups", "repeated")))), ""),
    readLines(file.path("R", "options.R"), warn = FALSE), "",
    readLines(file.path("R", "results.R"), warn = FALSE), "",
    readLines(file.path("R", "t_tests.R"), warn = FALSE), "",
    readLines(file.path("R", "analysis.R"), warn = FALSE), "",
    readLines(file.path("R", "plots.R"), warn = FALSE), "",
    "input_data <-", capture.output(dput(data)), "config <-", capture.output(dput(config)),
    "plot_options <-", capture.output(dput(plot_options)),
    "graph_kind <-", capture.output(dput(graph_kind)),
    "worksheet <-", capture.output(dput(worksheet)),
    "result <- run_analysis(input_data, config)",
    if (config$method %in% nested_methods) c("write_nested_bundle(result, 'linear_mixed_analysis.zip')", "writeLines(result$report_markdown, 'linear_mixed_analysis.md')", "writeLines(nested_html_report(result), 'linear_mixed_analysis.html')", "if (interactive() && requireNamespace('rstudioapi', quietly=TRUE) && rstudioapi::isAvailable()) rstudioapi::viewer(normalizePath('linear_mixed_analysis.html'))"),
    "writeLines(analysis_report(result, worksheet = worksheet), 'analysis_report.txt')",
    "tables <- result_export_tables(result, worksheet = worksheet, options = config)",
    "for (name in names(tables)) write.csv(tables[[name]], paste0(name, '.csv'), row.names = FALSE)",
    "figure <- make_graph(result, graph_kind, plot_options)",
    "save_result_figure(result, graph_kind, figure, 'figure.pdf', 'pdf')",
    "save_result_figure(result, graph_kind, figure, 'figure.png', 'png')",
    "writeLines(capture.output(sessionInfo()), 'session_info.txt')", "if (interactive()) draw_nested_graph(figure)")
}
