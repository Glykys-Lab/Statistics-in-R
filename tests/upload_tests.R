# Simulate real uploaded files and Shiny events with independent per-file states.
# Cover switching/editing files, isolated invalid uploads, CSV decimal conventions,
# Excel import, project restoration and readable errors when required inputs are absent.
# Temporary CSV fixtures are written under artifacts/ and removed after import checks.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
upload_app <- new.env(parent = globalenv())
sys.source("app.R", envir = upload_app)
dir.create("artifacts", showWarnings = FALSE)
paths <- vapply(1:4, function(i) tempfile(paste0("upload_", i, "_"), tmpdir = "artifacts", fileext = ".csv"), character(1))
first <- upload_app$column_demo("groups")
second <- first; second$Control <- second$Control + 3
same <- first; same$Control[1] <- 11
write.csv(first, paths[1], row.names = FALSE)
write.csv(second, paths[2], row.names = FALSE)
write.csv(same, paths[3], row.names = FALSE)
writeLines("Control,Control\n1,2\n3,4", paths[4])
uploads <- data.frame(name = c("first.csv", "second.csv", "same.csv", "broken.csv"),
  size = file.info(paths)$size, type = "text/csv", datapath = paths, stringsAsFactors = FALSE)
shiny::testServer(upload_app$server, {
  session$setInputs(table_type = "column", source = "upload", file = uploads, upload_index = "1",
    layout = "wide", wide_columns = c("Control", "Drug A", "Drug B"), wide_id = "Subject", family = "many",
    matching = "independent", distribution = "normal", variance = "unequal", repeated = "complete",
    method = "welch_anova", y = "response", group = "group", confidence = ".95", adjustment = "holm",
    sphericity = "gg", decimal = ".", separator = ",", table_name = "first",
    plot_style = "points", error_style = "ci", plot_title = "", x_label = "", y_label = "")
  stopifnot(identical(raw_data(), first), length(project()$tables) == 0)
  session$setInputs(cell_edit = list(row = 1, column = 2, value = "11", revision = sheet_revision(), nonce = 1))
  session$setInputs(upload_index = "2")
  stopifnot(identical(raw_data(), second), is.null(saved()), is.null(working_table_id()))
  session$setInputs(upload_index = "1")
  stopifnot(raw_data()$Control[1] == 11, !analysis_unchanged())
  session$setInputs(run = 1)
  first_id <- working_table_id(); first_analysis <- saved()$id
  stopifnot(analysis_unchanged(), nrow(result()$data) == 24)
  session$setInputs(upload_index = "2", table_name = "second")
  stopifnot(identical(raw_data(), second), !analysis_unchanged())
  session$setInputs(run = 2)
  second_id <- working_table_id()
  stopifnot(first_id != second_id, length(project()$tables) == 2, length(project()$analyses) == 2,
    identical(project()$tables[[first_id]]$data, same), identical(project()$tables[[second_id]]$data, second))
  session$setInputs(upload_index = "1", table_name = "first")
  stopifnot(working_table_id() == first_id, saved()$id == first_analysis, analysis_unchanged())
  session$setInputs(run = 3)
  stopifnot(length(project()$analyses) == 2)
  session$setInputs(upload_index = "3", table_name = "same")
  stopifnot(identical(raw_data(), same), !analysis_unchanged(), is.null(working_table_id()))
  session$setInputs(run = 4)
  stopifnot(length(project()$tables) == 3, length(project()$analyses) == 3)
  session$setInputs(upload_index = "1", table_name = "first")
  session$setInputs(cell_edit = list(row = 1, column = 2, value = "12", revision = sheet_revision(), nonce = 2))
  stopifnot(!analysis_unchanged(), project()$tables[[first_id]]$data$Control[1] == 12,
    project()$analyses[[first_analysis]]$input_data$response[1] == 11)
  session$setInputs(upload_index = "2")
  stopifnot(identical(raw_data(), second), analysis_unchanged())
  session$setInputs(upload_index = "1", table_name = "first", run = 5)
  stopifnot(raw_data()$Control[1] == 12, length(project()$tables) == 3, length(project()$analyses) == 4)
  session$setInputs(upload_index = "4")
  bad <- tryCatch(raw_data(), error = identity)
  stopifnot(inherits(bad, "error"), length(project()$tables) == 3, length(project()$analyses) == 4)
  session$setInputs(upload_index = "1")
  stopifnot(raw_data()$Control[1] == 12, analysis_unchanged())
  reopened <- upload_app$read_project(upload_app$project_json(project()))
  stopifnot(length(reopened$tables) == 3, length(reopened$analyses) == 4,
    reopened$analyses[[first_analysis]]$input_data$response[1] == 11)
  # Saved-table edits and analyses must still be visible when returning to an upload.
  session$setInputs(project_open = list(kind = "tables", key = first_id, nonce = 1))
  session$setInputs(source = "workspace", table_name = "first")
  session$setInputs(cell_edit = list(row = 1, column = 2, value = "13", revision = sheet_revision(), nonce = 3))
  session$setInputs(run = 6)
  workspace_analysis <- saved()$id
  session$setInputs(source = "upload")
  stopifnot(working_table_id() == first_id, raw_data()$Control[1] == 13,
    saved()$id == workspace_analysis, analysis_unchanged(), length(project()$analyses) == 5)
  session$setInputs(run = 7)
  stopifnot(length(project()$tables) == 3, length(project()$analyses) == 5)
  # Replacing the upload selection clears only transient upload edits and references.
  session$setInputs(file = uploads[1, , drop = FALSE], upload_index = "1")
  stopifnot(identical(raw_data(), first), is.null(saved()), is.null(working_table_id()), length(project()$analyses) == 5)
})
semicolon_path <- tempfile("decimal_upload_", tmpdir = "artifacts", fileext = ".csv")
write.table(first, semicolon_path, sep = ";", dec = ",", row.names = FALSE)
semicolon <- data.frame(name = "decimal.csv", datapath = semicolon_path)
stopifnot(identical(upload_app$read_uploaded_table(semicolon, separator = ";", decimal = ","), first))
excel_path <- readxl::readxl_example("datasets.xlsx")
excel_file <- data.frame(name = "datasets.xlsx", datapath = excel_path)
stopifnot(nrow(upload_app$read_uploaded_table(excel_file, sheet = 1)) > 0)
unlink(c(paths, semicolon_path))
shiny::testServer(upload_app$server, {
  notifications <- list()
  send_notification <- session$sendNotification
  session$sendNotification <- function(type, message) {
    notifications[[length(notifications) + 1L]] <<- list(type = type, message = message)
    send_notification(type, message)
  }
  session$setInputs(table_type = "column", source = "upload", layout = "wide", run = 1)
  expect_visible_error <- function(text) {
    stopifnot(grepl(text, analysis_error(), fixed = TRUE),
      grepl(text, output$status$html, fixed = TRUE),
      grepl('role="alert"', output$status$html, fixed = TRUE),
      is.null(saved()), length(project()$analyses) == 0L, length(project()$tables) == 0L)
    shown <- tail(Filter(function(n) identical(n$type, "show"), notifications), 1L)[[1L]]$message
    stopifnot(identical(shown$id, "analysis-error"), grepl(text, as.character(shown$html), fixed = TRUE))
  }
  expect_visible_error("Choose a CSV, TSV or Excel file using Browse")
  session$setInputs(file = uploads[FALSE, , drop = FALSE], run = 2)
  expect_visible_error("Choose a CSV, TSV or Excel file using Browse")
  session$setInputs(source = "paste", run = 3)
  expect_visible_error("Paste a table and click Use pasted table")
  session$setInputs(source = "workspace", run = 4)
  expect_visible_error("Open a saved project table")
  # Example tables can be analyzed without uploading a file.
  session$setInputs(source = "demo", demo = "groups", wide_columns = character(), run = 5)
  expect_visible_error("Select the measurement columns")
  session$setInputs(wide_columns = c("Control", "Drug A", "Drug B"), wide_id = "Subject",
    family = "many", matching = "independent", distribution = "normal", variance = "unequal",
    repeated = "complete", method = "", y = "response", group = "group", confidence = ".95",
    adjustment = "holm", sphericity = "gg", decimal = ".", run = 6)
  expect_visible_error("Select a data table and complete the analysis settings")
  session$setInputs(method = "welch_anova", table_name = "Example", run = 7)
  stopifnot(is.null(analysis_error()), is_current(), length(project()$analyses) == 1L,
    !grepl('role="alert"', output$status$html, fixed = TRUE),
    identical(tail(notifications, 1L)[[1L]]$type, "remove"))
})
cat("PASS: Analyze shows readable banner and notification errors for missing inputs, then recovers without an upload\n")
cat("PASS: Multiple uploads keep file selection, unsaved edits, tables, analyses and snapshots separate; invalid files are isolated\n")
cat("PASS: Single-file uploads, decimal commas, Excel and project round-trips still work\n")
