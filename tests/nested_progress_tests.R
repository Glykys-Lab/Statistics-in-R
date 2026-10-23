# Exercise real mixed-model runs and capture the progress messages sent to the
# browser. Progress must cover fitting AND saving, and always close after errors.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
progress_app <- new.env(parent = globalenv())
sys.source("app.R", envir = progress_app)

shiny::testServer(progress_app$server, {
  messages <- list()
  session$sendProgress <- function(type, message) {
    # Record whether the saved sheets are ready when the progress bar closes.
    messages[[length(messages) + 1L]] <<- list(type = type, message = message,
      ready = if (type == "close") !is.null(shiny::isolate(saved())) else FALSE)
  }
  session$setInputs(table_type = "nested", source = "demo", demo = "nested_repeated",
    nested_design = "repeated", method = "nested_repeated", nested_mouse = "Mouse",
    nested_slice = "Slice", nested_cell = "Cell", nested_columns = c("Baseline2", "Drug", "Wash"),
    nested_model_engine = "lmer", nested_slice_condition = "shared", nested_slice_variance = "shared",
    nested_check_variance_homogeneity = FALSE, nested_outlier_action = "none",
    nested_make_pub_figures = FALSE, nested_show_diagnostics = FALSE,
    nested_comparison_plot = "none", table_name = "Progress regression", decimal = ".")

  expect_progress <- function(success) {
    types <- vapply(messages, `[[`, character(1), "type")
    stopifnot(sum(types == "open") == 1L, sum(types == "close") == 1L,
      types[1] == "open", tail(types, 1L) == "close", tail(messages, 1L)[[1L]]$ready == success)
    updates <- lapply(Filter(function(m) m$type == "update", messages), `[[`, "message")
    values <- unlist(lapply(updates, `[[`, "value"))
    stopifnot(length(values) > 1L, all(diff(values) >= 0), values[1] == 0)
    if (success) {
      stopifnot(all(c(.05, .1, .15, .3, .45, .6, .7, .8, .9, .95, 1) %in% values),
        any(grepl("Fitting Model", vapply(updates, function(m) if (is.null(m$detail)) "" else m$detail, character(1)), fixed = TRUE)),
        is.null(analysis_error()), isTRUE(saved()$result$nested), is_current())
    } else stopifnot(tail(values, 1L) < 1, !is.null(analysis_error()), is.null(saved()))
  }

  session$setInputs(run = 1)
  expect_progress(TRUE)
  cat("PASS: Repeated-cell progress advances through fitting, figures, reports and saving\n")
  count <- length(project()$analyses)
  messages <- list()
  session$setInputs(run = 2)
  stopifnot(!length(messages), length(project()$analyses) == count)
  cat("PASS: An unchanged repeat click does not open a new progress display\n")

  session$setInputs(demo = "nested_between", nested_design = "between", method = "nested_between",
    nested_condition = "Condition", nested_value = "Value")
  messages <- list()
  session$setInputs(run = 3)
  expect_progress(TRUE)
  cat("PASS: Between-group progress advances and closes after saved results are ready\n")

  # A valid table with only one Mouse reaches candidate fitting, where the
  # random effect cannot be fitted. The failure must still dismiss progress.
  broken <- raw_data(); broken$Slice <- paste(broken$Mouse, broken$Slice, sep = "-")
  broken$Mouse <- "only-mouse"; edited(broken)
  messages <- list()
  count <- length(project()$analyses)
  session$setInputs(run = 4)
  expect_progress(FALSE)
  updates <- Filter(function(m) m$type == "update", messages)
  stopifnot(any(vapply(updates, function(m) identical(m$message$value, .3), logical(1))),
    length(project()$analyses) == count, grepl("could not be fitted", analysis_error(), fixed = TRUE))
  cat("PASS: A fitting failure closes progress, reports the error and saves no analysis\n")
})
