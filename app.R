# Application entry point: load the analysis modules and connect the Shiny UI to
# per-session data, analysis snapshots, project navigation, graphs and downloads.
# Launch from the project directory through launch.R or Shiny's Run App command.
# Calculation functions live in R/ so they can also run without the browser.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
required <- c("shiny", "ggplot2", "readxl", "survival", "nlme", "jsonlite", "nortest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Run source('setup.R') first. Missing packages: ", paste(missing, collapse = ", "))
library(shiny)
dir.create(".cache/sass", recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = normalizePath(".cache"))
source("R/analysis.R", local = TRUE)
source("R/data.R", local = TRUE)
source("R/plots.R", local = TRUE)
source("R/project.R", local = TRUE)
options(shiny.maxRequestSize = 25 * 1024^2)

source("R/ui.R", local = TRUE)

server <- function(input, output, session) {
  # Each browser session owns these values; edits must not overwrite an archived
  # result or leak into another session's worksheet.
  saved <- reactiveVal(NULL)
  analysis_error <- reactiveVal(NULL)
  edited <- reactiveVal(NULL)
  pasted <- reactiveVal(NULL)
  sheet_revision <- reactiveVal(0L)
  upload_states <- reactiveVal(list())
  source("R/workspace.inc", local = TRUE)
  source("R/options_server.inc", local = TRUE)
  source("R/nested_server.inc", local = TRUE)
  output$distribution_ui <- renderUI({
    choices <- c("Approximately Gaussian (parametric)" = "normal", "Rank-based (nonparametric)" = "rank")
    if (isTRUE(input$family %in% c("one", "multiple_t"))) choices <- choices[1]
    if (isTRUE(input$family %in% c("one", "two", "multiple_t")) || (identical(input$matching, "matched") && identical(input$family, "many"))) choices <- c(choices, "Lognormal (analyze log responses)" = "lognormal")
    selected <- pending_preferences()$distribution
    if (is.null(selected)) selected <- isolate(input$distribution)
    if (is.null(selected) || !selected %in% choices) selected <- "normal"
    radioButtons("distribution", "Distribution assumption", choices, selected = selected)
  })
  # Keep this input mounted even when its tab is hidden. Rebuilding it as an
  # output leaves a faded, stale control while hidden outputs are suspended.
  correction_update <- reactiveVal(NULL)
  observe({
    choices <- comparison_correction_choices(input$method, input_or("comparisons", "all"), input_or("error_control", "fwer"))
    selected <- pending_preferences()$adjustment
    if (is.null(selected)) selected <- isolate(input$adjustment)
    if (is.null(selected) || !selected %in% choices) selected <- unname(choices[1])
    update <- list(choices = choices, selected = selected)
    if (!identical(update, isolate(correction_update())) || !identical(selected, isolate(input$adjustment))) {
      updateSelectInput(session, "adjustment", choices = choices, selected = selected)
      correction_update(update)
    }
  })
  output$control_ui <- renderUI({
    d <- data(); req(input$group %in% names(d)); selectInput("control", "Control group", levels(group_column(d, input$group)))
  })
  output$demo_ui <- renderUI({
    choices <- switch(input$table_type,
      column = c("Three independent groups" = "groups", "Matched before / after" = "paired", "Repeated measures with missing values" = "missing", "Multiple t-tests across markers" = "multiple_tests"),
      xy = c("XY measurements" = "xy"),
      nested = c("Repeated cells within slices and mice" = "nested_repeated", "Between groups of mice" = "nested_between"),
      other = c("Three groups" = "groups", "Matched before / after" = "paired", "XY measurements" = "xy", "Dose-response" = "dose", "Category counts" = "contingency", "Time to event" = "survival", "Two-factor experiment" = "factorial"))
    selectInput("demo", "Example table", choices)
  })
  # Keep example selection available when the data sheet is hidden.
  outputOptions(output, "demo_ui", suspendWhenHidden = FALSE)
  observeEvent(input$table_type, {
    layout <- if (identical(input$source, "workspace") && !is.null(opened_table())) opened_table()$layout else if (input$table_type == "column") "wide" else "long"
    updateSelectInput(session, "layout", selected = layout)
    updateTabsetPanel(session, "tabs", selected = "data")
  })
  observeEvent(list(input$table_type, input$family, input$matching, input$distribution, input$variance, input$repeated, input$nested_design), {
    req(input$table_type)
    choices <- design_methods(input$table_type, if (is.null(input$family)) "many" else input$family, identical(input$matching, "matched"), if (is.null(input$distribution)) "normal" else input$distribution, if (is.null(input$variance)) "unequal" else input$variance, if (is.null(input$repeated)) "complete" else input$repeated)
    if (identical(input$table_type, "nested")) choices <- choices[choices == if (identical(input$nested_design, "between")) "nested_between" else "nested_repeated"]
    selected <- isolate(input$method)
    if (is.null(selected) || !selected %in% choices) selected <- unname(choices[1])
    updateSelectInput(session, "method", choices = choices, selected = selected)
  })
  observeEvent(input$apply_paste, {
    tryCatch({ pasted(paste_table(input$paste_text, input$decimal)); showNotification("Pasted table loaded.", type = "message") },
      error = function(e) showNotification(conditionMessage(e), type = "error"))
  })
  observeEvent(input$file, {
    upload_states(list())
    updateSelectInput(session, "upload_index", selected = "1")
  }, priority = 200)
  uploaded_file <- reactive({
    req(input$file, nrow(input$file) > 0)
    index <- suppressWarnings(as.integer(input$upload_index))
    if (length(index) != 1 || is.na(index) || !index %in% seq_len(nrow(input$file))) index <- 1L
    input$file[index, , drop = FALSE]
  })
  output$uploaded_file_ui <- renderUI({
    req(input$source == "upload", input$file, nrow(input$file) > 1)
    labels <- make.unique(input$file$name, sep = " #")
    selected <- isolate(input$upload_index)
    if (is.null(selected) || !selected %in% as.character(seq_len(nrow(input$file)))) selected <- "1"
    selectInput("upload_index", "File to analyze", setNames(as.character(seq_len(nrow(input$file))), labels), selected = selected)
  })
  file_is_excel <- reactive({ tolower(tools::file_ext(uploaded_file()$name)) %in% c("xlsx", "xls") })
  upload_key <- reactive({
    # Parsing settings identify an upload state too: the same file interpreted
    # with a different sheet or decimal mark must not reuse incompatible edits.
    file <- uploaded_file()
    if (file_is_excel()) return(paste(file$datapath, "excel", if (is.null(input$sheet)) "" else input$sheet, sep = "\r"))
    extension <- tolower(tools::file_ext(file$name))
    separator <- if (extension %in% c("tsv", "txt")) "\t" else if (is.null(input$separator)) "," else input$separator
    paste(file$datapath, separator,
      if (is.null(input$decimal)) "." else input$decimal,
      sep = "\r")
  })
  output$sheet_ui <- renderUI({
    req(input$source == "upload", file_is_excel())
    sheets <- tryCatch(readxl::excel_sheets(uploaded_file()$datapath), error = function(e) character())
    selectInput("sheet", "Excel worksheet", sheets)
  })
  base_data <- reactive({
    # Keep the source table separate from edits and the mapped analysis table.
    # This lets file switching restore edits while exports retain the original IDs.
    if (input$source == "workspace") { req(opened_table()); return(opened_table()$data) }
    if (input$source == "demo") {
      req(input$demo)
      if (input$table_type == "nested") return(nested_demo(identical(input$demo, "nested_between")))
      if (identical(input$demo, "multiple_tests")) return(multiple_t_demo())
      if (input$table_type == "column") { req(input$demo %in% c("groups", "paired", "missing")); return(column_demo(input$demo)) }
      return(demo_data(input$demo))
    }
    if (input$source == "blank") {
      req(input$blank_rows)
      rows <- max(3L, min(100L, as.integer(input$blank_rows)))
      if (input$table_type == "xy") return(data.frame(X = rep(NA_real_, rows), Y = rep(NA_real_, rows)))
      req(input$blank_groups)
      groups <- max(2L, min(12L, as.integer(input$blank_groups)))
      return(as.data.frame(setNames(rep(list(rep(NA_real_, rows)), groups), paste("Group", LETTERS[seq_len(groups)])), check.names = FALSE))
    }
    if (input$source == "paste") { req(pasted()); return(pasted()) }
    file <- uploaded_file()
    tryCatch({
      if (file_is_excel()) req(input$sheet)
      read_uploaded_table(file, if (is.null(input$separator)) "," else input$separator,
        if (is.null(input$decimal)) "." else input$decimal, input$sheet)
    }, error = function(e) validate(need(FALSE, conditionMessage(e))))
  })
  observeEvent(base_data(), {
    sheet_revision(isolate(sheet_revision()) + 1L)
    if (identical(input$source, "upload")) {
      state <- upload_states()[[upload_key()]]
      table <- if (is.null(state$table_id)) NULL else project()$tables[[state$table_id]]
      edited(if (is.null(table)) state$data else table$data)
      active_table_id(state$table_id); working_table_id(state$table_id)
      select_sheet("tables", state$table_id)
      restored <- if (is.null(state$analysis_id)) NULL else project()$analyses[[state$analysis_id]]
      saved(restored)
      archive_view(FALSE)
      updateTextInput(session, "table_name", value = if (is.null(state$table_id)) tools::file_path_sans_ext(uploaded_file()$name) else project()$tables[[state$table_id]]$name)
      if (!is.null(restored)) {
        updateTextInput(session, "analysis_name", value = restored$name)
        updateSelectInput(session, "graph_kind", selected = "raw")
      }
    } else {
      edited(NULL)
      if (input$source != "workspace") { active_table_id(NULL); working_table_id(NULL); selected_sheet(NULL) }
    }
  }, priority = 100)
  raw_data <- reactive({ d <- base_data(); if (is.null(edited())) d else edited() })
  observeEvent(input$cell_edit, {
    cell <- input$cell_edit
    req(identical(as.integer(cell$revision), sheet_revision()))
    tryCatch({
      edited(edit_cell(raw_data(), as.integer(cell$row), as.integer(cell$column), cell$value, input$decimal)); archive_view(FALSE)
      id <- working_table_id()
      if (identical(input$source, "workspace") && !is.null(id)) { p <- project(); p$tables[[id]]$data <- edited(); project(p) }
      if (identical(input$source, "upload")) {
        states <- upload_states(); key <- upload_key(); state <- states[[key]]
        state$data <- edited(); state$table_id <- id; states[[key]] <- state; upload_states(states)
        if (!is.null(id)) { p <- project(); p$tables[[id]]$data <- edited(); project(p) }
      }
    },
      error = function(e) { showNotification(conditionMessage(e), type = "error"); sheet_revision(sheet_revision() + 1L) })
  })
  output$wide_ui <- renderUI({
    d <- base_data(); numeric <- names(d)[vapply(d, is.numeric, logical(1))]
    previous <- isolate(input$wide_columns); id <- intersect(c("Subject", "subject", "id", "ID"), names(d))
    t <- if (identical(input$source, "workspace")) isolate(opened_table()) else NULL
    selected_id <- if (!is.null(t)) t$id_column else if (length(id)) id[1] else ""
    if (!is.null(t)) previous <- t$columns
    tagList(selectInput("wide_columns", "Measurement columns", names(d), selected = if (length(previous) && all(previous %in% numeric)) previous else setdiff(numeric, selected_id), multiple = TRUE),
      selectInput("wide_id", "Subject ID", c("Use row numbers" = "", names(d)), selected = selected_id))
  })
  # Analysis and restored projects need column selections on other tabs too.
  outputOptions(output, "wide_ui", suspendWhenHidden = FALSE)
  data <- reactive({
    d <- raw_data()
    if (identical(input$table_type, "nested") || (isTRUE(input$method %in% multiple_t_methods) && identical(input$t_test_mode, "endpoints"))) return(d)
    if (input$table_type != "xy" && input$layout == "wide") {
      req(input$wide_columns)
      tryCatch(wide_to_long(d, input$wide_columns, if (is.null(input$wide_id)) "" else input$wide_id), error = function(e) validate(need(FALSE, conditionMessage(e))))
    } else d
  })
  observeEvent(input$demo, {
    if (input$source != "demo") return()
    if (input$table_type == "nested") {
      updateRadioButtons(session, "nested_design", selected = if (input$demo == "nested_between") "between" else "repeated")
      return()
    }
    if (identical(input$demo, "multiple_tests")) {
      updateSelectInput(session, "family", selected = "multiple_t")
      updateSelectInput(session, "t_test_mode", selected = "endpoints")
      updateSelectInput(session, "layout", selected = "long")
      return()
    }
    if (input$table_type == "column") {
      updateSelectInput(session, "layout", selected = "wide")
      updateSelectInput(session, "family", selected = if (input$demo == "paired") "two" else "many")
      preferred <- preferences()
      updateRadioButtons(session, "matching", selected = if (!is.null(preferred$matching)) preferred$matching else if (input$demo %in% c("paired", "missing")) "matched" else "independent")
      updateRadioButtons(session, "repeated", selected = if (!is.null(preferred$repeated)) preferred$repeated else if (input$demo == "missing") "mixed" else "auto")
      return()
    }
    updateSelectInput(session, "layout", selected = "long")
    updateSelectInput(session, "method", selected = switch(input$demo, groups = "anova", paired = "paired_t", xy = "linear", dose = "dose", contingency = "fisher", survival = "survival", factorial = "two_way"))
  }, ignoreInit = TRUE)
  observeEvent(list(input$method, input$t_test_mode), {
    if (isTRUE(input$method %in% multiple_t_methods) && identical(input$t_test_mode, "endpoints")) updateSelectInput(session, "layout", selected = "long")
  })
  output$t_hypothesis_help <- renderText({
    req(input$method %in% t_test_methods)
    cfg <- current_config(); d <- data()
    groups <- if (identical(input$group, "")) "All" else levels(group_column(d, input$group))
    operator <- switch(cfg$alternative, greater = ">", less = "<", "!=")
    lognormal <- startsWith(cfg$method, "log_")
    target <- if (lognormal) cfg$null_ratio else cfg$null_value
    if (cfg$method %in% one_sample_methods) return(paste(if (lognormal) "Geometric mean of" else "Mean of", input_or("one_sample_group", groups[1]), operator, target))
    pairs <- if (cfg$method %in% multiple_t_methods && identical(cfg$t_test_mode, "columns")) {
      tryCatch(comparison_pairs(groups, cfg), error = function(e) list())
    } else list(if (isTRUE(cfg$reverse_comparisons)) rev(groups) else groups)
    if (!length(pairs)) return("Choose the pairs to test.")
    comparisons <- vapply(head(pairs, 6L), function(pair) paste(pair, collapse = if (lognormal) " / " else " - "), character(1))
    paste(if (lognormal) "Geometric-mean ratio:" else "Mean difference:", paste(paste(comparisons, operator, target), collapse = "; "), if (length(pairs) > 6L) paste("and", length(pairs) - 6L, "more comparisons") else "")
  })
  output$t_sample_group_ui <- renderUI({
    req(input$method %in% one_sample_methods)
    d <- data()
    groups <- if (identical(input$group, "")) "All" else levels(group_column(d, input$group))
    previous <- isolate(input$one_sample_group)
    selectInput("one_sample_group", "Sample to test", groups, selected = if (length(previous) == 1L && previous %in% groups) previous else groups[1])
  })
  outputOptions(output, "t_sample_group_ui", suspendWhenHidden = FALSE)
  output$mapping <- renderUI({
    if (identical(input$table_type, "nested")) return(nested_mapping_ui(raw_data(), input))
    d <- if (isTRUE(input$method %in% multiple_t_methods) && identical(input$t_test_mode, "endpoints")) raw_data() else if (input$table_type != "xy" && input$layout == "wide") {
      req(input$wide_columns)
      wide_to_long(base_data(), input$wide_columns, if (is.null(input$wide_id)) "" else input$wide_id)
    } else base_data()
    n <- names(d); nums <- n[vapply(d, is.numeric, logical(1))]; m <- input$method
    pick <- function(preferred, choices) { hits <- intersect(preferred, choices); if (length(hits)) hits[1] else if (length(choices)) choices[1] else character() }
    xy_method <- m %in% c("pearson", "spearman", "linear", "dose")
    x_selected <- pick(c(isolate(input$x), "concentration", "x", "X", "time"), nums)
    y_selected <- pick(c(isolate(input$y), if (m == "survival") "time" else c("response", "value", "y", "Y")), if (xy_method) setdiff(nums, x_selected) else nums)
    y <- selectInput("y", if (m == "survival") "Time" else "Response / Y", nums, selected = y_selected)
    group_selected <- pick(c(isolate(input$group), "group", "condition", "treatment"), n)
    g <- selectInput("group", if (m %in% c("chisq", "fisher")) "Row category" else "Group / condition", if (m == "summary" || m %in% one_sample_methods) c("All observations" = "", n) else n, selected = group_selected)
    if (m %in% t_test_methods) {
      paired <- m %in% c("paired_t", "log_paired_t") || (m %in% multiple_t_methods && identical(input$matching, "matched"))
      endpoint <- m %in% multiple_t_methods && identical(input$t_test_mode, "endpoints")
      return(tagList(y, g,
        if (m %in% one_sample_methods) uiOutput("t_sample_group_ui"),
        if (endpoint) selectInput("test_column", "Marker / time-point column", setdiff(n, c(y_selected, group_selected)), selected = pick(c(isolate(input$test_column), "endpoint", "marker", "timepoint", "time"), setdiff(n, c(y_selected, group_selected)))),
        if (paired) selectInput("id", "Matched subject ID", setdiff(n, c(y_selected, group_selected)), selected = pick(c(isolate(input$id), "subject", "Subject", "id"), setdiff(n, c(y_selected, group_selected))))))
    }
    if (m %in% c("chisq", "fisher")) return(tagList(g, selectInput("factor2", "Column category", n, selected = pick(c("outcome", "category"), setdiff(n, pick("group", n)))), selectInput("count", "Cell count", nums, selected = pick("count", nums))))
    if (xy_method) return(tagList(selectInput("x", if (m == "dose") "Concentration (original scale)" else "X", nums, selected = x_selected), y))
    tagList(y, g,
      if (m %in% c("paired_t", "paired_wilcox", "friedman", "rm_anova", "mixed", "log_rm_anova", "log_mixed", "rm_auto", "log_rm_auto")) selectInput("id", "Matched subject ID", n, selected = pick(c("subject", "Subject", "id", "sample"), n)),
      if (m == "two_way") selectInput("factor2", "Second factor", n, selected = pick(c("sex", "factor2"), setdiff(n, pick("group", n)))),
      if (m == "survival") selectInput("event", "Event indicator (0 / 1)", nums, selected = pick("event", nums)))
  })
  # Populate column mappings even when Analysis options is hidden so Analyze
  # can run from the data sheet without requiring a visit to that tab first.
  outputOptions(output, "mapping", suspendWhenHidden = FALSE)
  config <- reactive(list(method = input$method, y = input$y, group = input$group, x = input$x, id = input$id, factor2 = input$factor2, event = input$event, count = input$count, confidence = as.numeric(input_or("confidence", .95)), adjustment = input_or("adjustment", "sidak"), direction = input$direction, sphericity = input_or("sphericity", "gg"), alpha = input_or("alpha", .05), comparisons = input_or("comparisons", "all"), control = input$control, rout_q = as.numeric(input_or("rout_q", .01))))
  # Track only fields used by the selected method, so hidden controls cannot change results.
  current_config <- reactive({
    z <- config(); m <- z$method; req(length(m) == 1, nzchar(m))
    if (m %in% nested_methods) return(nested_config(input, raw_data()))
    if (m %in% t_test_methods) {
      z <- list(method = m, y = input$y, group = input$group, confidence = as.numeric(input_or("confidence", .95)),
        alpha = input_or("alpha", .05), alternative = input_or("t_alternative", "two.sided"),
        report_descriptive = isTRUE(input_or("report_descriptive", TRUE)), report_adjusted_p = isTRUE(input_or("report_adjusted_p", TRUE)))
      if (startsWith(m, "log_")) z$null_ratio <- input_or("t_ratio", 1) else z$null_value <- input_or("t_mu", 0)
      kind <- t_test_kind(m, list(t_test_kind = if (identical(input$matching, "matched")) "paired_t" else if (identical(input$variance, "equal")) "student" else "welch"))
      if (kind == "paired_t") z$id <- input$id
      if (m %in% one_sample_methods && !identical(input$group, "")) z$one_sample_group <- input$one_sample_group
      if (m %in% multiple_t_methods) {
        z$t_test_kind <- kind; z$t_test_mode <- input_or("t_test_mode", "endpoints")
        if (z$t_test_mode == "endpoints") z$test_column <- input$test_column
        z$comparisons <- if (z$t_test_mode == "endpoints") "all" else input_or("comparisons", "all")
        if (z$comparisons == "control") z$control <- input$control
        if (z$comparisons == "selected") z$selected_pairs <- as.character(input$selected_pairs)
        z$reverse_comparisons <- isTRUE(input$reverse_comparisons)
        z$adjustment <- if (identical(input$error_control, "none")) "none" else input_or("adjustment", "holm")
        if (identical(input$error_control, "fdr") && !z$adjustment %in% c("BKY", "BH", "BY")) z$adjustment <- "BH"
        if (z$adjustment %in% c("BKY", "BH", "BY")) z$fdr_q <- input_or("fdr_q", .05)
        z$diagnostic_graphs <- if (isTRUE(input_or("graph_comparison", TRUE))) "comparison" else character()
      }
      return(z)
    }
    keys <- c("method", "confidence", "alpha")
    if (m %in% c("chisq", "fisher")) keys <- c(keys, "group", "factor2", "count")
    else if (m %in% c("pearson", "spearman", "linear", "dose")) keys <- c(keys, "x", "y", if (m == "dose") "direction")
    else keys <- c(keys, "y", "group", if (m %in% c("paired_t", "paired_wilcox", "friedman", repeated_methods)) "id", if (m == "two_way") "factor2", if (m == "survival") "event")
    if (m %in% one_way_methods) keys <- c(keys, "adjustment", "comparisons", if (identical(z$comparisons, "control")) "control")
    if (m %in% repeated_methods) keys <- c(keys, "sphericity")
    if (m == "outliers") keys <- c(keys, "rout_q")
    z <- z[keys]
    if (m %in% one_way_methods && identical(input$error_control, "none")) z$adjustment <- "none"
    if (m %in% one_way_methods && identical(input$error_control, "fdr") && !z$adjustment %in% c("BKY", "BH", "BY")) z$adjustment <- "BKY"
    modifyList(z, advanced_config())
  })
  is_current <- reactive({ s <- saved(); !is.null(s) && (archive_view() || (identical(s$table_id, working_table_id()) && identical(s$config, current_config()) && identical(s$input_data, data()))) })
  # Compare live inputs even while viewing an archived result.
  analysis_unchanged <- reactive({
    s <- saved()
    !is.null(s) && identical(s$table_id, working_table_id()) && identical(s$config, current_config()) && identical(s$input_data, data())
  })
  observe({
    unchanged <- tryCatch(analysis_unchanged(), error = function(e) FALSE)
    updateActionButton(session, "run", disabled = unchanged)
  })
  result <- reactive({ req(saved()); validate(need(is_current(), "Data or analysis settings changed. Run the analysis again to view or export results.")); saved()$result })
  observeEvent(input$run, {
    # Analyze is the calculation boundary. Save inputs and settings with the result
    # so later worksheet changes can be detected instead of silently changing it.
    # Reject queued repeat clicks before clearing the result or creating saved sheets.
    if (isTRUE(tryCatch(analysis_unchanged(), error = function(e) FALSE))) return()
    saved(NULL); analysis_error(NULL)
    progress <- NULL
    update_progress <- function(value = NULL, detail = NULL) {
      if (!is.null(progress)) progress$set(value = value, detail = detail)
    }
    tryCatch({
      # req() silently stops reactive outputs, but Analyze needs a visible error.
      if (identical(input$source, "upload")) {
        assert(!is.null(input$file) && nrow(input$file) > 0L,
          "Choose a CSV, TSV or Excel file using Browse before clicking Run analysis.")
        if (file_is_excel()) assert(isTruthy(input$sheet), "Choose an Excel worksheet before clicking Run analysis.")
      }
      if (identical(input$source, "paste")) assert(!is.null(pasted()),
        "Paste a table and click Use pasted table before clicking Run analysis.")
      if (identical(input$source, "workspace")) assert(!is.null(opened_table()),
        "Open a saved project table or choose another data source before clicking Run analysis.")
      if (!input$table_type %in% c("xy", "nested") && identical(input$layout, "wide")) assert(isTruthy(input$wide_columns),
        "Select the measurement columns to analyze.")
      d <- data(); cfg <- current_config()
      if (input$table_type == "column" && input$family == "many") assert(nlevels(group_column(d, cfg$group)) >= 3, "Select at least three group columns, or choose the two-group comparison.")
      if (cfg$method %in% nested_methods) {
        progress <- shiny::Progress$new(session, style = "notification")
        progress$set(message = "Running linear mixed analysis", value = 0,
          detail = "Starting analysis; progress follows completed stages")
      }
      res <- run_analysis(d, cfg, progress = update_progress)
      update_progress(.95, "Saving analysis and preparing graph sheets")
      table_id <- save_current_table()
      p <- project(); name <- paste(names(analysis_methods)[match(res$method, analysis_methods)], "of", p$tables[[table_id]]$name)
      added <- project_add(p, "analyses", list(name = name, table_id = table_id, worksheet = p$tables[[table_id]], result = res, input_data = d, config = cfg, time = as.character(Sys.time())))
      p <- added$project; analysis_id <- added$id
      kinds <- if (isTRUE(res$nested)) names(res$graphs) else unique(c("raw", "normal_qq", "lognormal_qq", if (is.null(cfg$diagnostic_graphs)) "comparison" else cfg$diagnostic_graphs))
      for (kind in kinds) {
        possible <- tryCatch({ make_graph(res, kind); TRUE }, error = function(e) FALSE)
        if (possible) p <- project_add(p, "graphs", list(name = paste(names(graph_kinds)[match(kind, graph_kinds)], "of", name), analysis_id = analysis_id, kind = kind, options = list()))$project
      }
      project(p); archive_view(FALSE); saved(p$analyses[[analysis_id]])
      removeNotification("analysis-error")
      select_sheet("analyses", analysis_id)
      states <- upload_states()
      for (key in names(states)) {
        if (identical(states[[key]]$table_id, table_id)) states[[key]]$analysis_id <- analysis_id
      }
      upload_states(states)
      updateTextInput(session, "analysis_name", value = name)
      updateSelectInput(session, "graph_kind", selected = "raw")
      updateTabsetPanel(session, "tabs", selected = "results")
      update_progress(1, "Analysis complete")
    }, error = function(e) {
      message <- trimws(conditionMessage(e))
      if (!nzchar(message)) message <- "Select a data table and complete the analysis settings before clicking Run analysis."
      analysis_error(message)
      showNotification(message, id = "analysis-error", type = "error", duration = NULL)
    }, finally = {
      if (!is.null(progress)) progress$close()
    })
  })
  output$method_help <- renderText(method_note(input$method))
  output$sheet_title <- renderText(if (identical(input$source, "workspace")) opened_table()$name else switch(input$table_type, column = "Column data table", xy = "XY data table", nested = "Linear mixed analysis data table", "Data table"))
  output$data_info <- renderText(paste(nrow(raw_data()), "rows |", ncol(raw_data()), "columns"))
  output$table_help <- renderUI({
    if (identical(input$table_type, "nested")) return(p("Repeated models: one row per Mouse / Slice / Cell and one numeric column per condition. Between-group models: Mouse, Slice, Cell, Condition and a numeric measurement column. ID combinations identify cells within slices and mice; IDs must be present even when a measurement is blank. Blank measurements are omitted individually."))
    if (input$table_type == "xy") p("Enter X in one numeric column and Y in another. Choose simple linear regression, Pearson correlation, or Spearman correlation. If you have several Y columns, select the one to analyze.")
    else p("Each measurement column is a group or condition. For matched data, values in the same row belong to the same subject. For independent data, rows do not imply matching. Subject IDs are optional in wide tables. Missing values are audited; complete-case matched tests exclude incomplete subjects, while mixed-effects analysis retains their observed values.")
  })
  output$worksheet <- renderUI({
    base_data(); revision <- sheet_revision()
    d <- head(isolate(raw_data()), 100)
    tags$table(class = "worksheet", tags$thead(tags$tr(tags$th(class = "row-number", ""), lapply(seq_along(d), function(j) tags$th(span(class = "column-letter", LETTERS[(j - 1) %% 26 + 1]), names(d)[j])))),
      tags$tbody(lapply(seq_len(nrow(d)), function(i) tags$tr(tags$th(class = "row-number", i), lapply(seq_along(d), function(j) tags$td(tags$input(
        class = "worksheet-cell", type = "text", value = if (is.na(d[i, j])) "" else as.character(d[i, j]),
        `data-row` = i, `data-column` = j, `data-revision` = revision, `aria-label` = paste(names(d)[j], "row", i))))))))
  })
  output$preview <- renderTable(head(data(), 100), striped = TRUE, bordered = FALSE, spacing = "s")
  output$status <- renderUI({
    if (is.null(saved()) && !is.null(analysis_error())) return(div(class = "status-bar stale", role = "alert", analysis_error()))
    if (is.null(saved())) return(div(class = "status-bar", span(class = "status-dot"), "Data worksheet", span(class = "status-detail", "Choose Analysis to run calculations")))
    if (!is_current()) return(div(class = "status-bar stale", "Settings changed \u2014 run the analysis again."))
    a <- result()$tables$data_audit
    div(class = "status-bar", span(class = "status-dot"), if (archive_view()) "Saved analysis snapshot" else "Analysis complete", span(class = "status-detail", paste(a$used_rows, if (isTRUE(result()$nested)) "measurements used /" else "rows used /", a$excluded_rows, "excluded")))
  })
  output$results_header <- renderUI({ r <- result(); tagList(span(class = "eyebrow", "ANALYSIS RESULTS"), h3(names(analysis_methods)[match(r$method, analysis_methods)])) })
  output$notes <- renderUI({ r <- result(); tagList(div(class = "method-note", lapply(r$notes, tags$p)), if (length(r$warnings)) div(class = "warning-note", h4("Review before interpreting"), lapply(r$warnings, tags$p))) })
  output$cleaned_action <- renderUI({ if (result()$method == "outliers") actionButton("save_cleaned", "Create cleaned data table") })
  display_sheets <- reactive(result_sheets(result(), worksheet = saved()$worksheet))
  output$results_tables <- renderUI({
    sheets <- display_sheets(); additional <- attr(sheets, "additional")
    panels <- lapply(names(sheets), function(n) tabPanel(sheets[[n]]$title, value = n, uiOutput(paste0("sheet_", n))))
    if (length(additional)) panels <- c(panels, list(tabPanel("Additional results", value = "additional", do.call(tabsetPanel, c(list(id = "additional_result_sheets"),
      lapply(names(additional), function(n) tabPanel(gsub("_", " ", tools::toTitleCase(n)), value = n, div(class = "table-scroll result-sheet", tableOutput(paste0("table_", n))))))))))
    do.call(tabsetPanel, c(list(id = "result_sheets", selected = if (length(sheets)) names(sheets)[1] else "additional"), panels))
  })
  observeEvent(saved(), {
    req(saved())
    sheets <- result_sheets(saved()$result, worksheet = saved()$worksheet)
    choices <- result_export_choices(saved()$result, worksheet = saved()$worksheet)
    updateSelectInput(session, "export_table", choices = choices, selected = unname(choices[1]))
    for (name in names(sheets)) local({
      n <- name
      output[[paste0("sheet_", n)]] <- renderUI(result_sheet_ui(display_sheets()[[n]], display_options(), names(analysis_methods)[match(result()$method, analysis_methods)]))
    })
    for (name in names(saved()$result$tables)) local({
      n <- name
      output[[paste0("table_", n)]] <- renderTable({
        table <- display_tables()[[n]]; req(table)
        opts <- display_options()
        format_result_table(table, opts$output_digits, opts$p_style, opts$p_digits)
      }, digits = 6, striped = TRUE, spacing = "s", na = "\u2014")
    })
  })
  plot_options <- reactive({
    z <- list(style = input$plot_style, error = input$error_style, title = input$plot_title, x_label = input$x_label, y_label = input$y_label,
      point_size = input$point_size, symbol = if (is.null(input$symbol)) NULL else as.numeric(input$symbol), palette_name = input$palette_name,
      connect = input$connect, varied_symbols = input$varied_symbols, legend = input$legend, x_angle = if (is.null(input$x_angle)) NULL else as.numeric(input$x_angle), show_stars = input$show_stars)
    z[!vapply(z, is.null, logical(1))]
  })
  figure <- reactive(make_graph(result(), if (is.null(input$graph_kind)) "raw" else input$graph_kind, plot_options()))
  output$figure <- renderPlot(draw_nested_graph(figure()), res = 110)
  output$diagnostic_help <- renderUI({ r <- result(); if (is.null(diagnostic_plot(r))) p("Residual diagnostics are available for ANOVA, mixed-effects, linear regression, and dose-response fits.") else p("Look for changing spread, curvature, and outliers. For mixed-effects models, these are conditional residuals. This plot does not establish that model assumptions hold.") })
  output$diagnostics <- renderPlot({ p <- diagnostic_plot(result()); if (!is.null(p)) draw_nested_graph(p) }, res = 100)
  output$normality_figure <- renderPlot({ p <- normality_plot(result()); if (!is.null(p)) draw_nested_graph(p) }, res = 100)
  output$settings <- renderPrint({ result(); print(saved()$config) })
  output$download_data <- downloadHandler(filename = "input_data.csv", content = function(file) write.csv(data(), file, row.names = FALSE))
  output$download_table <- downloadHandler(filename = function() paste0(input$export_table, ".csv"), content = function(file) { tables <- result_export_tables(result(), worksheet = saved()$worksheet, options = display_options()); req(input$export_table %in% names(tables)); write.csv(tables[[input$export_table]], file, row.names = FALSE) })
  output$download_report <- downloadHandler(filename = "analysis_report.txt", content = function(file) writeLines(c(analysis_report(result(), display_options(), saved()$worksheet), "", "SETTINGS", capture.output(dput(saved()$config))), file))
  output$download_script <- downloadHandler(filename = "reproduce_analysis.R", content = function(file) { result(); writeLines(reproduction_script(saved()$input_data, modifyList(saved()$config, display_options()), plot_options(), if (is.null(input$graph_kind)) "raw" else input$graph_kind, saved()$worksheet), file) })
  output$download_session <- downloadHandler(filename = "session_info.txt", content = function(file) writeLines(capture.output(sessionInfo()), file))
  output$download_png <- downloadHandler(filename = "figure.png", content = function(file) save_result_figure(result(), input$graph_kind, figure(), file, "png"))
  output$download_svg <- downloadHandler(filename = "figure.svg", content = function(file) save_result_figure(result(), input$graph_kind, figure(), file, "svg"))
}

shinyApp(ui, server)
