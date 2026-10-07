if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
required <- c("shiny", "ggplot2", "readxl", "survival", "nlme", "jsonlite")
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
  saved <- reactiveVal(NULL)
  analysis_error <- reactiveVal(NULL)
  edited <- reactiveVal(NULL)
  pasted <- reactiveVal(NULL)
  sheet_revision <- reactiveVal(0L)
  source("R/workspace.inc", local = TRUE)
  observeEvent(input$import_data, {
    updateSelectInput(session, "source", selected = "upload")
    updateTabsetPanel(session, "tabs", selected = "data")
  })
  output$distribution_ui <- renderUI({
    choices <- c("Approximately Gaussian (parametric)" = "normal", "Rank-based (nonparametric)" = "rank")
    if (identical(input$matching, "matched") && identical(input$family, "many")) choices <- c(choices, "Lognormal (analyze log responses)" = "lognormal")
    selected <- isolate(input$distribution); if (is.null(selected) || !selected %in% choices) selected <- "normal"
    radioButtons("distribution", "Distribution assumption", choices, selected = selected)
  })
  output$adjustment_ui <- renderUI({
    choices <- c("Sidak (adjusted p-values and intervals)" = "sidak", "Holm (p-values; unadjusted intervals)" = "holm", "Bonferroni (p-values and intervals)" = "bonferroni", "Benjamini-Hochberg (p-values only)" = "BH", "None" = "none")
    if (identical(input$method, "anova")) choices <- c(choices, "Tukey (all pairs only)" = "tukey")
    selected <- isolate(input$adjustment); if (is.null(selected) || !selected %in% choices) selected <- "sidak"
    selectInput("adjustment", "Multiple-comparison correction", choices, selected = selected)
  })
  output$control_ui <- renderUI({
    d <- data(); req(input$group %in% names(d)); selectInput("control", "Control group", levels(group_column(d, input$group)))
  })
  output$demo_ui <- renderUI({
    choices <- switch(input$table_type,
      column = c("Three independent groups" = "groups", "Matched before / after" = "paired", "Repeated measures with missing values" = "missing"),
      xy = c("XY measurements" = "xy"),
      other = c("Three groups" = "groups", "Matched before / after" = "paired", "XY measurements" = "xy", "Dose-response" = "dose", "Category counts" = "contingency", "Time to event" = "survival", "Two-factor experiment" = "factorial"))
    selectInput("demo", "Example table", choices)
  })
  observeEvent(input$table_type, {
    layout <- if (identical(input$source, "workspace") && !is.null(opened_table())) opened_table()$layout else if (input$table_type == "column") "wide" else "long"
    updateSelectInput(session, "layout", selected = layout)
    updateTabsetPanel(session, "tabs", selected = "data")
  })
  observeEvent(list(input$table_type, input$family, input$matching, input$distribution, input$variance, input$repeated), {
    req(input$table_type)
    choices <- design_methods(input$table_type, if (is.null(input$family)) "many" else input$family, identical(input$matching, "matched"), if (is.null(input$distribution)) "normal" else input$distribution, if (is.null(input$variance)) "unequal" else input$variance, if (is.null(input$repeated)) "complete" else input$repeated)
    selected <- isolate(input$method)
    if (is.null(selected) || !selected %in% choices) selected <- unname(choices[1])
    updateSelectInput(session, "method", choices = choices, selected = selected)
  })
  observeEvent(input$apply_paste, {
    tryCatch({ pasted(paste_table(input$paste_text, input$decimal)); showNotification("Pasted table loaded.", type = "message") },
      error = function(e) showNotification(conditionMessage(e), type = "error"))
  })
  file_is_excel <- reactive({ req(input$file); tolower(tools::file_ext(input$file$name)) %in% c("xlsx", "xls") })
  output$sheet_ui <- renderUI({
    req(input$source == "upload", input$file, file_is_excel())
    sheets <- tryCatch(readxl::excel_sheets(input$file$datapath), error = function(e) character())
    selectInput("sheet", "Excel worksheet", sheets)
  })
  base_data <- reactive({
    if (input$source == "workspace") { req(opened_table()); return(opened_table()$data) }
    if (input$source == "demo") {
      req(input$demo)
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
    req(input$file)
    tryCatch({
      if (file_is_excel()) {
        req(input$sheet)
        validate_import(readxl::read_excel(input$file$datapath, sheet = input$sheet, .name_repair = "minimal"))
      } else {
        sep <- if (tolower(tools::file_ext(input$file$name)) %in% c("tsv", "txt")) "\t" else input$separator
        validate_import(read.table(input$file$datapath, header = TRUE, sep = sep, dec = input$decimal, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"), comment.char = "", quote = "\""))
      }
    }, error = function(e) validate(need(FALSE, conditionMessage(e))))
  })
  observeEvent(base_data(), {
    edited(NULL); sheet_revision(isolate(sheet_revision()) + 1L)
    if (input$source != "workspace") { active_table_id(NULL); working_table_id(NULL) }
  }, priority = 100)
  raw_data <- reactive({ d <- base_data(); if (is.null(edited())) d else edited() })
  observeEvent(input$cell_edit, {
    cell <- input$cell_edit
    req(identical(as.integer(cell$revision), sheet_revision()))
    tryCatch({
      edited(edit_cell(raw_data(), as.integer(cell$row), as.integer(cell$column), cell$value, input$decimal)); archive_view(FALSE)
      id <- working_table_id()
      if (identical(input$source, "workspace") && !is.null(id)) { p <- project(); p$tables[[id]]$data <- edited(); project(p) }
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
  data <- reactive({
    d <- raw_data()
    if (input$table_type != "xy" && input$layout == "wide") {
      req(input$wide_columns)
      tryCatch(wide_to_long(d, input$wide_columns, if (is.null(input$wide_id)) "" else input$wide_id), error = function(e) validate(need(FALSE, conditionMessage(e))))
    } else d
  })
  observeEvent(input$demo, {
    if (input$source != "demo") return()
    if (input$table_type == "column") {
      updateSelectInput(session, "layout", selected = "wide")
      updateSelectInput(session, "family", selected = if (input$demo == "paired") "two" else "many")
      updateRadioButtons(session, "matching", selected = if (input$demo %in% c("paired", "missing")) "matched" else "independent")
      updateSelectInput(session, "repeated", selected = if (input$demo == "missing") "mixed" else "complete")
      return()
    }
    updateSelectInput(session, "layout", selected = "long")
    updateSelectInput(session, "method", selected = switch(input$demo, groups = "anova", paired = "paired_t", xy = "linear", dose = "dose", contingency = "fisher", survival = "survival", factorial = "two_way"))
  }, ignoreInit = TRUE)
  output$mapping <- renderUI({
    d <- if (input$table_type != "xy" && input$layout == "wide") {
      req(input$wide_columns)
      wide_to_long(base_data(), input$wide_columns, if (is.null(input$wide_id)) "" else input$wide_id)
    } else base_data()
    n <- names(d); nums <- n[vapply(d, is.numeric, logical(1))]; m <- input$method
    pick <- function(preferred, choices) { hits <- intersect(preferred, choices); if (length(hits)) hits[1] else if (length(choices)) choices[1] else character() }
    xy_method <- m %in% c("pearson", "spearman", "linear", "dose")
    x_selected <- pick(c(isolate(input$x), "concentration", "x", "X", "time"), nums)
    y_selected <- pick(c(isolate(input$y), if (m == "survival") "time" else c("response", "value", "y", "Y")), if (xy_method) setdiff(nums, x_selected) else nums)
    y <- selectInput("y", if (m == "survival") "Time" else "Response / Y", nums, selected = y_selected)
    g <- selectInput("group", if (m %in% c("chisq", "fisher")) "Row category" else "Group / condition", if (m == "summary") c("All observations" = "", n) else n, selected = pick(c("group", "condition", "treatment"), n))
    if (m %in% c("chisq", "fisher")) return(tagList(g, selectInput("factor2", "Column category", n, selected = pick(c("outcome", "category"), setdiff(n, pick("group", n)))), selectInput("count", "Cell count", nums, selected = pick("count", nums))))
    if (xy_method) return(tagList(selectInput("x", if (m == "dose") "Concentration (original scale)" else "X", nums, selected = x_selected), y))
    tagList(y, g,
      if (m %in% c("paired_t", "paired_wilcox", "friedman", "rm_anova", "mixed", "log_rm_anova", "log_mixed")) selectInput("id", "Matched subject ID", n, selected = pick(c("subject", "Subject", "id", "sample"), n)),
      if (m == "two_way") selectInput("factor2", "Second factor", n, selected = pick(c("sex", "factor2"), setdiff(n, pick("group", n)))),
      if (m == "survival") selectInput("event", "Event indicator (0 / 1)", nums, selected = pick("event", nums)))
  })
  config <- reactive(list(method = input$method, y = input$y, group = input$group, x = input$x, id = input$id, factor2 = input$factor2, event = input$event, count = input$count, confidence = as.numeric(input$confidence), adjustment = input$adjustment, direction = input$direction, sphericity = input$sphericity, alpha = if (is.null(input$alpha)) .05 else input$alpha, comparisons = if (is.null(input$comparisons)) "all" else input$comparisons, control = input$control, rout_q = if (is.null(input$rout_q)) .01 else as.numeric(input$rout_q)))
  # Track only fields used by the selected method, so hidden controls cannot change results.
  current_config <- reactive({
    z <- config(); m <- z$method
    keys <- c("method", "confidence", "alpha")
    if (m %in% c("chisq", "fisher")) keys <- c(keys, "group", "factor2", "count")
    else if (m %in% c("pearson", "spearman", "linear", "dose")) keys <- c(keys, "x", "y", if (m == "dose") "direction")
    else keys <- c(keys, "y", "group", if (m %in% c("paired_t", "paired_wilcox", "friedman", "rm_anova", "mixed", "log_rm_anova", "log_mixed")) "id", if (m == "two_way") "factor2", if (m == "survival") "event")
    if (m %in% c("anova", "welch_anova", "kruskal", "rm_anova", "mixed", "log_rm_anova", "log_mixed")) keys <- c(keys, "adjustment", "comparisons", if (identical(z$comparisons, "control")) "control")
    if (m %in% c("rm_anova", "mixed", "log_rm_anova", "log_mixed")) keys <- c(keys, "sphericity")
    if (m == "outliers") keys <- c(keys, "rout_q")
    z[keys]
  })
  is_current <- reactive({ s <- saved(); !is.null(s) && (archive_view() || (identical(s$config, current_config()) && identical(s$input_data, data()))) })
  # Compare live inputs even while viewing an archived result.
  analysis_unchanged <- reactive({
    s <- saved()
    !is.null(s) && identical(s$config, current_config()) && identical(s$input_data, data())
  })
  observe({
    unchanged <- tryCatch(analysis_unchanged(), error = function(e) FALSE)
    updateActionButton(session, "run", disabled = unchanged)
  })
  result <- reactive({ req(saved()); validate(need(is_current(), "Data or analysis settings changed. Run the analysis again to view or export results.")); saved()$result })
  observeEvent(input$run, {
    # Reject queued repeat clicks before clearing the result or creating saved sheets.
    if (isTRUE(tryCatch(analysis_unchanged(), error = function(e) FALSE))) return()
    saved(NULL); analysis_error(NULL)
    tryCatch({
      d <- data(); cfg <- current_config()
      if (input$table_type == "column" && input$family == "many") assert(nlevels(group_column(d, cfg$group)) >= 3, "Select at least three group columns, or choose the two-group comparison.")
      res <- run_analysis(d, cfg)
      table_id <- save_current_table()
      p <- project(); name <- paste(names(analysis_methods)[match(res$method, analysis_methods)], "of", p$tables[[table_id]]$name)
      added <- project_add(p, "analyses", list(name = name, table_id = table_id, worksheet = p$tables[[table_id]], result = res, input_data = d, config = cfg, time = as.character(Sys.time())))
      p <- added$project; analysis_id <- added$id
      for (kind in c("raw", "normal_qq", "lognormal_qq", "comparison")) {
        possible <- tryCatch({ make_graph(res, kind); TRUE }, error = function(e) FALSE)
        if (possible) p <- project_add(p, "graphs", list(name = paste(names(graph_kinds)[match(kind, graph_kinds)], "of", name), analysis_id = analysis_id, kind = kind, options = list()))$project
      }
      project(p); archive_view(FALSE); saved(p$analyses[[analysis_id]])
      updateTextInput(session, "analysis_name", value = name)
      updateSelectInput(session, "graph_kind", selected = "raw")
      updateSelectInput(session, "export_table", choices = names(res$tables), selected = if ("test" %in% names(res$tables)) "test" else names(res$tables)[1])
      updateTabsetPanel(session, "tabs", selected = "results")
    }, error = function(e) { analysis_error(conditionMessage(e)); showNotification(conditionMessage(e), type = "error", duration = NULL) })
  })
  output$method_help <- renderText(method_note(input$method))
  output$sheet_title <- renderText(if (identical(input$source, "workspace")) opened_table()$name else switch(input$table_type, column = "Column data table", xy = "XY data table", "Data table"))
  output$data_info <- renderText(paste(nrow(raw_data()), "rows |", ncol(raw_data()), "columns"))
  output$table_help <- renderUI({
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
  for (pair in list(c("show_data", "data"), c("show_results", "results"), c("show_graph", "figure"), c("nav_data", "data"), c("nav_results", "results"), c("nav_graph", "figure"))) local({
    button <- pair[1]; tab <- pair[2]
    observeEvent(input[[button]], updateTabsetPanel(session, "tabs", selected = tab))
  })
  output$preview <- renderTable(head(data(), 100), striped = TRUE, bordered = FALSE, spacing = "s")
  output$status <- renderUI({
    if (is.null(saved()) && !is.null(analysis_error())) return(div(class = "status-bar stale", analysis_error()))
    if (is.null(saved())) return(div(class = "status-bar", span(class = "status-dot"), "Data worksheet", span(class = "status-detail", "Enter data and click Analyze")))
    if (!is_current()) return(div(class = "status-bar stale", "Settings changed \u2014 run the analysis again."))
    a <- result()$tables$data_audit
    div(class = "status-bar", span(class = "status-dot"), if (archive_view()) "Saved analysis snapshot" else "Analysis complete", span(class = "status-detail", paste(a$used_rows, "rows used /", a$excluded_rows, "excluded")))
  })
  output$results_header <- renderUI({ r <- result(); tagList(span(class = "eyebrow", "ANALYSIS RESULTS"), h3(names(analysis_methods)[match(r$method, analysis_methods)])) })
  output$notes <- renderUI({ r <- result(); tagList(div(class = "method-note", lapply(r$notes, tags$p)), if (length(r$warnings)) div(class = "warning-note", h4("Review before interpreting"), lapply(r$warnings, tags$p))) })
  output$cleaned_action <- renderUI({ if (result()$method == "outliers") actionButton("save_cleaned", "Create cleaned data table") })
  output$results_tables <- renderUI({
    r <- result()
    do.call(tabsetPanel, c(list(id = "result_sheets", selected = if ("outlier_summary" %in% names(r$tables)) "outlier_summary" else if ("test" %in% names(r$tables)) "test" else if ("distribution_tests" %in% names(r$tables)) "distribution_tests" else names(r$tables)[1]),
      lapply(names(r$tables), function(n) tabPanel(gsub("_", " ", tools::toTitleCase(n)), value = n, div(class = "table-scroll result-sheet", tableOutput(paste0("table_", n)))))))
  })
  observeEvent(saved(), {
    req(saved())
    for (name in names(saved()$result$tables)) local({
      n <- name
      output[[paste0("table_", n)]] <- renderTable({
        table <- result()$tables[[n]]
        for (col in intersect(names(table), c("p_value", "adjusted_p_value", "p adj", "Pr(>F)", "Pr(>|t|)"))) {
          table[[col]] <- ifelse(is.na(table[[col]]), NA_character_, format.pval(table[[col]], digits = 4, eps = 1e-15))
        }
        table
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
  output$figure <- renderPlot(print(figure()), res = 110)
  output$diagnostic_help <- renderUI({ r <- result(); if (is.null(diagnostic_plot(r))) p("Residual diagnostics are available for ANOVA, mixed-effects, linear regression, and dose-response fits.") else p("Look for changing spread, curvature, and outliers. For mixed-effects models, these are conditional residuals. This plot does not establish that model assumptions hold.") })
  output$diagnostics <- renderPlot({ p <- diagnostic_plot(result()); if (!is.null(p)) print(p) }, res = 100)
  output$normality_figure <- renderPlot({ p <- normality_plot(result()); if (!is.null(p)) print(p) }, res = 100)
  output$settings <- renderPrint({ result(); print(saved()$config) })
  output$download_data <- downloadHandler(filename = "input_data.csv", content = function(file) write.csv(data(), file, row.names = FALSE))
  output$download_table <- downloadHandler(filename = function() paste0(input$export_table, ".csv"), content = function(file) { r <- result(); req(input$export_table %in% names(r$tables)); write.csv(r$tables[[input$export_table]], file, row.names = FALSE) })
  output$download_report <- downloadHandler(filename = "analysis_report.txt", content = function(file) writeLines(c(analysis_report(result()), "", "SETTINGS", capture.output(dput(saved()$config))), file))
  output$download_script <- downloadHandler(filename = "reproduce_analysis.R", content = function(file) { result(); writeLines(reproduction_script(saved()$input_data, saved()$config, plot_options(), if (is.null(input$graph_kind)) "raw" else input$graph_kind), file) })
  output$download_session <- downloadHandler(filename = "session_info.txt", content = function(file) writeLines(capture.output(sessionInfo()), file))
  output$download_png <- downloadHandler(filename = "figure.png", content = function(file) ggplot2::ggsave(file, figure(), width = 7, height = 5, dpi = 300, bg = "white"))
  output$download_pdf <- downloadHandler(filename = "figure.pdf", content = function(file) ggplot2::ggsave(file, figure(), width = 7, height = 5, device = "pdf", bg = "white"))
}

shinyApp(ui, server)
