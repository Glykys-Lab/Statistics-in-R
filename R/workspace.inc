# Sourced inside the Shiny server: project navigation and immutable analysis snapshots.
# Own the session's project records and selected-sheet state. Loading a worksheet
# restores editable inputs; opening an archived result preserves its saved analysis.
# Transformations create derived tables and provenance instead of replacing sources.
project <- reactiveVal(new_project())
opened_table <- reactiveVal(NULL)
active_table_id <- reactiveVal(NULL)
selected_sheet <- reactiveVal(NULL)
select_sheet <- function(kind, id) {
  selected_sheet(if (is.null(id)) NULL else list(kind = kind, id = id))
}
working_table_id <- reactiveVal(NULL)
archive_view <- reactiveVal(FALSE)

notify_error <- function(e) showNotification(conditionMessage(e), type = "error", duration = NULL)
guard <- function(code) tryCatch(force(code), error = notify_error)

load_table <- function(id) {
  # Defer selector updates until Shiny has mounted the controls for the new table.
  t <- project()$tables[[id]]; req(t)
  active_table_id(id); working_table_id(id); opened_table(t); edited(NULL); saved(NULL); archive_view(FALSE)
  select_sheet("tables", id)
  updateRadioButtons(session, "table_type", selected = t$type)
  updateSelectInput(session, "source", selected = "workspace")
  updateSelectInput(session, "layout", selected = t$layout)
  updateTextInput(session, "table_name", value = t$name)
  session$onFlushed(function() {
    updateSelectInput(session, "wide_columns", selected = t$columns)
    updateSelectInput(session, "wide_id", selected = t$id_column)
    restore_nested_config(t$nested_config)
  }, once = TRUE)
  updateTabsetPanel(session, "tabs", selected = "data")
}

save_current_table <- function(copy = FALSE, name = input$table_name, d = raw_data(), provenance = NULL) {
  assert(!is.null(d) && nrow(d) > 0, "Enter or load a table first.")
  p <- isolate(project()); id <- isolate(working_table_id())
  record <- list(name = if (is.null(name) || !nzchar(trimws(name))) "Untitled table" else trimws(name),
    type = input$table_type, layout = if (input$table_type == "nested") if (identical(input$nested_design, "between")) "long" else "wide" else if (input$table_type == "xy") "long" else input$layout,
    columns = names(d)[names(d) %in% input$wide_columns], id_column = if (is.null(input$wide_id)) "" else input$wide_id,
    data = d, provenance = provenance, nested_config = if (identical(input$table_type, "nested")) isolate(current_config()) else NULL)
  if (!copy && !is.null(id) && id %in% names(p$tables)) {
    record$id <- id; p$tables[[id]] <- record
  } else { added <- project_add(p, "tables", record); p <- added$project; id <- added$id }
  project(p); active_table_id(id); working_table_id(id)
  select_sheet("tables", id)
  if (!copy && identical(input$source, "upload")) {
    states <- upload_states(); key <- upload_key(); state <- states[[key]]
    state$data <- d; state$table_id <- id; states[[key]] <- state; upload_states(states)
  }
  id
}

observeEvent(input$save_table, guard({ id <- save_current_table(); load_table(id); showNotification("Table saved in this project.") }))
observeEvent(input$copy_table, guard({ id <- save_current_table(TRUE, paste(input$table_name, "copy")); load_table(id) }))
observeEvent(input$source, { if (!input$source %in% c("workspace", "upload")) { active_table_id(NULL); working_table_id(NULL); selected_sheet(NULL) } }, ignoreInit = TRUE)
# Open project information without adding another workspace navigation tab.
show_project_details <- function() {
  p <- project()
  showModal(modalDialog(title = "Project information",
    textInput("project_name", "Project title", p$name),
    textAreaInput("project_info", "Notes / experiment information", p$info, rows = 8),
    footer = tagList(modalButton("Close"), actionButton("save_info", "Save information", class = "btn-primary")),
    easyClose = TRUE))
}
observeEvent(input$project_details, show_project_details())
observeEvent(input$save_info, {
  p <- project(); p$name <- input$project_name; p$info <- input$project_info; project(p)
  removeModal()
  showNotification("Project information saved.")
})

project_button <- function(kind, id, name) {
  active <- identical(selected_sheet(), list(kind = kind, id = id))
  tags$button(type = "button",
  class = paste("project-link", if (active) "selected" else ""), `data-kind` = kind, `data-key` = id,
  `aria-pressed` = if (active) "true" else "false",
  title = name, span(class = "sheet-icon", switch(kind, tables = "\u25a6", analyses = "\u25a4", graphs = "\u2301", layouts = "\u25a3", "\u24d8")), name)
}

output$project_tree <- renderUI({
  p <- project(); search <- if (is.null(input$project_search)) "" else input$project_search
  tagList(lapply(c("tables", "analyses", "graphs", "layouts"), function(section) {
    entries <- p[[section]]
    if (nzchar(search)) entries <- Filter(function(e) grepl(tolower(search), tolower(e$name), fixed = TRUE), entries)
    tags$details(open = "open", class = "project-section", tags$summary(paste(switch(section, tables = "Data Tables", analyses = "Results", graphs = "Graphs", layouts = "Layouts"), length(entries))),
      lapply(entries, function(e) project_button(section, e$id, e$name)),
      if (!length(entries)) div(class = "hint", if (nzchar(search)) "No matches" else "No saved sheets yet"))
  }))
})
output$family_tree <- renderUI({
  p <- project(); id <- active_table_id(); t <- if (is.null(id)) NULL else p$tables[[id]]
  if (is.null(t)) return(p(class = "hint", "Save or open a table to view related sheets."))
  analyses <- Filter(function(a) identical(a$table_id, id), p$analyses)
  graphs <- Filter(function(g) g$analysis_id %in% names(analyses), p$graphs)
  tagList(project_button("tables", id, t$name), div(class = "family-children",
    lapply(analyses, function(a) project_button("analyses", a$id, a$name)),
    lapply(graphs, function(g) project_button("graphs", g$id, g$name))))
})

open_analysis <- function(id) {
  # Archive view intentionally bypasses live-input freshness checks: display the
  # stored snapshot while retaining the user's separate working worksheet.
  a <- project()$analyses[[id]]; req(a)
  saved(a); archive_view(TRUE); active_table_id(a$table_id)
  select_sheet("analyses", id)
  updateTextInput(session, "analysis_name", value = a$name)
  updateTabsetPanel(session, "tabs", selected = "results")
}
observeEvent(input$project_open, guard({
  event <- input$project_open; p <- project()
  switch(event$kind,
    tables = load_table(event$key),
    analyses = open_analysis(event$key),
    graphs = {
      g <- p$graphs[[event$key]]; req(g); open_analysis(g$analysis_id)
      updateSelectInput(session, "graph_kind", selected = g$kind)
      updateTextInput(session, "graph_name", value = g$name)
      ids <- c(style = "plot_style", error = "error_style", title = "plot_title")
      options <- modifyList(lapply(formals(make_plot)[-1], eval), g$options)
      for (n in names(options)) {
        value <- options[[n]]; id <- if (n %in% names(ids)) ids[[n]] else n
        if (n %in% c("connect", "varied_symbols", "show_stars")) updateCheckboxInput(session, id, value = value)
        else if (n == "point_size") updateNumericInput(session, id, value = value)
        else if (n %in% c("title", "x_label", "y_label")) updateTextInput(session, id, value = value)
        else updateSelectInput(session, id, selected = value)
      }
      updateTabsetPanel(session, "tabs", selected = "figure")
      select_sheet("graphs", event$key)
    },
    layouts = {
      l <- p$layouts[[event$key]]; req(l)
      updateTextInput(session, "layout_name", value = l$name)
      updateSelectInput(session, "layout_graphs", selected = l$graph_ids)
      updateSelectInput(session, "layout_mode", selected = l$mode)
      updateNumericInput(session, "layout_columns", value = l$columns)
      updateTabsetPanel(session, "tabs", selected = "layouts")
      select_sheet("layouts", event$key)
    }, info = show_project_details())
}))
observeEvent(input$rename_analysis, guard({
  req(saved()); assert(nzchar(trimws(input$analysis_name)), "Enter an analysis name.")
  p <- project(); a <- saved(); a$name <- trimws(input$analysis_name)
  p$analyses[[a$id]]$name <- a$name; project(p); saved(a)
}))

output$baseline_ui <- renderUI({
  req(input$layout == "wide", input$wide_columns)
  selectInput("baseline", "Baseline column", input$wide_columns, selected = input$wide_columns[1])
})
observeEvent(input$transform_table, guard({
  assert(input$table_type == "column" && input$layout == "wide", "Baseline transformations use Column tables with one group per column.")
  source_id <- save_current_table()
  t <- normalize_wide(raw_data(), input$wide_columns, input$baseline, input$transform_operation)
  id <- save_current_table(TRUE, input$derived_name, t$data, list(source_table = source_id, operation = input$transform_operation, baseline = input$baseline, note = t$note))
  load_table(id); showNotification(paste(t$note, sum(t$audit$invalid_values_blanked), "undefined values blanked."), duration = 10)
}))
observeEvent(input$save_cleaned, guard({
  r <- result(); assert(r$method == "outliers", "Open an outlier analysis first.")
  p <- project(); a <- saved(); t <- a$worksheet
  if (is.null(t)) t <- p$tables[[a$table_id]]
  # Reconstruct from the immutable analysis input, even if its source worksheet was subsequently edited.
  if (t$layout == "wide") {
    d <- t$data
    assert(nrow(a$input_data) == nrow(d)*length(t$columns), "Source shape changed; export the cleaned-data result instead.")
    for (j in seq_along(t$columns)) d[[t$columns[j]]] <- a$result$tables$cleaned_data[[a$config$y]][((j-1)*nrow(d)+1):(j*nrow(d))]
  } else d <- r$tables$cleaned_data
  t$data <- d; t$name <- paste(t$name, "cleaned"); t$provenance <- list(analysis = a$id, Q = r$config$rout_q)
  added <- project_add(p, "tables", t); project(added$project); load_table(added$id)
}))

observeEvent(input$save_graph, guard({
  figure(); p <- project(); req(saved()$id)
  name <- if (is.null(input$graph_name) || !nzchar(trimws(input$graph_name))) paste(names(graph_kinds)[match(input$graph_kind, graph_kinds)], saved()$name) else trimws(input$graph_name)
  added <- project_add(p, "graphs", list(name = name, analysis_id = saved()$id, kind = input$graph_kind, options = plot_options()))
  project(added$project); select_sheet("graphs", added$id); showNotification("Graph saved.")
}))
observeEvent(project(), {
  graphs <- project()$graphs
  choices <- setNames(as.character(names(graphs)), vapply(graphs, function(g) g$name, character(1)))
  updateSelectInput(session, "layout_graphs", choices = choices, selected = intersect(isolate(input$layout_graphs), names(graphs)))
})
layout_entries <- reactive({
  # Saved graph settings point to the exact analysis snapshot used to create them.
  ids <- input$layout_graphs; req(length(ids))
  p <- project(); assert(all(ids %in% names(p$graphs)), "Select saved graphs.")
  lapply(ids, function(id) { g <- p$graphs[[id]]; g$result <- p$analyses[[g$analysis_id]]$result; g })
})
draw_layout <- function() {
  entries <- layout_entries()
  if (input$layout_mode == "overlay") print(overlay_graphs(entries, input$layout_name))
  else draw_graph_layout(lapply(entries, function(g) make_graph(g$result, g$kind, g$options)), input$layout_columns, input$layout_name)
}
output$layout_figure <- renderPlot(draw_layout(), res = 100)
observeEvent(input$save_layout, guard({
  entries <- layout_entries()
  if (input$layout_mode == "overlay") overlay_graphs(entries)
  p <- project(); added <- project_add(p, "layouts", list(name = input$layout_name, graph_ids = input$layout_graphs, mode = input$layout_mode, columns = as.integer(input$layout_columns)))
  project(added$project); select_sheet("layouts", added$id); showNotification("Layout saved.")
}))
output$download_layout_png <- downloadHandler(filename = "layout.png", content = function(file) {
  grDevices::png(file, width = 10, height = 7, units = "in", res = 300); on.exit(grDevices::dev.off()); draw_layout()
})
output$download_layout_pdf <- downloadHandler(filename = "layout.pdf", content = function(file) {
  grDevices::cairo_pdf(file, width = 10, height = 7); on.exit(grDevices::dev.off()); draw_layout()
})
output$download_layout_svg <- downloadHandler(filename = "layout.svg", content = function(file) {
  svglite::svglite(file, width = 10, height = 7); on.exit(grDevices::dev.off()); draw_layout()
})
output$download_project <- downloadHandler(filename = "labstats_project.json", content = function(file) writeLines(project_json(project()), file, useBytes = TRUE))
observeEvent(input$project_file, guard({
  req(input$project_file); restored <- read_project(paste(readLines(input$project_file$datapath, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  project(restored); saved(NULL); edited(NULL); upload_states(list()); opened_table(NULL); active_table_id(NULL); working_table_id(NULL); archive_view(FALSE)
  selected_sheet(NULL)
  updateTextInput(session, "project_name", value = restored$name); updateTextAreaInput(session, "project_info", value = restored$info)
  if (length(restored$tables)) load_table(names(restored$tables)[1])
  else updateSelectInput(session, "source", selected = "demo")
  showNotification("Project opened; saved analyses recalculated.")
}))
