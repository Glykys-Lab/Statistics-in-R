# Sourced within the server; reactive mappings remain separate from Column/XY.
# Build nested ID/condition selectors, restore saved controls, and expose report
# and figure downloads. Pending selections survive hidden Shiny output suspension.
pending_nested_control <- reactiveVal(NULL)
pending_nested_pairs <- reactiveVal(NULL)
nested_pair_keys <- function(pairs) vapply(lapply(pairs, unlist), function(p) paste(vapply(p, utils::URLencode, character(1), reserved = TRUE), collapse = "::"), character(1))
nested_mapping_ui <- function(d, input) {
  n <- names(d); nums <- n[vapply(d, is.numeric, logical(1))]
  pick <- function(field, preferred, choices = n) {
    previous <- isolate(input[[field]])
    if (length(previous) == 1 && previous %in% choices) previous else {
      hit <- intersect(preferred, choices); if (length(hit)) hit[1] else if (length(choices)) choices[1] else character()
    }
  }
  mouse <- pick("nested_mouse", c("Mouse", "mouse"))
  slice <- pick("nested_slice", c("Slice", "slice"), setdiff(n, mouse))
  cell <- pick("nested_cell", c("Cell", "cell"), setdiff(n, c(mouse, slice)))
  ids <- c(mouse, slice, cell)
  choices <- setdiff(nums, ids)
  previous <- isolate(input$nested_columns)
  tagList(selectInput("nested_mouse", "Mouse ID column", n, selected = ids[1]),
    selectInput("nested_slice", "Slice ID column", n, selected = ids[2]),
    selectInput("nested_cell", "Cell ID column", n, selected = ids[3]),
    if (!identical(input$nested_design, "between")) selectInput("nested_columns", "Condition measurement columns", setdiff(n, ids), selected = if (length(previous) && all(previous %in% choices)) previous else choices, multiple = TRUE)
    else tagList(selectInput("nested_condition", "Condition column", n, selected = pick("nested_condition", c("Condition", "condition", "group"))),
      selectInput("nested_value", "Measurement column", choices, selected = pick("nested_value", c("Value", "value", "response"), choices))))
}
nested_conditions <- reactive({
  d <- raw_data()
  if (identical(input$nested_design, "between")) {
    req(input$nested_condition %in% names(d)); unique(as.character(d[[input$nested_condition]]))
  } else names(d)[names(d) %in% input$nested_columns]
})
output$nested_control_ui <- renderUI({
  groups <- nested_conditions(); req(length(groups))
  selected <- pending_nested_control()
  if (is.null(selected)) selected <- isolate(input$nested_control)
  selectInput("nested_control", "Control condition", groups, selected = if (length(selected) == 1L && selected %in% groups) selected else groups[1])
})
output$nested_pairs_ui <- renderUI({
  groups <- nested_conditions(); req(length(groups) >= 2L)
  pairs <- unlist(lapply(combn(groups, 2, simplify = FALSE), function(p) list(p, rev(p))), recursive = FALSE)
  values <- vapply(pairs, function(p) paste(vapply(p, utils::URLencode, character(1), reserved = TRUE), collapse = "::"), character(1))
  selected <- pending_nested_pairs()
  if (is.null(selected)) selected <- isolate(input$nested_pairs)
  selectInput("nested_pairs", "Pairs to compare (first minus second)", setNames(values, vapply(pairs, paste, collapse = " minus ", character(1))), selected = intersect(selected, values), multiple = TRUE)
})
observeEvent(input$nested_control, {
  if (identical(input$nested_control, pending_nested_control())) pending_nested_control(NULL)
}, ignoreNULL = FALSE)
observeEvent(input$nested_pairs, {
  if (identical(sort(as.character(input$nested_pairs)), sort(pending_nested_pairs()))) pending_nested_pairs(NULL)
}, ignoreNULL = FALSE)
output$nested_report_view <- renderUI({
  r <- result(); if (!isTRUE(r$nested)) return(NULL)
  tags$details(tags$summary("Full linear mixed analysis report: figures, decisions and writing template"),
    tags$iframe(title = "Linear mixed analysis report", srcdoc = nested_html_report(r), sandbox = "allow-scripts", style = "width:100%;height:780px;border:1px solid #cddff1;margin-top:12px"))
})
output$publication_preview <- renderImage({
  # Reuse the generated PNG at its publication aspect ratio instead of stretching
  # a millimetre-sized grid panel to the ordinary plot viewport.
  r <- result(); req(isTRUE(r$nested), identical(input$graph_kind, "publication"))
  bytes <- r$artifacts$nested_analysis_figure.png
  validate(need(!is.null(bytes), "Enable publication figures and run the linear mixed analysis first."))
  file <- tempfile(fileext = ".png"); writeBin(bytes, file)
  height <- r$publication_dimensions$height
  list(src = file, contentType = "image/png", width = 480, height = round(480*height/r$config$nested$pub_width_mm), alt = "Publication SuperPlot preview")
}, deleteFile = TRUE)
output$nested_exports <- renderUI({
  r <- result(); if (!isTRUE(r$nested)) return(NULL)
  tagList(h4("Linear mixed analysis report and publication figures"), div(class = "download-row", downloadButton("download_nested_html", "HTML report"), downloadButton("download_nested_note", "Obsidian Markdown note"), downloadButton("download_nested_bundle", "All tables and figures (ZIP)")))
})
observe({
  r <- saved()$result
  kinds <- if (isTRUE(r$nested)) unique(c(names(r$graphs), "normal_qq", "lognormal_qq")) else setdiff(unname(graph_kinds), c("publication", "condition_residual_qq"))
  choices <- graph_kinds[graph_kinds %in% kinds]
  selected <- isolate(input$graph_kind)
  if (is.null(selected) || !selected %in% kinds) selected <- "raw"
  updateSelectInput(session, "graph_kind", choices = choices, selected = selected)
})
output$download_nested_html <- downloadHandler(filename = "linear_mixed_analysis_report.html", content = function(file) {
  r <- result(); req(isTRUE(r$nested)); writeLines(nested_html_report(r), file, useBytes = TRUE)
})
output$download_nested_note <- downloadHandler(filename = "linear_mixed_analysis.md", content = function(file) {
  r <- result(); req(isTRUE(r$nested)); writeLines(r$report_markdown, file, useBytes = TRUE)
})
output$download_nested_bundle <- downloadHandler(filename = "linear_mixed_analysis.zip", content = function(file) {
  r <- result(); req(isTRUE(r$nested)); write_nested_bundle(r, file)
})
restore_nested_config <- function(cfg) {
  # Retain choices until their dynamic widgets send the restored values back;
  # updates sent while a hidden tab is suspended can otherwise be lost.
  if (is.null(cfg) || !cfg$method %in% nested_methods) return()
  pending_nested_control(if (length(cfg$control)) as.character(unlist(cfg$control)) else NULL)
  pending_nested_pairs(nested_pair_keys(cfg$nested$custom_pairs))
  updateRadioButtons(session, "nested_design", selected = if (cfg$method == "nested_between") "between" else "repeated")
  for (n in c("mouse", "slice", "cell", "columns", "condition", "value")) updateSelectInput(session, paste0("nested_", n), selected = unlist(cfg[[n]]))
  updateTextAreaInput(session, "nested_condition_order", value = paste(unlist(cfg$condition_order), collapse = "\n"))
  defaults <- nested_defaults()
  for (n in intersect(names(cfg$nested), names(defaults))) {
    v <- cfg$nested[[n]]
    if (n == "custom_pairs") next
    if (n %in% c("pub_key_labels", "pub_condition_labels", "pub_mouse_colours")) {
      v <- unlist(v); v <- gsub("\n", "\\n", v, fixed = TRUE)
      updateTextAreaInput(session, paste0("nested_", n), value = paste(paste(names(v), v, sep = " = "), collapse = "\n"))
    } else if (is.logical(defaults[[n]])) updateCheckboxInput(session, paste0("nested_", n), value = isTRUE(v))
    else if (is.numeric(defaults[[n]])) updateNumericInput(session, paste0("nested_", n), value = v)
    else if (n %in% c("slice_condition", "slice_variance", "model_engine", "df_method", "comparison_type", "p_adjust", "outlier_action", "outlier_unit", "outlier_level", "comparison_plot", "summary_display", "pub_formats", "pub_summary", "pub_error_bar", "pub_cell_layout", "change_measure")) updateSelectInput(session, paste0("nested_", n), selected = as.character(unlist(v)))
    else updateTextInput(session, paste0("nested_", n), value = if (is.null(v)) "" else v)
  }
  session$onFlushed(function() {
    keys <- nested_pair_keys(cfg$nested$custom_pairs)
    updateSelectInput(session, "nested_pairs", selected = keys)
    if (length(cfg$control)) updateSelectInput(session, "nested_control", selected = unlist(cfg$control))
  }, once = TRUE)
}
