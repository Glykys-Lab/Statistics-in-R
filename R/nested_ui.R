# UI constructors and configuration parsing for Mouse/Slice/Cell nested models.
# Settings mirror nested_defaults(); the server handles column mappings, selected
# comparisons and restoration. Label text may use literal backslash-n line breaks.
nested_control <- function(name, label, choices = NULL, min = NA, max = NA, step = NA) {
  value <- nested_defaults()[[name]]; id <- paste0("nested_", name)
  if (!is.null(choices)) return(selectInput(id, label, choices, selected = value, multiple = length(value) > 1L || name == "pub_formats"))
  if (is.logical(value)) return(checkboxInput(id, label, value))
  if (is.numeric(value)) return(numericInput(id, label, value, min = min, max = max, step = step))
  textInput(id, label, if (is.null(value)) "" else value)
}
nested_options_ui <- function() tagList(
  h3("Linear mixed analysis: mouse / slice / cell"),
  p("Use these models when cells share slices and mice. Repeated data have one row per cell and a column per condition; between-group data have one measurement and condition per row. Map the ID columns in Analysis."),
  tabsetPanel(id = "nested_tabs",
    tabPanel("Design & fitting", value = "fit",
      radioButtons("nested_design", "Experimental design", c("Same cells followed across conditions (wide table)" = "repeated", "Different groups, one measurement per cell (long table)" = "between")),
      nested_control("log_transform", "Analyze log10 responses and report ratios / geometric means"),
      p(class = "hint", "Log10 requires strictly positive measurements. This path uses Satterthwaite inference rather than the single-subject path's Geisser-Greenhouse options."),
      conditionalPanel("input.nested_design == 'repeated'", nested_control("slice_condition", "Slice x condition random effect", c("Automatic: shared or separate SD" = "auto", "Shared SD" = "shared", "Separate SD per condition" = "per_condition", "Omit slice x condition" = "none"))),
      conditionalPanel("input.nested_design == 'between'", nested_control("slice_variance", "Slice variability within mice", c("Automatic" = "auto", "Shared SD" = "shared", "Separate SD per condition" = "per_condition"))),
      nested_control("model_engine", "Residual / cell noise", c("Automatic: compare shared and unequal noise" = "auto", "lmer: shared residual SD" = "lmer", "nlme: separate residual SD per condition" = "nlme")),
      nested_control("variance_alpha", "Threshold for automatic variance tests", min = .001, max = .5, step = .01),
      nested_control("df_method", "Satterthwaite degrees of freedom", c("Each mean / comparison uses its own df" = "per_comparison", "Smallest df within each set of means / contrasts" = "smallest")),
      nested_control("check_variance_homogeneity", "Include candidate variance tables and likelihood-ratio tests"),
      p(class = "hint", "Automatic variance choices follow the supplied scripts' conventional likelihood-ratio approximations. Fitting checks identify candidate-model and selected-model warnings.")),
    tabPanel("Comparisons", value = "comparisons",
      nested_control("comparison_type", "Which condition comparisons?", c("Every pair (Tukey by default)" = "all_pairs", "Each versus a control (Dunnett approximation)" = "vs_first", "Consecutive conditions (Holm by default)" = "consecutive", "Selected pairs (Holm by default)" = "custom")),
      conditionalPanel("input.nested_comparison_type == 'vs_first'", uiOutput("nested_control_ui")),
      conditionalPanel("input.nested_comparison_type == 'custom'", uiOutput("nested_pairs_ui")),
      nested_control("p_adjust", "Comparison correction", c("Use the comparison family's default" = "auto", "Holm p-values / Bonferroni intervals" = "holm", "Bonferroni" = "bonferroni", "No correction" = "none")),
      nested_control("compare_to_zero", "Also test each condition mean against zero"),
      p(class = "hint", "With log10 enabled, zero on the model scale corresponds to a geometric mean of 1. These mean tests follow the scripts and are unadjusted; condition comparisons use the selected family correction. Intervals are 95%."),
      textAreaInput("nested_condition_order", "Optional condition order (one name per line; blank preserves table order)", "", rows = 3)),
    tabPanel("Outliers", value = "outliers",
      nested_control("outlier_action", "ROUT cell outliers", c("Flag and keep measurements" = "flag", "Remove flagged cell measurements from this analysis" = "remove", "Do not test outliers" = "none")),
      nested_control("outlier_unit", "Units to check", c("Cells" = "cell", "Slices (flag only)" = "slice", "Cells and slices" = "both")),
      nested_control("outlier_level", "Values to test", c("Linear mixed model residuals: remove mouse / slice offsets" = "residuals", "Raw values within each condition" = "raw")),
      nested_control("outlier_Q", "ROUT false discovery rate Q (0.01 = 1%)", min = .001, max = .1, step = .001),
      p(class = "hint", "Removal affects this analysis snapshot; the source table is preserved. Slice outliers are always flagged and kept. Tables and reports record the action and each flagged measurement.")),
    tabPanel("Figures", value = "figures",
      nested_control("summary_display", "Summary over the cells", c("Model mean and 95% CI" = "model", "Raw median and interquartile range" = "median_iqr")),
      nested_control("show_slice_means", "Show slice means on the between-group notebook plot"),
      nested_control("show_p_on_plot", "Draw adjusted p-values on the main plot"),
      nested_control("p_threshold", "Draw comparisons below this p-value", min = .001, max = 1, step = .01),
      nested_control("comparison_plot", "Comparison interval plot", c("Create a separate plot" = "separate", "Do not create" = "none")),
      nested_control("show_diagnostics", "Create residual and condition diagnostics"),
      nested_control("y_label", "Notebook y-axis title"),
      nested_control("change_measure", "Dashed reference line at zero", c("Automatic: change values or tests against zero" = "auto", "Always" = "TRUE", "Never" = "FALSE")),
      nested_control("y_include_zero", "Include zero on the raw-scale y-axis"),
      nested_control("make_pub_figures", "Create publication SuperPlot with key and comparison strip"),
      tags$details(open = "open", tags$summary("Publication dimensions and summaries"),
        nested_control("pub_formats", "Formats included in the download bundle", c("SVG" = "svg", "PDF" = "pdf", "EPS" = "eps")),
        fluidRow(column(4, nested_control("pub_width_mm", "Width (mm)", min = 10, step = 1)), column(4, nested_control("pub_main_height_mm", "Main panel height (mm)", min = 10, step = 1)), column(4, nested_control("pub_comp_height_mm", "Comparison strip height (mm; 0 hides it)", min = 0, step = 1))),
        nested_control("pub_summary", "Model summary marker", c("Mean and bars" = "mean_ci", "Mean only" = "mean", "Hide summary" = "none")),
        nested_control("pub_error_bar", "Publication error bars", c("95% model CI" = "ci", "Model SEM" = "sem")),
        nested_control("pub_cell_layout", "Cell layout", c("Beeswarm" = "beeswarm", "Jitter" = "jitter")),
        nested_control("pub_y_title", "Publication y-axis title (blank uses notebook title)"), nested_control("pub_x_title", "Publication x-axis title"),
        textAreaInput("nested_pub_condition_labels", "Condition display labels: original = display (one per line; use \\n for a line break)", "", rows = 3),
        textAreaInput("nested_pub_mouse_colours", "Mouse colours: mouse ID = #hex (blank uses Okabe-Ito)", "", rows = 3),
        nested_control("pub_micro", "Convert u-prefixed units to the micro sign"), nested_control("pub_ion_superscript", "Raise ion charges as superscripts")),
      tags$details(tags$summary("Publication key, fonts, strokes and symbols"),
        nested_control("pub_key_width_mm", "Key width in mm (blank fits labels; 0 hides key)"),
        textAreaInput("nested_pub_key_labels", "Key labels: cell / slice / mouse = text", "cell = Cell\nslice = Slice\nmouse = Colored, Mouse", rows = 3),
        nested_control("pub_font", "Font family"), nested_control("pub_dash", "Dash / gap pattern"),
        fluidRow(lapply(c("pub_key_nudge_mm", "pub_key_pt", "pub_title_pt", "pub_tick_pt", "pub_p_pt", "pub_small_pt", "pub_axis_pt", "pub_dot_pt", "pub_ci_pt", "pub_mean_pt", "pub_slice_line_pt", "pub_cell_size", "pub_slice_size", "pub_cell_alpha", "pub_summary_offset", "pub_mean_width", "pub_ci_cap_width"), function(n) column(4, nested_control(n, gsub("_", " ", sub("^pub_", "", n)), step = .05)))))),
    tabPanel("Reports & export", value = "reports",
      nested_control("num_decimals", "Numeric decimal places in script reports", min = 0, max = 10, step = 1),
      nested_control("p_decimals", "P-value decimal places in script reports", min = 1, max = 10, step = 1),
      p("Results include decisions, candidate models, model means, comparisons, degrees-of-freedom sources, counts of mice / slices / cells, variance tests, normality, fitting checks, outliers and a writing template."),
      p("The Export & R code tab provides a self-contained HTML report with section navigation, an Obsidian Markdown note and a ZIP containing figures and result tables. The download uses your browser's chosen folder.")))
)

nested_named_lines <- function(text) {
  # Convert human-readable name/value lines into named vectors while rejecting
  # duplicate keys that would make figure labels or mouse colours ambiguous.
  if (is.null(text) || !nzchar(trimws(text))) return(NULL)
  lines <- trimws(strsplit(text, "\n", fixed = TRUE)[[1]]); lines <- lines[nzchar(lines)]
  pairs <- strsplit(lines, "=", fixed = TRUE)
  assert(all(lengths(pairs) == 2L), "Use one 'name = value' entry per line for figure labels and colours.")
  keys <- vapply(pairs, function(p) trimws(p[1]), character(1))
  assert(all(nzchar(keys)) && !anyDuplicated(keys), "Figure labels and colours need unique nonempty names.")
  setNames(vapply(pairs, function(p) gsub("\\n", "\n", trimws(p[2]), fixed = TRUE), character(1)), keys)
}

nested_config <- function(input, data) {
  # Normalize optional fields and pair encodings without moving the control to
  # the front of the worksheet; experimental sequence and reference are separate.
  settings <- nested_defaults()
  for (n in setdiff(names(settings), "custom_pairs")) {
    v <- input[[paste0("nested_", n)]]
    if (!is.null(v)) settings[n] <- list(v)
  }
  for (n in c("pub_y_title", "pub_x_title")) if (!nzchar(trimws(settings[[n]] %||% ""))) settings[n] <- list(NULL)
  for (n in c("pub_key_labels", "pub_condition_labels", "pub_mouse_colours")) settings[n] <- list(nested_named_lines(input[[paste0("nested_", n)]]))
  if (is.null(settings$pub_key_labels)) settings$pub_key_labels <- nested_defaults()$pub_key_labels
  if (nzchar(settings$pub_key_width_mm %||% "")) {
    settings$pub_key_width_mm <- suppressWarnings(as.numeric(settings$pub_key_width_mm))
    assert(is.finite(settings$pub_key_width_mm) && settings$pub_key_width_mm >= 0, "Key width must be blank or a nonnegative number.")
  } else settings["pub_key_width_mm"] <- list(NULL)
  if (identical(settings$change_measure, "TRUE")) settings$change_measure <- TRUE
  if (identical(settings$change_measure, "FALSE")) settings$change_measure <- FALSE
  settings$custom_pairs <- lapply(input$nested_pairs, function(k) vapply(strsplit(k, "::", fixed = TRUE)[[1]], utils::URLdecode, character(1)))
  order <- trimws(strsplit(input$nested_condition_order %||% "", "\n", fixed = TRUE)[[1]]); order <- order[nzchar(order)]
  method <- if (identical(input$nested_design, "between")) "nested_between" else "nested_repeated"
  list(method = method, mouse = input$nested_mouse, slice = input$nested_slice, cell = input$nested_cell,
    columns = input$nested_columns, condition = input$nested_condition, value = input$nested_value,
    condition_order = order, control = if (settings$comparison_type == "vs_first") input$nested_control else NULL,
    source_name = input$table_name, nested = settings)
}

`%||%` <- function(x, fallback) if (is.null(x)) fallback else x
