# Adapter for the user's nested-model reference scripts. Their model calculations
# stay in separate files; no personal input/output paths are used by LabStats.
# Validate hierarchical mappings/settings, execute a reference in an isolated
# environment with temporary paths, and adapt its outputs to LabStats result records.
# Also handle named-setting JSON round trips, report image embedding and bundles.
nested_methods <- c("nested_repeated", "nested_between")
nested_packages <- c("lme4", "lmerTest", "emmeans", "gridExtra", "tidyr", "dplyr", "nlme", "ggplot2", "commonmark", "svglite", "ggbeeswarm", "zip")
nested_defaults <- function() list(log_transform = FALSE, slice_condition = "auto", slice_variance = "auto",
  model_engine = "auto", variance_alpha = .05, check_variance_homogeneity = TRUE,
  compare_to_zero = FALSE, df_method = "per_comparison", comparison_type = "all_pairs", p_adjust = "auto",
  custom_pairs = list(), outlier_action = "flag", outlier_unit = "cell", outlier_Q = .01, outlier_level = "residuals",
  show_p_on_plot = TRUE, p_threshold = .05, comparison_plot = "separate", show_slice_means = TRUE,
  num_decimals = 2L, p_decimals = 3L, summary_display = "model", y_label = "Response",
  show_diagnostics = TRUE, change_measure = "auto", y_include_zero = TRUE, make_pub_figures = TRUE,
  pub_formats = c("svg", "pdf"), pub_width_mm = 80, pub_key_width_mm = NULL, pub_key_nudge_mm = 1,
  pub_key_labels = c(cell = "Cell", slice = "Slice", mouse = "Colored, Mouse"), pub_key_pt = 4,
  pub_main_height_mm = 50, pub_comp_height_mm = 20, pub_font = "Arial", pub_title_pt = 8, pub_tick_pt = 7,
  pub_p_pt = 5, pub_small_pt = 5, pub_axis_pt = .25, pub_dot_pt = .25, pub_ci_pt = .5, pub_mean_pt = 1,
  pub_slice_line_pt = .25, pub_cell_size = .8, pub_slice_size = 1.8, pub_cell_alpha = .4,
  pub_cell_layout = "beeswarm", pub_summary = "mean_ci", pub_error_bar = "ci", pub_summary_offset = 0,
  pub_mean_width = .25, pub_ci_cap_width = .1, pub_micro = TRUE, pub_ion_superscript = TRUE, pub_dash = "22",
  pub_y_title = NULL, pub_x_title = NULL, pub_condition_labels = NULL, pub_mouse_colours = NULL)

nested_reference_text <- function(method) {
  key <- if (method == "nested_between") "between-groups" else "repeated"
  if (exists("embedded_nested_references", inherits = TRUE)) return(get("embedded_nested_references", inherits = TRUE)[[key]])
  readLines(file.path("R", paste0("nested-reference-", key, ".inc")), warn = FALSE, encoding = "UTF-8")
}

nested_input <- function(data, cfg) {
  # Canonical IDs scope reused slice/cell labels by their parents. Repeated input
  # is wide; between-group input is long and a slice cannot span conditions.
  ids <- unlist(cfg[c("mouse", "slice", "cell")], use.names = FALSE)
  assert(length(ids) == 3L && all(ids %in% names(data)) && !anyDuplicated(ids), "Select three different Mouse, Slice and Cell ID columns.")
  out <- data[ids]; names(out) <- c("Mouse", "Slice", "Cell")
  for (id in names(out)) {
    out[[id]] <- trimws(as.character(out[[id]]))
    assert(all(!is.na(out[[id]]) & nzchar(out[[id]])), paste(id, "IDs cannot be blank."))
  }
  if (cfg$method == "nested_repeated") {
    columns <- as.character(unlist(cfg$columns))
    columns <- names(data)[names(data) %in% columns]
    assert(length(columns) > 0 && !any(columns %in% ids), "Choose the condition measurement columns, excluding IDs.")
    for (n in columns) numeric_column(data, n)
    assert(!anyDuplicated(out), "Each Mouse / Slice / Cell combination must occupy one row in a repeated table.")
    out <- cbind(out, data[columns])
    all_values <- unlist(data[columns], use.names = FALSE)
    conditions <- columns
  } else {
    assert(length(cfg$condition) == 1 && cfg$condition %in% names(data) && !cfg$condition %in% ids, "Select the Condition column.")
    assert(length(cfg$value) == 1 && cfg$value %in% names(data) && !cfg$value %in% c(ids, cfg$condition), "Select the measurement column.")
    out$Condition <- trimws(as.character(data[[cfg$condition]]))
    assert(all(!is.na(out$Condition) & nzchar(out$Condition)), "Condition names cannot be blank.")
    out$Value <- numeric_column(data, cfg$value)
    assert(!anyDuplicated(out[c("Mouse", "Slice", "Cell")]), "Each cell must be measured once in a between-group table.")
    slice_conditions <- split(out$Condition, interaction(out$Mouse, out$Slice, drop = TRUE))
    assert(all(vapply(slice_conditions, function(x) length(unique(x)) == 1L, logical(1))), "A slice appears under several conditions. Choose the repeated-measures design.")
    all_values <- out$Value; conditions <- unique(out$Condition)
  }
  order <- as.character(unlist(cfg$condition_order))
  if (length(order)) {
    assert(setequal(order, conditions) && !anyDuplicated(order), "Condition order must list every condition exactly once.")
    conditions <- order
    if (cfg$method == "nested_repeated") out <- out[c("Mouse", "Slice", "Cell", order)]
  }
  assert(any(!is.na(all_values)), "Enter at least one measurement before running a linear mixed analysis.")
  out
}

run_nested_analysis <- function(data, config, progress = function(value = NULL, detail = NULL) NULL) {
  # Reference scripts expect named SETTINGS variables and CSV files. Override
  # both before their execution; capture warnings and outputs for this snapshot.
  progress(.05, "Checking data and model settings")
  missing <- nested_packages[!vapply(nested_packages, requireNamespace, logical(1), quietly = TRUE)]
  assert(!length(missing), paste("Run source('setup.R') to install linear mixed analysis packages:", paste(missing, collapse = ", ")))
  assert(is.null(config$nested) || (is.list(config$nested) && all(names(config$nested) %in% names(nested_defaults()))), "Unknown linear mixed analysis settings.")
  settings <- modifyList(nested_defaults(), if (is.null(config$nested)) list() else config$nested, keep.null = TRUE)
  for (n in c("custom_pairs", "pub_key_labels", "pub_condition_labels", "pub_mouse_colours")) if (n %in% names(config$nested)) settings[n] <- list(config$nested[[n]])
  # JSON round-trips may represent vectors as lists.
  for (n in c("pub_formats", "pub_key_labels", "pub_condition_labels", "pub_mouse_colours")) if (!is.null(settings[[n]])) settings[[n]] <- unlist(settings[[n]])
  settings$custom_pairs <- lapply(settings$custom_pairs, function(p) as.character(unlist(p)))
  d <- nested_input(data, config)
  if (!is.null(settings$pub_mouse_colours)) {
    colours <- settings$pub_mouse_colours
    assert(!is.null(names(colours)) && all(names(colours) %in% unique(d$Mouse)), "Custom colours must name Mouse IDs present in the table.")
    valid <- tryCatch({ grDevices::col2rgb(colours); TRUE }, error = function(e) FALSE)
    assert(valid, "Use valid colour names or hexadecimal colours for mice.")
    mice <- sort(unique(d$Mouse))
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    filled <- setNames(rep(palette, length.out = length(mice)), mice)
    filled[names(colours)] <- colours
    settings$pub_mouse_colours <- filled
  }
  assert(!is.null(names(settings$pub_key_labels)) && all(names(settings$pub_key_labels) %in% c("cell", "slice", "mouse")), "Key labels must use cell, slice or mouse as their names.")
  labels <- nested_defaults()$pub_key_labels
  labels[names(settings$pub_key_labels)] <- settings$pub_key_labels
  settings$pub_key_labels <- labels
  values <- if (config$method == "nested_repeated") unlist(d[setdiff(names(d), c("Mouse", "Slice", "Cell"))], use.names = FALSE) else d$Value
  observed <- if (config$method == "nested_repeated") vapply(d[setdiff(names(d), c("Mouse", "Slice", "Cell"))], function(x) sum(!is.na(x)), integer(1)) else tapply(!is.na(d$Value), d$Condition, sum)
  assert(all(observed > 0), "Every selected condition must contain at least one observed measurement.")
  conditions <- if (config$method == "nested_repeated") names(observed) else unique(d$Condition)
  if (identical(settings$comparison_type, "custom")) {
    valid <- vapply(settings$custom_pairs, function(p) length(p) == 2L && all(p %in% conditions) && p[1] != p[2], logical(1))
    assert(length(valid) > 0 && all(valid), "Select pairs of two different observed conditions.")
    keys <- vapply(settings$custom_pairs, function(p) paste(sort(p), collapse = "\r"), character(1))
    assert(!anyDuplicated(keys), "Choose only one direction for each condition pair.")
  }
  control <- if (length(config$control)) as.character(config$control) else conditions[1]
  assert(length(control) == 1L && control %in% conditions, "Choose a control condition present in the selected measurements.")
  if (isTRUE(settings$log_transform)) assert(all(values[!is.na(values)] > 0), "Log10 linear mixed analysis requires every observed value to be positive. Zeros and negatives cannot be logged.")
  assert(settings$variance_alpha > 0 && settings$variance_alpha < 1 && settings$outlier_Q > 0 && settings$outlier_Q < 1, "Variance alpha and ROUT Q must lie between 0 and 1.")
  assert(settings$p_threshold > 0 && settings$p_threshold < 1 && settings$pub_width_mm > 0 && settings$pub_main_height_mm > 0 && settings$pub_comp_height_mm >= 0, "Choose valid significance and positive figure dimensions.")
  assert(settings$num_decimals %in% 0:10 && settings$p_decimals %in% 1:10, "Choose 0-10 numeric and 1-10 p-value decimal places.")
  progress(.1, "Preparing measurements")
  work <- tempfile("labstats-nested-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  # Font/layout measurement must not create Rplots.pdf in the project folder.
  grDevices::cairo_pdf(file.path(work, "layout-measurement.pdf"))
  device <- grDevices::dev.cur()
  on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE, after = FALSE)
  csv <- "nested_analysis.csv"
  write.csv(d, file.path(work, csv), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  overrides <- c(settings, list(r_csv_dir = work, csv_name = csv, results_dir = work, plot_dir = work, pub_dir = work,
    experiments_dir = work, interactive_prompts = FALSE, pick_csv_each_run = FALSE, pub_ask_titles = FALSE,
    show_plots = FALSE, save_outputs = TRUE, html_report = FALSE, value_col = "Value",
    condition_order = if (config$method == "nested_between") if (length(config$condition_order)) as.character(unlist(config$condition_order)) else unique(d$Condition) else NULL))
  env <- new.env(parent = globalenv()); env$labstats_overrides <- overrides
  # Keep reference variables out of the app environment and restore emmeans
  # options after the run so subsequent analyses retain their own inference settings.
  env$labstats_normality <- nested_normality
  env$labstats_progress <- progress
  env$labstats_control <- control
  if (identical(settings$change_measure, "auto") && length(config$source_name) && grepl("change", config$source_name, ignore.case = TRUE)) env$labstats_overrides$change_measure <- TRUE
  emm_options <- emmeans::emm_options()
  on.exit(do.call(emmeans::emm_options, emm_options), add = TRUE)
  warnings <- character()
  console <- tryCatch(capture.output(withCallingHandlers(eval(parse(text = nested_reference_text(config$method)), env),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))), error = function(e) {
      checks <- if (exists("fit_warnings", env, inherits = FALSE)) env$fit_warnings else character()
      stop(paste(c(conditionMessage(e), tail(unique(c(checks, warnings)), 4L)), collapse = "\n"), call. = FALSE)
    })
  progress(.9, "Collecting result tables and figures")
  get_optional <- function(n) if (exists(n, env, inherits = FALSE)) get(n, env, inherits = FALSE) else NULL
  # Preserve numeric tables alongside display-formatted copies and retain artifact
  # bytes before deleting the temporary directory used by the reference script.
  objects <- c(analysis_decisions = "decisions_tbl", candidate_models = "models_tbl", model_means = "estimates_tbl",
    multiple_comparisons = "comparisons_tbl", degrees_of_freedom = "df_breakdown_tbl", descriptives = "descriptives_tbl",
    slice_variances = "slice_sd_tbl", cell_variances = "cell_sd_tbl", slice_model_tests = "slice_test", residual_model_tests = "variance_test",
    random_effects = "random_tbl", fixed_effects = "fixed_tbl", normality = "normality_tbl", residual_normality = "resid_normality_tbl",
    condition_normality = "condition_normality_tbl", cell_outliers = "outlier_tbl", slice_outliers = "slice_outlier_tbl")
  tables <- lapply(objects, function(n) { x <- get_optional(n); if (is.null(x)) NULL else as.data.frame(x) })
  tables <- tables[!vapply(tables, is.null, logical(1))]
  if (!isTRUE(settings$check_variance_homogeneity)) tables <- tables[!names(tables) %in% c("slice_variances", "cell_variances", "slice_model_tests", "residual_model_tests")]
  tables$fitting_checks <- data.frame(Check = get_optional("fit_check_lines"))
  tables$writing_template <- data.frame(Section = c("Methods", "Results", "Figure legend", "Checklist"),
    Text = vapply(list(get_optional("methods_text"), c(get_optional("means_text"), get_optional("results_text")), get_optional("legend_text"), get_optional("caution_text")), paste, collapse = "\n", character(1)))
  long <- env$df_long
  tables$data_audit <- data.frame(input_rows = length(values), used_rows = nrow(long), excluded_rows = length(values) - nrow(long), missing_values = sum(is.na(values)), removed_outliers = if (is.null(get_optional("outliers_removed"))) 0L else env$outliers_removed)
  tables$analyzed_measurements <- long[c("Mouse", "Slice", "Cell", "Condition", "Value")]
  graph_objects <- c(raw = "p", comparison = "p_comp", publication = "pub_figure", residual = "residual_diagnostics", residual_qq = "p_qq_resid", condition_residual_qq = "resid_by_condition", assumption_qq = "condition_qq")
  graphs <- lapply(graph_objects, get_optional); graphs <- graphs[!vapply(graphs, is.null, logical(1))]
  if (!isTRUE(settings$show_diagnostics)) graphs <- graphs[!names(graphs) %in% c("residual", "residual_qq", "condition_residual_qq", "assumption_qq")]
  if (!is.null(graphs$publication)) graphs$publication <- grid::grobTree(graphs$publication,
    vp = grid::viewport(width = grid::unit(settings$pub_width_mm, "mm"), height = grid::unit(env$pub_total_height_mm, "mm")))
  paths <- list.files(work, full.names = TRUE)
  paths <- paths[!grepl("\\.csv$", paths) & basename(paths) != "layout-measurement.pdf"]
  artifacts <- setNames(lapply(paths, function(p) readBin(p, "raw", n = file.info(p)$size)), basename(paths))
  config$nested <- settings
  if (settings$log_transform && settings$compare_to_zero) {
    warnings <- c(warnings, "On the log10 scale, comparison to zero tests geometric means against 1, not original-scale means against 0.")
    tables$mean_tests <- data.frame(Condition = as.character(env$emm_df$Condition), null_geometric_mean = 1, p_value = env$emm_df$p.value)
  }
  env$count_cols_all <- c(env$count_cols_all, "input_rows", "used_rows", "excluded_rows", "missing_values", "removed_outliers")
  formatted <- lapply(tables, function(t) if (is.data.frame(t) && ncol(t)) env$fmt_tbl(t) else t)
  notes <- c(env$model_desc, get_optional("df_note"), get_optional("sd_note"), get_optional("outlier_note"), get_optional("slice_outlier_note"),
    "Linear mixed analyses follow the supplied reference scripts. Control comparisons use the emmeans Dunnett approximation (dunnettx); Holm p-values use Bonferroni intervals. Candidate variance tests are conventional likelihood-ratio approximations.",
    "Cells share mouse and slice effects. Report mice / slices / cells separately. Missing measurements are omitted individually; no values are imputed.")
  list(method = config$method, config = config, tables = tables, formatted_tables = formatted, data = data.frame(y = long$Response, group = long$Condition),
    column_data = data.frame(y = long$Value, group = long$Condition), fit = env$model, nested = TRUE,
    graphs = graphs, publication_dimensions = list(width = settings$pub_width_mm, height = env$pub_total_height_mm), artifacts = artifacts, console = console, report_markdown = get_optional("note_lines"), notes = notes,
    warnings = unique(c(warnings, get_optional("fit_warnings"), get_optional("pub_save_errors"))))
}

pack_nested_config <- function(cfg) {
  # JSON named lists preserve label/colour keys that named atomic vectors can lose.
  if (is.null(cfg) || !cfg$method %in% nested_methods) return(cfg)
  for (n in c("pub_key_labels", "pub_condition_labels", "pub_mouse_colours")) {
    v <- cfg$nested[[n]]
    if (!is.null(v) && !is.list(v)) cfg$nested[[n]] <- as.list(v)
  }
  cfg
}

nested_demo <- function(between = FALSE) {
  # A fixed seed makes the hierarchy and one intentionally blank measurement
  # reproducible for examples, regression checks and UI troubleshooting.
  set.seed(109)
  d <- expand.grid(Cell = paste0("C", 1:6), Slice = paste0("S", 1:3), Mouse = paste0("M", 1:9), stringsAsFactors = FALSE)
  mouse <- rep(rnorm(9, sd = 1.4), each = 18); slice <- rep(rnorm(27, sd = .9), each = 6)
  baseline <- 12 + mouse + slice + rnorm(nrow(d), sd = .8)
  if (between) {
    d$Condition <- rep(c("Control", "Drug", "Wash"), each = 54)
    d$Value <- baseline + rep(c(0, 2, 1), each = 54); d$Value[5] <- NA_real_
  } else {
    d$Baseline2 <- baseline
    d$Drug <- baseline + 2 + rep(rnorm(27, sd = 1.1), each = 6) + rnorm(nrow(d), sd = .7)
    d$Wash <- baseline + .8 + rep(rnorm(27, sd = .6), each = 6) + rnorm(nrow(d), sd = .7)
    d$Drug[5] <- NA_real_
  }
  d[c("Mouse", "Slice", "Cell", setdiff(names(d), c("Mouse", "Slice", "Cell")))]
}

draw_nested_graph <- function(plot) {
  if (inherits(plot, "grob")) { grid::grid.newpage(); grid::grid.draw(plot) } else print(plot)
}

nested_normality <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  sw <- if (n >= 3 && n <= 5000 && length(unique(x)) > 1) shapiro.test(x) else NULL
  da <- dagostino_pearson(x)
  data.frame(n = n, SW_W = if (is.null(sw)) NA_real_ else unname(sw$statistic),
    SW_p = if (is.null(sw)) NA_real_ else sw$p.value, DAgo_K2 = unname(da["K2"]), DAgo_p = unname(da["p_value"]))
}

nested_html_report <- function(result) {
  # Turn Obsidian figure links into self-contained browser images. Base64 must be
  # one unbroken URI; wrapped encoder output is interpreted as text by Markdown.
  text <- paste(result$report_markdown, collapse = "\n")
  # Obsidian metadata belongs in the Markdown note, not the visible HTML report.
  text <- sub("(?s)^---\\r?\\n.*?\\r?\\n---(?:\\r?\\n|$)", "", text, perl = TRUE)
  text <- gsub("!\\[\\[([^]]+)\\]\\]", "![Figure](\\1)", text)
  for (n in names(result$artifacts)[grepl("\\.png$", names(result$artifacts))]) {
    encoded <- gsub("[[:space:]]", "", jsonlite::base64_enc(result$artifacts[[n]]))
    text <- gsub(paste0("(", n, ")"), paste0("(data:image/png;base64,", encoded, ")"), text, fixed = TRUE)
  }
  text <- gsub("\\[\\[([^]]+)\\]\\]", "\\1", text)
  # Keep data labels as text rather than active HTML inside the report frame.
  text <- gsub("&", "&amp;", text, fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE); text <- gsub(">", "&gt;", text, fixed = TRUE)
  body <- commonmark::markdown_html(text, extensions = TRUE)
  paste0('<!doctype html><html><head><meta charset="utf-8"><title>Linear mixed analysis report</title><style>body{font:15px system-ui;max-width:1100px;margin:30px auto;padding:0 24px;color:#25364a;background:white}[hidden]{display:none!important}img{max-width:100%}table{border-collapse:collapse;display:block;overflow:auto}td,th{border:1px solid #ddd;padding:6px}button,select{padding:8px;margin:4px}.report-section{margin:24px 0}</style></head><body>',
    '<nav><button onclick="move(-1)">Previous</button><button onclick="move(1)">Next</button><select id="sections" onchange="show(this.value)"></select><button onclick="allSections()">Show all</button></nav>',
    body, '<script>const nodes=[...document.body.children].filter(n=>n.tagName!=="NAV"&&n.tagName!=="SCRIPT");let groups=[],current=[];for(const n of nodes){if(/^H[12]$/.test(n.tagName)&&current.length){groups.push(current);current=[]}current.push(n)}if(current.length)groups.push(current);let index=0;const s=document.getElementById("sections");groups.forEach((g,i)=>s.add(new Option(g.find(n=>/^H[12]$/.test(n.tagName))?.textContent||"Introduction",i)));function show(i){index=Math.max(0,Math.min(groups.length-1,+i));groups.forEach((g,j)=>g.forEach(n=>n.hidden=j!==index));s.value=index}function move(d){show(index+d)}function allSections(){groups.forEach(g=>g.forEach(n=>n.hidden=false))}document.addEventListener("keydown",e=>{if(e.target.matches("input,textarea,select"))return;if(e.key==="ArrowRight")move(1);if(e.key==="ArrowLeft")move(-1)});show(0);</script></body></html>')
}

write_nested_bundle <- function(result, file) {
  file <- file.path(normalizePath(dirname(file), mustWork = TRUE), basename(file))
  work <- tempfile("labstats-bundle-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  for (n in names(result$artifacts)) writeBin(result$artifacts[[n]], file.path(work, basename(n)))
  for (n in names(result$tables)) write.csv(result$tables[[n]], file.path(work, paste0(n, ".csv")), row.names = FALSE)
  writeLines(nested_html_report(result), file.path(work, "report.html"), useBytes = TRUE)
  zip::zipr(file, list.files(work), root = work)
}

save_result_figure <- function(result, kind, plot, file, format) {
  if (isTRUE(result$nested) && identical(kind, "publication")) {
    artifact <- result$artifacts[[paste0("nested_analysis_figure.", format)]]
    if (!is.null(artifact)) { writeBin(artifact, file); return(invisible(file)) }
    s <- result$config$nested
    ggplot2::ggsave(file, plot, width = s$pub_width_mm, height = result$publication_dimensions$height, units = "mm", device = if (format == "pdf") grDevices::cairo_pdf else format, dpi = 600, bg = "white")
  } else ggplot2::ggsave(file, plot, width = 7, height = 5, device = if (format == "pdf") grDevices::cairo_pdf else format, dpi = 300, bg = "white")
}
