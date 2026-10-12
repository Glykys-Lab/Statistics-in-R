# Focused server checks sourced by server_tests.R; app_env supplies the server.
# Follow advanced UI choices through configuration, automatic method selection,
# diagnostic graph creation, FDR comparisons and validation of restored defaults.
shiny::testServer(app_env$server, {
  session$setInputs(table_type = "column", source = "demo", demo = "groups", layout = "wide",
    wide_columns = c("Drug B", "Control", "Drug A"), wide_id = "Subject", family = "many",
    matching = "matched", distribution = "normal", variance = "unequal", repeated = "auto",
    method = "rm_auto", y = "response", group = "group", id = "subject", confidence = ".95",
    alpha = .05, adjustment = "sidak", error_control = "fwer", comparisons = "all", sphericity = "gg",
    decimal = ".", random_effect = "remove", glm_initial = TRUE, rank_test = "dunn",
    report_descriptive = TRUE, report_adjusted_p = TRUE, graph_comparison = TRUE,
    report_goodness = TRUE, report_aicc = TRUE, report_effects = TRUE,
    residual_normality = TRUE, graph_residual = TRUE, graph_homoscedasticity = TRUE,
    graph_residual_qq = TRUE, graph_residual_heatmap = TRUE, graph_differences = TRUE)
  session$setInputs(run = 1)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(result()$method == "rm_anova", analysis_unchanged(),
    identical(levels(result()$data$group), c("Control", "Drug A", "Drug B")))
  saved_kinds <- vapply(project()$graphs, `[[`, character(1), "kind")
  stopifnot(all(c("comparison", "residual", "homoscedasticity", "residual_qq", "residual_heatmap", "differences") %in% saved_kinds),
    all(c("goodness_of_fit", "effect_sizes", "model_comparison", "residual_normality") %in% names(display_tables())))
  # Formatting must stay live without making the saved calculation stale.
  session$setInputs(output_digits = 3, p_style = "decimal", p_digits = 4)
  stopifnot(analysis_unchanged(), display_options()$output_digits == 3)
  count <- length(project()$analyses)
  session$setInputs(run = 11)
  stopifnot(length(project()$analyses) == count)
  pair <- paste(vapply(c("Control", "Drug B"), utils::URLencode, character(1), reserved = TRUE), collapse = "::")
  session$setInputs(comparisons = "selected", selected_pairs = pair, error_control = "fdr", adjustment = "BKY", fdr_q = .1,
    reverse_comparisons = TRUE, report_descriptive = FALSE, report_adjusted_p = FALSE)
  stopifnot(!analysis_unchanged())
  session$setInputs(run = 2)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(nrow(result()$tables$posthoc) == 1, result()$tables$posthoc$group_1 == "Drug B",
    result()$config$adjustment == "BKY", "discovery" %in% names(display_tables()$posthoc),
    !"summary" %in% names(display_tables()), !"adjusted_p_value" %in% names(display_tables()$posthoc))
  session$setInputs(demo = "missing", method = "rm_auto", repeated = "auto", comparisons = "none", report_descriptive = TRUE)
  session$setInputs(run = 3)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(result()$method == "mixed", result()$tables$data_audit$used_rows == 22,
    is.null(result()$tables$posthoc), analysis_unchanged())
  id <- saved()$id
  stopifnot(!any(vapply(Filter(function(g) g$analysis_id == id, project()$graphs), `[[`, character(1), "kind") == "comparison"))
  # Browser preferences cannot inject data mappings, invalid enums or numbers.
  restore_preferences(list(random_effect = "full", repeated = "mixed", matching = "matched", alpha = .05,
    confidence = .95, output_digits = 50, p_style = "invalid", y = "private-data-column", glm_initial = TRUE))
  stopifnot(setequal(names(preferences()), c("random_effect", "repeated", "matching", "alpha", "confidence", "glm_initial")))
  cat("PASS: Advanced server options, diagnostic graph creation, display formatting, selected FDR contrasts, automatic missing-data selection and preference validation\n")
})
