# Verify that visible t-test controls become the intended saved configuration.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
if (!exists("app_env")) {
  app_env <- new.env(parent = globalenv())
  sys.source("app.R", envir = app_env)
}
shiny::testServer(app_env$server, {
  session$setInputs(table_type = "column", source = "demo", demo = "groups", layout = "wide",
    wide_columns = c("Control", "Drug A"), wide_id = "Subject", family = "two", matching = "independent",
    distribution = "normal", variance = "unequal", method = "welch", y = "response", group = "group",
    confidence = ".95", alpha = .05, decimal = ".", t_alternative = "greater", t_mu = 0,
    report_descriptive = TRUE, report_adjusted_p = TRUE, plot_style = "points", error_style = "ci",
    plot_title = "", x_label = "", y_label = "")
  session$setInputs(run = 1)
  if (!is.null(analysis_error())) stop(analysis_error())
  z <- split(data()$response, data()$group)
  stopifnot(result()$config$alternative == "greater", is_current(), analysis_unchanged(),
    isTRUE(all.equal(result()$tables$test$p_value, t.test(z[[1]], z[[2]], alternative = "greater")$p.value)))
  stopifnot(grepl("Control - Drug A > 0", output$t_hypothesis_help, fixed = TRUE))
  session$setInputs(t_ratio = 3)
  stopifnot(analysis_unchanged())
  session$setInputs(t_alternative = "less")
  stopifnot(!is_current(), !analysis_unchanged())
  session$setInputs(run = 2)
  stopifnot(is_current(), result()$config$alternative == "less")
  session$setInputs(t_mu = -1)
  stopifnot(!analysis_unchanged())
  session$setInputs(run = 3)
  stopifnot(result()$config$null_value == -1)
  session$setInputs(family = "one", method = "one_sample", one_sample_group = "Drug A", t_mu = 10)
  session$setInputs(run = 4)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(result()$tables$test$group_1 == "Drug A", result()$tables$test$n1 == 8,
    result()$config$null_value == 10, result()$config$alternative == "less")
  session$setInputs(family = "two", method = "log_paired_t", distribution = "lognormal", matching = "matched", id = "subject", t_ratio = 1)
  session$setInputs(run = 5)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(result()$config$scale == "log", result()$tables$test$n_pairs == 8,
    "geometric_mean_ratio" %in% names(result()$tables$test), result()$config$null_ratio == 1,
    !"null_value" %in% names(current_config()))

  session$setInputs(demo = "multiple_tests", family = "multiple_t", method = "multiple_t", t_test_mode = "endpoints",
    layout = "long", matching = "independent", distribution = "normal", variance = "unequal",
    test_column = "endpoint", t_alternative = "two.sided", t_mu = 0,
    error_control = "fwer", adjustment = "holm", graph_comparison = TRUE)
  session$setInputs(run = 6)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(nrow(data()) == 36, result()$config$t_test_mode == "endpoints", result()$config$test_column == "endpoint",
    nrow(result()$tables$posthoc) == 3, all(result()$tables$posthoc$family_size == 3),
    result()$config$t_test_kind == "welch", analysis_unchanged())
  session$setInputs(error_control = "fdr", adjustment = "BH", fdr_q = .1)
  stopifnot(!analysis_unchanged())
  session$setInputs(run = 7)
  stopifnot("discovery" %in% names(result()$tables$posthoc), result()$config$fdr_q == .1)
  session$setInputs(method = "log_multiple_t", distribution = "lognormal", matching = "matched", t_alternative = "less", id = "subject")
  session$setInputs(run = 8)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(all(result()$tables$posthoc$n_pairs == 6), result()$config$t_test_kind == "paired_t", result()$config$scale == "log")
  archived <- saved()$id
  session$setInputs(demo = "groups", layout = "wide", wide_columns = c("Control", "Drug A", "Drug B"),
    method = "multiple_t", distribution = "normal", matching = "independent", variance = "equal", t_test_mode = "columns",
    comparisons = "selected", selected_pairs = "Control::Drug%20B", error_control = "fwer", adjustment = "bonferroni",
    t_alternative = "greater")
  session$setInputs(run = 9)
  if (!is.null(analysis_error())) stop(analysis_error())
  stopifnot(result()$config$t_test_mode == "columns", result()$config$t_test_kind == "student",
    nrow(result()$tables$posthoc) == 1, result()$tables$posthoc$group_2 == "Drug B")
  open_analysis(archived)
  stopifnot(is_current(), result()$config$t_test_mode == "endpoints", result()$config$alternative == "less")
  restore_preferences(list(t_alternative = "less", t_test_mode = "columns", t_mu = 1234, test_column = "private", t_ratio = 5))
  stopifnot(setequal(names(preferences()), c("t_alternative", "t_test_mode")))
  cat("PASS: T-test controls preserve tails, hypotheses, sample selection, lognormal pairing, both multiple-test workflows and archived configurations\n")
})
