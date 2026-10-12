# Nested-model regressions; run standalone from the project root or through run_tests.R.
# Fit deterministic Mouse/Slice/Cell examples and compare estimates/contrasts to
# directly specified models. Check invalid IDs/pairs, omissions, ROUT, report images,
# physical-size figure exports, project settings and standalone R reproducibility.
# Files produced by these checks stay under artifacts/nested-tests/.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
source("R/analysis.R"); source("R/data.R"); source("R/plots.R"); source("R/project.R")
dir.create("artifacts/nested-tests", recursive = TRUE, showWarnings = FALSE)
nested_checks <- 0L
expect_nested <- function(ok, label) { stopifnot(isTRUE(ok)); nested_checks <<- nested_checks + 1L; cat("PASS:", label, "\n") }
expect_bad_nested <- function(expr, text) {
  e <- tryCatch(force(expr), error = identity)
  expect_nested(inherits(e, "error") && grepl(text, conditionMessage(e), fixed = TRUE), text)
}
d <- nested_demo()
cfg <- list(method = "nested_repeated", mouse = "Mouse", slice = "Slice", cell = "Cell", columns = c("Baseline2", "Drug", "Wash"),
  nested = list(model_engine = "lmer", slice_condition = "shared", outlier_action = "none", pub_formats = c("svg", "pdf")))
r <- run_analysis(d, cfg)
expect_nested(r$tables$data_audit$used_rows == nrow(d)*3L - 1L && r$tables$data_audit$missing_values == 1L,
  "Nested repeated models omit only the missing measurement")
z <- r$tables$analyzed_measurements; z$Condition <- factor(z$Condition, levels = c("Baseline2", "Drug", "Wash")); z$Response <- z$Value
reference <- lmerTest::lmer(Response ~ Condition + (1 | Mouse/Slice/Cell) + (1 | Mouse:Slice:Condition), data = z, control = lme4::lmerControl(optimizer = "bobyqa"))
reference_means <- emmeans::emmeans(reference, ~Condition, lmer.df = "satterthwaite")
reference_pairs <- as.data.frame(emmeans::contrast(reference_means, "revpairwise", adjust = "tukey"))
expect_nested(max(abs(lme4::fixef(r$fit) - lme4::fixef(reference))) < 1e-7 &&
  max(abs(r$tables$multiple_comparisons$p_adj - reference_pairs$p.value)) < 1e-7,
  "Repeated estimates and Tukey probabilities match a directly fitted nested model")
expect_nested(all(c("Mouse", "Slice:Mouse", "Cell:Slice:Mouse", "Mouse:Slice:Condition") %in% names(lme4::getME(r$fit, "flist"))), "Mouse, slice, cell and slice x condition terms are fitted")
expect_nested(all(c("analysis_decisions", "candidate_models", "degrees_of_freedom", "random_effects", "writing_template") %in% names(r$tables)), "Decision, candidate model, df source, variance and writing tables are present")
expect_nested(all(is.finite(r$tables$residual_normality$DAgo_p)), "Shared D'Agostino-Pearson helper supplies residual normality")
expect_nested(identical(as.character(r$tables$model_means$Condition), c("Baseline2", "Drug", "Wash")), "Nested tables preserve condition order")
expect_nested(all(c("nested_analysis_figure.svg", "nested_analysis_figure.pdf") %in% names(r$artifacts)) && inherits(make_graph(r, "publication"), "grob"), "Publication SuperPlot and editable SVG / vector PDF are generated")
html <- nested_html_report(r)
expect_nested(grepl("data:image/png;base64,", html, fixed = TRUE) && grepl("Previous", html) && grepl("Show all", html), "HTML report embeds figures and supports section navigation")
expect_nested(grepl('<img src="data:image/png;base64,', html, fixed = TRUE) && !grepl('![Figure](data:', html, fixed = TRUE), "Report figures render as image elements instead of wrapped base64 text")
expect_nested(!grepl('<p>tags:', html, fixed = TRUE), "Obsidian metadata is excluded from the visible HTML report")
decision_frame <- sheet_frame(result_sheets(r)$analysis_decisions)
expect_nested(identical(names(decision_frame), names(r$formatted_tables$analysis_decisions)) && nzchar(decision_frame[1, 1]), "Nested result sheets and arranged exports have no empty statistic column")
write_nested_bundle(r, "artifacts/nested-tests/report.zip")
bundle <- zip::zip_list("artifacts/nested-tests/report.zip")$filename
expect_nested(all(c("report.html", "nested_analysis.md", "model_means.csv", "nested_analysis_figure.svg") %in% bundle), "Portable ZIP contains report, Obsidian note, numeric tables and figures")
bad <- d; bad$Drug[1] <- 0
expect_bad_nested(run_analysis(bad, modifyList(cfg, list(nested = list(log_transform = TRUE)))), "every observed value to be positive")
bad <- d; bad$Mouse[1] <- NA
expect_bad_nested(run_analysis(bad, cfg), "IDs cannot be blank")
bad <- rbind(d, d[1, ])
expect_bad_nested(run_analysis(bad, cfg), "must occupy one row")
bad <- d; bad$Drug <- NA_real_
expect_bad_nested(run_analysis(bad, cfg), "Every selected condition")
expect_bad_nested(run_analysis(d, modifyList(cfg, list(nested = list(comparison_type = "custom", custom_pairs = list(c("Drug", "Drug")))))), "two different observed conditions")
expect_bad_nested(run_analysis(d, modifyList(cfg, list(nested = list(comparison_type = "custom", custom_pairs = list(c("Drug", "Wash"), c("Wash", "Drug")))))), "only one direction")
control_cfg <- modifyList(cfg, list(control = "Drug", nested = list(comparison_type = "vs_first", comparison_plot = "none", pub_formats = "svg", show_diagnostics = FALSE, check_variance_homogeneity = FALSE)))
control_result <- run_analysis(d, control_cfg)
control_reference <- as.data.frame(emmeans::contrast(reference_means, "trt.vs.ctrl", ref = 2, adjust = "dunnettx"))
expect_nested(identical(as.character(control_result$tables$model_means$Condition), c("Baseline2", "Drug", "Wash")) &&
  max(abs(control_result$tables$multiple_comparisons$p_adj - control_reference$p.value)) < 1e-7,
  "Choosing a middle control preserves graph order and matches direct Dunnett contrasts")
expect_nested(is.null(control_result$graphs$comparison) && control_result$publication_dimensions$height == 70,
  "Publication dimensions remain correct when the notebook comparison plot is disabled")
save_result_figure(control_result, "publication", make_graph(control_result, "publication"), "artifacts/nested-tests/publication.pdf", "pdf")
expect_nested(file.info("artifacts/nested-tests/publication.pdf")$size > 1000, "PDF export works when only SVG was requested in the publication bundle")
expect_nested(is.null(control_result$graphs$residual) && is.null(control_result$tables$cell_variances), "Diagnostic and variance-report switches hide the requested outputs")
cfg_log <- modifyList(cfg, list(nested = list(log_transform = TRUE, comparison_type = "custom", custom_pairs = list(c("Drug", "Baseline2")), p_adjust = "holm", df_method = "smallest", compare_to_zero = TRUE, make_pub_figures = FALSE)))
logged <- run_analysis(d, cfg_log)
expect_nested(nrow(logged$tables$multiple_comparisons) == 1L && "Ratio" %in% names(logged$tables$multiple_comparisons) &&
  all(logged$tables$model_means$Geometric_mean > 0), "Log10 geometric means and selected ratio comparison are returned")
expect_nested(all(logged$tables$mean_tests$null_geometric_mean == 1) && any(grepl("against 1", logged$warnings)), "Log-scale zero tests explicitly report their original-scale null of one")
b <- nested_demo(TRUE)
cfg_between <- list(method = "nested_between", mouse = "Mouse", slice = "Slice", cell = "Cell", condition = "Condition", value = "Value",
  nested = list(model_engine = "nlme", slice_variance = "per_condition", comparison_type = "vs_first", df_method = "smallest", outlier_action = "none", make_pub_figures = FALSE))
between <- run_analysis(b, cfg_between)
expect_nested(inherits(between$fit, "lme") && !is.null(between$fit$modelStruct$varStruct) && nrow(between$tables$multiple_comparisons) == 2L,
  "Between-group nlme fits per-condition slice and residual SDs and control comparisons")
expect_nested(length(unique(between$tables$multiple_comparisons$DF)) == 1L && max(between$tables$multiple_comparisons$DF) < 30,
  "Between-group df come from the lmer reference and can be shared across comparisons")
bad <- b; bad$Condition[1] <- "Drug"
expect_bad_nested(run_analysis(bad, cfg_between), "slice appears under several conditions")
cfg_auto <- modifyList(cfg, list(nested = list(model_engine = "auto", slice_condition = "auto", comparison_type = "consecutive", outlier_action = "none", make_pub_figures = FALSE)))
automatic <- run_analysis(d, cfg_auto)
expect_nested(nrow(automatic$tables$multiple_comparisons) == 2L && sum(automatic$tables$candidate_models$Used == "<- used") == 1L && nrow(automatic$tables$residual_model_tests) > 0,
  "Automatic candidate selection and consecutive comparisons are supported")
# A clearly extreme cell is audited by raw ROUT; slice flags remain in the data.
outlier_data <- d; outlier_data$Drug[1] <- 1000
flag_cfg <- modifyList(cfg, list(nested = list(outlier_action = "flag", outlier_unit = "both", outlier_level = "raw", make_pub_figures = FALSE)))
flagged <- run_analysis(outlier_data, flag_cfg)
removed <- run_analysis(outlier_data, modifyList(flag_cfg, list(nested = list(outlier_action = "remove"))))
expect_nested(nrow(flagged$tables$cell_outliers) >= 1L && flagged$tables$data_audit$removed_outliers == 0L && removed$tables$data_audit$removed_outliers >= 1L && outlier_data$Drug[1] == 1000,
  "ROUT flags cells and optionally excludes them from the model without editing the source")
saved_cfg <- r$config
saved_cfg$nested$pub_condition_labels <- c(Baseline2 = "Baseline", Drug = "Drug", Wash = "Wash")
saved_cfg$nested$pub_mouse_colours <- c(M1 = "#0072B2")
p <- new_project(); t <- project_add(p, "tables", list(name = "Nested", type = "nested", layout = "wide", columns = cfg$columns, id_column = "", data = d, nested_config = saved_cfg)); p <- t$project
a <- project_add(p, "analyses", list(name = "Repeated nested", table_id = t$id, worksheet = p$tables[[t$id]], input_data = d, config = saved_cfg, result = r, time = "2026-10-09")); p <- a$project
reopened <- read_project(project_json(p))
expect_nested(max(abs(reopened$analyses[[a$id]]$result$tables$multiple_comparisons$p_adj - r$tables$multiple_comparisons$p_adj)) < 1e-10,
  "Nested settings and results survive project JSON recalculation")
expect_nested(identical(unlist(reopened$analyses[[a$id]]$config$nested$pub_key_labels), nested_defaults()$pub_key_labels) &&
  identical(names(reopened$analyses[[a$id]]$config$nested$pub_mouse_colours), "M1") &&
  identical(names(reopened$analyses[[a$id]]$config$nested$pub_condition_labels), c("Baseline2", "Drug", "Wash")), "Named key labels, condition labels and mouse colours survive JSON")
script <- reproduction_script(d, cfg, graph_kind = "publication", worksheet = p$tables[[t$id]])
writeLines(script, "artifacts/nested-tests/reproduce.R", useBytes = TRUE)
path <- normalizePath("artifacts/nested-tests/reproduce.R"); old <- getwd(); setwd("artifacts/nested-tests")
env <- new.env(parent = globalenv()); source(path, local = env, encoding = "UTF-8"); setwd(old)
expect_nested(max(abs(env$result$tables$multiple_comparisons$p_adj - r$tables$multiple_comparisons$p_adj)) < 1e-10 && file.exists("artifacts/nested-tests/nested_analysis.zip"),
  "Exported nested script runs independently of the project directory and reproduces results")
cat(nested_checks, "nested-model checks passed.\n")
