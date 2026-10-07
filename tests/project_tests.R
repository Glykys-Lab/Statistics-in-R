source("R/project.R")

# Published SciPy skew/kurtosis reference values on a nonsymmetric sample.
dp <- dagostino_pearson(c(2, 8, 0, 4, 1, 9, 9, 0))
check(near(dp["skew_z"], .44626385374196975) && near(dp["kurtosis_z"], -1.850682505457718), "D'Agostino-Pearson transformations match numerical reference")
check(near(dp["p_value"], pchisq(sum(dp[c("skew_z", "kurtosis_z")]^2), 2, lower.tail = FALSE)), "Omnibus normality uses chi-square with two df")
check(near(dagostino_pearson(0:19)["kurtosis_z"], -1.7058104152122062), "Kurtosis transformation matches published 20-value reference")
check(near(dagostino_pearson(c(2, 8, 0, 4, 1, 9, 9, 0)*3 + 90), dp), "Normality test is invariant to location and positive scale")
check(all(is.na(dagostino_pearson(1:7))) && all(is.na(dagostino_pearson(rep(2, 10)))), "Normality sample-size and constant-data limits are explicit")
nd <- run_analysis(d, c(cfg, method = "normality"))
log_sw <- nd$tables$distribution_tests
log_sw <- log_sw[log_sw$group == "Control" & log_sw$distribution == "Lognormal" & log_sw$test == "Shapiro-Wilk", ]
check(near(log_sw$p_value, shapiro.test(log(d$response[d$group == "Control"]))$p.value), "Lognormal Shapiro test uses logged column values")
nonpositive <- d; nonpositive$response[1] <- 0
z <- run_analysis(nonpositive, c(cfg, method = "normality"))$tables$distribution_tests
check(all(is.na(z$p_value[z$group == "Control" & z$distribution == "Lognormal"])), "Zero values make lognormal tests unavailable without dropping them")

sidak <- run_analysis(four, c(rm_cfg, adjustment = "sidak", comparisons = "control", control = "A"))
post <- sidak$tables$posthoc
reference_pair <- t.test(scores[, 1], scores[, 2], paired = TRUE, conf.level = .95^(1/3))
check(nrow(post) == 3 && all(post$group_1 == "A"), "Control comparisons form exactly the selected family")
check(near(post$adjusted_p_value, 1 - (1-post$p_value)^3), "Sidak p-values use the selected family size")
check(near(post$ci_low[1], reference_pair$conf.int[1]) && near(post$ci_high[1], reference_pair$conf.int[2]), "Sidak intervals match independently adjusted paired t intervals")
ordinary <- run_analysis(d, c(cfg, method = "anova", adjustment = "sidak"))
po <- ordinary$tables$posthoc; mse <- summary(aov(response ~ group, d))[[1]]["Residuals", "Mean Sq"]
check(near(po$se[1], sqrt(mse*(1/8 + 1/8))), "Ordinary Sidak comparisons use ANOVA pooled error")
check(fails(run_analysis(d, c(cfg, method = "anova", adjustment = "tukey", comparisons = "control", control = "Control")), "all pairwise"), "Tukey cannot silently replace a selected control family")
log_result <- run_analysis(four |> transform(response = exp(response/10)), modifyList(rm_cfg, list(method = "log_rm_anova", adjustment = "sidak")))
log_reference <- run_analysis(four |> transform(response = response/10), c(rm_cfg, adjustment = "sidak"))
check(near(log_result$tables$test$p_value, log_reference$tables$test$p_value), "Lognormal repeated ANOVA equals ANOVA on logged responses")
check(near(log_result$tables$posthoc$geometric_mean_ratio, exp(log_reference$tables$posthoc$estimate)), "Lognormal comparisons report geometric-mean ratios")
check(fails(run_analysis(nonpositive, modifyList(rm_cfg, list(method = "log_rm_anova"))), "positive"), "Lognormal models reject nonpositive observed outcomes")
check(unname(design_methods(matched = TRUE, distribution = "lognormal", repeated = "mixed")) == "log_mixed", "Design routes lognormal missing-data analysis")

gg_mixed <- run_analysis(four, c(mixed_cfg, adjustment = "sidak"))
check(near(gg_mixed$tables$test$epsilon, expected_epsilon), "Complete mixed-model GG uses the complete covariance epsilon")
check(near(gg_mixed$tables$test$p_value, pf(gg_mixed$tables$test$statistic, 3*expected_epsilon, gg_mixed$tables$test$denominator_df, lower.tail = FALSE)), "Mixed GG applies epsilon to both F-test degrees of freedom")
independent_fit <- nlme::gls(y ~ group, data = gg_mixed$data, method = "REML")
matching <- gg_mixed$tables$matching
check(near(matching$statistic[1], 2*(as.numeric(logLik(gg_mixed$fit))-as.numeric(logLik(independent_fit)))), "Matching-effectiveness LRT compares identical fixed-effect REML fits")
check(near(matching$p_value[1], pchisq(matching$statistic[1], 1, lower.tail = FALSE)), "Conventional matching test uses chi-square one df")
check(near(gg_epsilon(em_covariance(scores)$covariance), expected_epsilon, 1e-6), "Gaussian EM reproduces complete-data covariance epsilon")
check(mixed_result$tables$test$epsilon >= .5 && mixed_result$tables$test$epsilon <= 1 && mixed_result$tables$data_audit$used_rows == 22, "Missing-data GG remains bounded and retains original observations")

outlier_input <- data.frame(response = c(1:20, 150), group = "A")
outlier_result <- run_analysis(outlier_input, c(cfg, method = "outliers", rout_q = .01))
check(identical(which(outlier_result$tables$outlier_details$outlier), 21L), "ROUT detects a known gross outlier")
check(is.na(outlier_result$tables$cleaned_data$response[21]) && outlier_input$response[21] == 150, "ROUT cleaned sheet blanks outlier while preserving source")
check(identical(rout_column(outlier_input$response*2+50, .01)$outlier, rout_column(outlier_input$response, .01)$outlier), "ROUT flags are invariant to positive scale and translation")
check(!any(rout_column(rep(4, 12), .01)$outlier), "Constant columns have no estimable ROUT scale")
check(identical(significance_summary(c(.2, .02, .002, .0002, .00001)), c("ns", "*", "**", "***", "****")), "Significance stars use conventional thresholds")

wide_transform <- data.frame(Subject = c("s1", "s2", "s3"), Baseline = c(2, 0, NA), After = c(6, 8, 9))
tr <- normalize_wide(wide_transform, c("Baseline", "After"), "Baseline", "ratio")
check(tr$data$After[1] == 3 && is.na(tr$data$After[2]) && is.na(tr$data$After[3]) && identical(tr$data$Subject, wide_transform$Subject), "Ratio tables preserve matching and blank undefined baseline divisions")
check(near(normalize_wide(wide_transform, c("Baseline", "After"), "Baseline", "change")$data$After[1:2], c(4, 8)), "Change-from-baseline keeps a zero baseline valid")
check(normalize_wide(wide_transform, c("Baseline", "After"), "Baseline", "percent_change")$data$After[1] == 200, "Percent-change normalization uses the correct denominator")
check(is.na(normalize_wide(wide_transform, c("Baseline", "After"), "Baseline", "log")$data$Baseline[2]), "Log transformation blanks nonpositive values in derived table")

p <- new_project("Experiment 1"); p$info <- "Notes: \u03b1 = 0.05"
add <- project_add(p, "tables", list(name = "Repeated data", type = "column", layout = "wide", data = column_demo(), columns = c("Control", "Drug A", "Drug B"), id_column = "Subject"))
p <- add$project; tid <- add$id
add <- project_add(p, "analyses", list(name = "RM result", table_id = tid, worksheet = p$tables[[tid]], input_data = d, result = rm_result, config = rm_cfg, time = "2026-10-05"))
p <- add$project; aid <- add$id
add <- project_add(p, "graphs", list(name = "Intervals", analysis_id = aid, kind = "comparison", options = list(title = "Confidence intervals")))
p <- add$project; gid <- add$id
p <- project_add(p, "layouts", list(name = "Layout 1", graph_ids = gid, columns = 1L, mode = "grid"))$project
restored <- read_project(project_json(p))
check(identical(restored$tables[[tid]]$data, p$tables[[tid]]$data) && identical(restored$info, p$info), "JSON project round-trip preserves data types, names, missing cells and notes")
check(near(restored$analyses[[aid]]$result$tables$test$p_value, rm_result$tables$test$p_value), "Reopening project recalculates immutable analysis snapshots")
check(length(restored$graphs) == 1 && length(restored$layouts) == 1, "Projects retain saved graph and layout references")
p_bad <- p; p_bad$graphs[[gid]]$analysis_id <- "absent"
check(fails(read_project(project_json(p_bad)), "reference"), "Project import rejects broken graph references before replacing workspace")
check(identical(unpack_frame(pack_frame(data.frame(X = c(1, NA), note = c("a", NA), flag = c(TRUE, NA))))$flag, c(TRUE, NA)), "JSON round-trip preserves logical missing values")

for (kind in c("raw", "normal_qq", "lognormal_qq", "comparison", "residual", "assumption_qq")) ggplot2::ggplot_build(make_graph(if (kind == "lognormal_qq") log_result else sidak, kind))
ggplot2::ggplot_build(make_plot(sidak, varied_symbols = TRUE, palette_name = "prism", symbol = 1, x_angle = 45))
ggplot2::ggplot_build(make_graph(log_result, "comparison"))
ggplot2::ggplot_build(overlay_graphs(list(list(name = "A", kind = "raw", result = sidak), list(name = "B", kind = "raw", result = sidak)), "Overlay"))
grDevices::pdf("artifacts/panel_layout.pdf", width = 10, height = 7)
draw_graph_layout(list(make_graph(sidak), make_graph(sidak, "comparison")), 2, "Panel layout"); grDevices::dev.off()
check(TRUE, "New QQ, adjusted CI, geometric ratio, formatted symbols, overlay and panel figures render")
writeLines(project_json(p), "artifacts/example_project.json")
