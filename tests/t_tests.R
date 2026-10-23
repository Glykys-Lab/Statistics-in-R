# Compare t-test workflows against independent calls to stats::t.test/p.adjust.
# Run directly or through run_tests.R; include persistence and standalone exports.
if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
if (!exists("t_test_methods")) source("R/analysis.R")
source("R/data.R"); source("R/plots.R"); source("R/project.R")
if (!exists("checks")) checks <- 0L
if (!exists("check")) check <- function(ok, description) {
  if (!isTRUE(ok)) stop("FAILED: ", description)
  checks <<- checks + 1L; cat("PASS:", description, "\n")
}
if (!exists("near")) near <- function(x, y, tolerance = 1e-9) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
if (!exists("fails")) fails <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  !is.null(error) && grepl(pattern, conditionMessage(error), ignore.case = TRUE)
}
local({
  x <- c(8, 10, 9, 12, 11, 13); y <- c(11, 12, 14, 13, 15, 14)
  d <- data.frame(subject = rep(paste0("S", 1:6), 2), group = rep(c("A", "B"), each = 6), response = c(x, y))
  base <- list(y = "response", group = "group", id = "subject")
  for (method in c("welch", "student", "paired_t", "log_welch", "log_student", "log_paired_t")) {
    lognormal <- startsWith(method, "log_")
    a <- if (lognormal) log(x) else x; b <- if (lognormal) log(y) else y
    for (alternative in c("two.sided", "greater", "less")) {
      cfg <- modifyList(base, list(method = method, alternative = alternative))
      # Explicit factor order makes the reference contrast A minus B after shuffling.
      shuffled <- d[c(9, 1, 12, 6, 8, 4, 10, 2, 11, 5, 7, 3), ]; shuffled$group <- factor(shuffled$group, levels = c("A", "B"))
      result <- run_analysis(shuffled, cfg)
      ref <- stats::t.test(a, b, paired = grepl("paired_t", method), var.equal = grepl("student", method), alternative = alternative)
      t <- result$tables$test
      check(near(t$p_value, ref$p.value) && near(t$statistic, ref$statistic) && near(t$df, ref$parameter) && near(c(t$ci_low, t$ci_high), ref$conf.int), paste(method, alternative, "matches t.test statistic, df, probability and interval"))
      if (lognormal) check(near(t$geometric_mean_ratio, exp(mean(a) - mean(b))) && near(c(t$original_ci_low, t$original_ci_high), exp(ref$conf.int)), paste(method, alternative, "reports original-scale ratios and bounds"))
    }
  }
  for (method in c("one_sample", "log_one_sample")) for (alternative in c("two.sided", "less", "greater")) {
    lognormal <- startsWith(method, "log_")
    cfg <- list(method = method, y = "response", alternative = alternative, null_value = 9, null_ratio = 9)
    result <- run_analysis(data.frame(response = c(x, NA)), cfg)
    ref <- stats::t.test(if (lognormal) log(x) else x, mu = if (lognormal) log(9) else 9, alternative = alternative)
    t <- result$tables$test
    check(near(t$p_value, ref$p.value) && near(t$estimate, unname(ref$estimate)) && near(c(t$ci_low, t$ci_high), ref$conf.int) && t$n1 == 6, paste(method, alternative, "matches one-sample reference and omits missing values"))
    check(t$null_value == 9 && any(grepl("Hypothesized mean", vapply(result_sheets(result)$tabular_results$rows, `[[`, character(1), "label"))), "One-sample results retain their original-scale hypothesis")
  }
  check(near(run_analysis(d, modifyList(base, list(method = "one_sample", one_sample_group = "B", null_value = 12)))$tables$test$p_value, t.test(y, mu = 12)$p.value), "One-sample group selection does not pool other columns")
  check(fails(run_analysis(d, modifyList(base, list(method = "one_sample"))), "one group"), "One-sample API refuses accidental pooling of several groups")
  partial <- d; partial$response[2] <- NA
  result <- run_analysis(partial, modifyList(base, list(method = "log_paired_t")))
  check(result$tables$test$n_pairs == 5 && result$tables$data_audit$used_rows == 10 && result$tables$data_audit$incomplete_subjects == 1 && near(result$tables$test$p_value, t.test(log(x[-2]), log(y[-2]), paired = TRUE)$p.value), "Lognormal pairing omits whole pairs and audits them")
  check(near(result$normality_samples$value, log(x[-2]) - log(y[-2])) && inherits(ggplot2::ggplot_build(normality_plot(result)), "ggplot_built"), "Lognormal paired QQ diagnostics use the tested log ratios")
  partial$subject[2] <- "S1"
  check(fails(run_analysis(partial, modifyList(base, list(method = "log_paired_t"))), "Duplicate"), "Duplicate pairs are rejected even when a duplicated response is missing")
  for (value in c(0, -1)) {
    bad <- d; bad$response[1] <- value
    check(fails(run_analysis(bad, modifyList(base, list(method = "log_welch"))), "positive"), "Lognormal tests reject zero and negative measurements")
  }
  check(fails(run_analysis(d, modifyList(base, list(method = "log_welch", null_ratio = 0))), "positive"), "Lognormal null ratios must be positive")
  check(fails(run_analysis(d, modifyList(base, list(method = "welch", alternative = "auto"))), "hypothesis"), "Invalid tail settings are rejected")
  check(fails(run_analysis(d, modifyList(base, list(method = "welch", confidence = 1))), "between"), "Invalid t-test confidence levels are rejected")
  shifted <- run_analysis(d, modifyList(base, list(method = "welch", null_value = -2, alternative = "greater")))
  check(near(shifted$tables$test$p_value, t.test(x, y, mu = -2, alternative = "greater")$p.value), "A nonzero hypothesized difference is respected")

  endpoints <- multiple_t_demo()
  cfg <- list(method = "multiple_t", y = "response", group = "group", t_test_mode = "endpoints", test_column = "endpoint", adjustment = "holm")
  result <- run_analysis(endpoints, cfg)
  refs <- lapply(split(endpoints, factor(endpoints$endpoint, levels = unique(endpoints$endpoint))), function(z) t.test(z$response[z$group == "Control"], z$response[z$group == "Treatment"]))
  raw <- vapply(refs, `[[`, numeric(1), "p.value")
  check(identical(result$tables$posthoc$test, unique(endpoints$endpoint)) && near(result$tables$posthoc$p_value, raw) && near(result$tables$posthoc$adjusted_p_value, p.adjust(raw, "holm")), "Endpoint tests use separate variance estimates and one Holm family")
  for (adjustment in c("none", "bonferroni", "BH", "BY")) {
    result <- run_analysis(endpoints, modifyList(cfg, list(adjustment = adjustment, fdr_q = .1)))
    check(near(result$tables$posthoc$adjusted_p_value, p.adjust(raw, adjustment)), paste(adjustment, "family adjustment matches p.adjust"))
    if (adjustment %in% c("BH", "BY")) check(identical(result$tables$posthoc$discovery, unname(p.adjust(raw, adjustment) <= .1)), "FDR discovery decisions use Q")
  }
  for (alternative in c("two.sided", "greater", "less")) {
    result <- run_analysis(endpoints, modifyList(cfg, list(adjustment = "bonferroni", alternative = alternative)))
    z <- endpoints[endpoints$endpoint == "Marker A", ]
    ref <- t.test(z$response[z$group == "Control"], z$response[z$group == "Treatment"], alternative = alternative, conf.level = 1 - .05/3)
    check(near(c(result$tables$posthoc$ci_low[1], result$tables$posthoc$ci_high[1]), ref$conf.int), paste("Bonferroni", alternative, "interval matches independently adjusted confidence level"))
    check(inherits(ggplot2::ggplot_build(make_graph(result, "comparison")), "ggplot_built"), paste(alternative, "comparison graph renders including unbounded limits"))
  }
  result <- run_analysis(endpoints, modifyList(cfg, list(method = "log_multiple_t", adjustment = "sidak", alternative = "less")))
  z <- endpoints[endpoints$endpoint == "Marker A", ]
  ref <- t.test(log(z$response[z$group == "Control"]), log(z$response[z$group == "Treatment"]), conf.level = .95^(1/3), alternative = "less")
  check(near(c(result$tables$posthoc$ratio_ci_low[1], result$tables$posthoc$ratio_ci_high[1]), exp(ref$conf.int)) && result$tables$posthoc$ratio_ci_low[1] == 0, "Sidak one-sided lognormal intervals exponentiate the unbounded limit to zero")
  check(inherits(ggplot2::ggplot_build(make_graph(result, "comparison")), "ggplot_built"), "Lognormal one-sided comparison graph renders with boundary arrows")
  paired_endpoints <- endpoints; paired_endpoints$response[2] <- NA
  result <- run_analysis(paired_endpoints, modifyList(cfg, list(t_test_kind = "paired_t", id = "subject")))
  check(identical(result$tables$posthoc$n_pairs, c(5L, 6L, 6L)) && near(result$tables$posthoc$p_value[1], t.test(endpoints$response[c(1,3:6)], endpoints$response[c(7,9:12)], paired = TRUE)$p.value), "Each endpoint matches subjects locally and keeps complete pairs at other endpoints")
  check(nrow(result$normality_samples) == 17 && length(unique(result$normality_samples$sample)) == 3, "Multiple paired QQ diagnostics keep endpoint differences separate")
  bad <- endpoints; bad$subject[2] <- bad$subject[1]
  check(fails(run_analysis(bad, modifyList(cfg, list(t_test_kind = "paired_t", id = "subject"))), "Duplicate"), "Endpoint paired tests reject duplicate IDs within an endpoint and group")
  bad <- endpoints; bad$group[1] <- "Third"
  check(fails(run_analysis(bad, cfg), "exactly two"), "Endpoint workflow requires a common pair of groups")
  sparse <- endpoints[!endpoints$endpoint %in% "Marker C", ]
  sparse <- rbind(sparse, data.frame(endpoint = "Marker C", subject = c("S1", "S1"), group = c("Control", "Treatment"), response = c(1, 2)))
  result <- run_analysis(sparse, cfg)
  p <- c(raw[1:2], NA_real_)
  check(nrow(result$tables$posthoc) == 3 && all(result$tables$posthoc$family_size == 3) && is.na(result$tables$posthoc$p_value[3]) && is.na(result$tables$posthoc$significant[3]) && near(result$tables$posthoc$adjusted_p_value, p.adjust(p, "holm", n = 3)), "Unavailable endpoints remain in the planned family and never get a false valid decision")

  cols <- rbind(d, data.frame(subject = paste0("S", 1:6), group = "C", response = c(20, 23, 19, 25, 21, 22)))
  cols$response[1] <- NA
  pair_key <- paste(vapply(c("A", "C"), utils::URLencode, character(1), reserved = TRUE), collapse = "::")
  cfg_cols <- modifyList(base, list(method = "multiple_t", t_test_mode = "columns", comparisons = "selected", selected_pairs = pair_key, adjustment = "bonferroni"))
  result <- run_analysis(cols, cfg_cols)
  check(nrow(result$tables$posthoc) == 1 && result$tables$posthoc$group_1 == "A" && result$tables$posthoc$group_2 == "C" && near(result$tables$posthoc$p_value, t.test(x[-1], c(20,23,19,25,21,22))$p.value), "Selected-column workflow tests the intended pair without pooling unrelated groups")
  reverse <- run_analysis(cols, modifyList(cfg_cols, list(reverse_comparisons = TRUE, alternative = "greater")))
  check(reverse$tables$posthoc$group_1 == "C" && near(reverse$tables$posthoc$p_value, t.test(c(20,23,19,25,21,22), x[-1], alternative = "greater")$p.value), "Reversing a one-tailed contrast changes direction explicitly")
  control <- run_analysis(cols, modifyList(cfg_cols, list(comparisons = "control", control = "B", adjustment = "holm")))
  check(nrow(control$tables$posthoc) == 2 && all(control$tables$posthoc$group_1 == "B"), "Multiple t-tests support a selected control group")
  matched <- run_analysis(cols, modifyList(cfg_cols, list(comparisons = "all", t_test_kind = "paired_t")))
  check(identical(matched$tables$posthoc$n_pairs, c(5L, 5L, 6L)), "Column pairs keep their own complete subjects rather than dropping subjects across every condition")
  equal <- run_analysis(cols, modifyList(cfg_cols, list(t_test_kind = "student")))
  check(near(equal$tables$posthoc$p_value, t.test(x[-1], c(20,23,19,25,21,22), var.equal = TRUE)$p.value), "Equal-variance multiple t-tests match pooled two-sample references")
  check(fails(run_analysis(cols, modifyList(cfg_cols, list(adjustment = "tukey"))), "Tukey"), "Tukey cannot be selected as a t-test family correction")
  check(identical(unname(design_methods(family = "one", distribution = "lognormal")), "log_one_sample") && identical(unname(design_methods(family = "two", distribution = "lognormal", matched = TRUE)), "log_paired_t") && identical(unname(design_methods(family = "multiple_t", distribution = "lognormal")), "log_multiple_t"), "Study-design routing exposes the new normal and lognormal t-test families")

  result <- run_analysis(endpoints, modifyList(cfg, list(method = "log_multiple_t", adjustment = "BH", alternative = "less", fdr_q = .1)))
  sheet <- result_sheets(result)$tabular_results
  check(sheet$title == "Multiple t-tests" && any(vapply(sheet$rows, `[[`, character(1), "label") == "Geometric-mean ratio") && any(vapply(sheet$rows, `[[`, character(1), "label") == "FDR discovery"), "Result sheets identify geometric ratios, one-sided limits and FDR discoveries")
  hidden <- result; hidden$config$report_adjusted_p <- FALSE
  check(!any(vapply(result_sheets(hidden)$tabular_results$rows, `[[`, character(1), "label") == "Adjusted P value") && !"adjusted_p_value" %in% names(visible_result_tables(hidden)$posthoc), "Adjusted-p reporting choices remain consistent in sheets and raw exports")
  project <- new_project()
  added <- project_add(project, "tables", list(name = "Endpoints", type = "column", layout = "long", data = endpoints, columns = character(), id_column = "subject")); project <- added$project
  added <- project_add(project, "analyses", list(name = "Log endpoint tests", table_id = added$id, input_data = endpoints, config = result$config, result = result, time = as.character(Sys.time()))); project <- added$project
  reopened <- read_project(project_json(project))
  restored <- reopened$analyses[[added$id]]$result
  check(restored$config$alternative == "less" && restored$config$test_column == "endpoint" && restored$config$scale == "log" && near(restored$tables$posthoc$adjusted_p_value, result$tables$posthoc$adjusted_p_value), "Project JSON preserves t-test tails, endpoint mappings, log scales and correction families")
  out <- file.path("artifacts", "t-tests"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (config in list(modifyList(base, list(method = "log_paired_t", alternative = "less")), list(method = "log_one_sample", y = "response", null_ratio = 9, alternative = "greater"), result$config)) {
    input <- if (config$method == "log_multiple_t") endpoints else if (config$method == "log_one_sample") data.frame(response = x) else d
    expected <- run_analysis(input, config)
    script <- file.path(out, paste0(config$method, ".R")); writeLines(reproduction_script(input, config), script, useBytes = TRUE)
    path <- normalizePath(script); old <- getwd(); setwd(out)
    env <- new.env(parent = globalenv())
    tryCatch(source(path, local = env, encoding = "UTF-8"), finally = setwd(old))
    a <- if (config$method == "log_multiple_t") env$result$tables$posthoc$p_value else env$result$tables$test$p_value
    b <- if (config$method == "log_multiple_t") expected$tables$posthoc$p_value else expected$tables$test$p_value
    check(near(a, b) && env$result$config$alternative == config$alternative, paste("Standalone", config$method, "export runs outside the project and preserves its tail"))
  }
})
