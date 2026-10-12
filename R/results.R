# Read-only result sheets. Calculations remain numeric in result$tables.
# Describe spreadsheet rows, headers, cell types and highlights separately from
# Shiny markup. The same sheet definitions drive browser tables and CSV exports.
# Labeled statistical sheets have a row-label column; plain data frames do not.
result_sheet <- function(title, columns, rows, note = "", headers = columns) list(title = title, columns = columns, headers = headers, rows = rows, note = note)
sheet_row <- function(label = "", values = list(), kind = "data", types = "text", highlight = integer()) {
  list(label = label, values = as.list(values), kind = kind, types = rep(types, length.out = length(values)), highlight = highlight)
}
sheet_number <- function(x, digits = 4L) {
  if (!length(x) || is.na(x)) return("\u2014")
  trimws(format(signif(x, digits), digits = digits, trim = TRUE))
}
sheet_cell <- function(x, type = "text", options = list()) {
  if (is.null(x) || !length(x)) return("")
  if (is.na(x)) return("\u2014")
  if (type == "decision") return(if (isTRUE(x)) "Yes" else "No")
  if (type == "count") return(as.character(x))
  digits <- if (is.null(options$output_digits)) 4L else options$output_digits
  if (type == "number") return(sheet_number(x, digits))
  if (type == "interval") {
    limits <- suppressWarnings(as.numeric(strsplit(x, " to ", fixed = TRUE)[[1]]))
    return(paste(vapply(limits, sheet_number, character(1), digits = digits), collapse = " to "))
  }
  if (type == "p") {
    style <- if (is.null(options$p_style)) "prism" else options$p_style
    pd <- if (is.null(options$p_digits)) 6L else options$p_digits
    if (style == "stars") return(significance_summary(x))
    if (style == "decimal") return(formatC(x, format = "f", digits = pd))
    if (style == "scientific") return(formatC(x, format = "e", digits = pd - 1L))
    return(if (x < .0001) "<0.0001" else trimws(format.pval(x, digits = 4, eps = .0001)))
  }
  as.character(x)
}
sheet_frame <- function(sheet, options = list()) {
  # Apply display formatting only at this boundary; raw numerical exports use
  # result$tables and retain the precision needed for downstream calculations.
  n <- length(sheet$columns)
  rows <- lapply(sheet$rows, function(row) {
    values <- vapply(seq_len(n), function(i) if (i <= length(row$values)) sheet_cell(row$values[[i]], row$types[i], options) else "", character(1))
    if (isFALSE(sheet$statistic_column)) values else c(row$label, values)
  })
  if (!length(rows)) return(data.frame())
  x <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- if (isFALSE(sheet$statistic_column)) sheet$columns else c("Statistic / comparison", sheet$columns); x
}
descriptive_sheet <- function(result) {
  table <- result$tables$descriptive_statistics
  if (is.null(table)) table <- result$tables$original_summary
  if (is.null(table)) table <- result$tables$summary
  rows <- list()
  fields <- list("Number of values" = "n", "Number of missing values" = "missing", "Spacer1" = "blank",
    "Minimum" = "minimum", "25% Percentile" = "q1", "Median" = "median", "75% Percentile" = "q3", "Maximum" = "maximum",
    "Spacer2" = "blank", "Mean" = "mean", "Std. Deviation" = "sd", "Std. Error of Mean" = "sem", "Spacer3" = "blank")
  fields[[paste("Lower", result$config$confidence*100, "% CI")]] <- "ci_low"
  fields[[paste("Upper", result$config$confidence*100, "% CI")]] <- "ci_high"
  for (label in names(fields)) {
    field <- fields[[label]]
    if (field == "blank") rows <- c(rows, list(sheet_row(kind = "blank")))
    else if (field %in% names(table)) rows <- c(rows, list(sheet_row(trimws(label), table[[field]], types = if (field %in% c("n", "missing")) "count" else "number")))
  }
  result_sheet("Descriptive statistics", as.character(table$group), rows,
    "Descriptives use all entered measurements, excluding blank cells. Matched complete-case analyses may use fewer subjects; see Data audit in Additional results.")
}
normality_sheet <- function(result) {
  table <- result$tables$distribution_tests; counts <- result$tables$distribution_counts
  groups <- unique(as.character(table$group)); rows <- list()
  for (distribution in c("Gaussian", "Lognormal")) {
    rows <- c(rows, list(sheet_row(if (distribution == "Gaussian") "Test for normal distribution" else "Test for lognormal distribution", kind = "section")))
    for (test in c("D'Agostino-Pearson", "Shapiro-Wilk")) {
      x <- table[table$distribution == distribution & table$test == test, ]; x <- x[match(groups, x$group), ]
      rows <- c(rows, list(sheet_row(if (test == "D'Agostino-Pearson") "D'Agostino & Pearson test" else "Shapiro-Wilk test", kind = "section"),
        sheet_row(if (test == "Shapiro-Wilk") "W" else "K2", x$statistic, types = "number"),
        sheet_row("P value", x$p_value, types = "p", highlight = seq_along(groups)),
        sheet_row(paste0("Passed ", if (distribution == "Gaussian") "normality" else "lognormality", " test (alpha=", result$config$alpha, ")?"), x$not_rejected, types = "decision"),
        sheet_row("P value summary", significance_summary(x$p_value)), sheet_row(kind = "blank")))
    }
  }
  counts <- counts[match(groups, counts$group), ]
  rows <- c(rows, list(sheet_row("Number of values", counts$n, types = "count"),
    sheet_row("Impossible values in lognormal distributions", kind = "section"), sheet_row("Number of zeros", counts$zeros, types = "count"),
    sheet_row("Number of negative values", counts$negatives, types = "count")))
  unavailable <- unique(table$note[is.na(table$p_value)])
  result_sheet("Tabular results", groups, rows, paste(c("A 'Yes' means the test did not reject at the selected alpha; it does not establish the distribution.", unavailable), collapse = " "))
}
comparison_sheet <- function(result) {
  x <- result$tables$posthoc; cfg <- result$config; fdr <- cfg$adjustment %in% c("BKY", "BH", "BY")
  if (all(c("diff", "lwr", "upr", "p adj") %in% names(x))) {
    x$estimate <- x$diff; x$ci_low <- x$lwr; x$ci_high <- x$upr; x$adjusted_p_value <- x[["p adj"]]
  }
  predicted <- "mean_type" %in% names(x) && any(x$mean_type == "Predicted (LS)")
  estimate_label <- if (predicted) "Predicted (LS) mean diff." else "Mean diff."
  log_scale <- identical(cfg$scale, "log")
  if (log_scale) estimate_label <- paste(estimate_label, "(ln scale)")
  columns <- c(estimate_label, paste0(cfg$confidence*100, "% CI of diff."), if (fdr) "Discovery?" else "Below threshold?", "Summary")
  show_adjusted <- !identical(cfg$report_adjusted_p, FALSE)
  if (show_adjusted) columns <- c(columns, if (fdr) "Adjusted q / P value" else "Adjusted P value")
  columns <- c(columns, "Group 1", "Group 2")
  # The detail section has its own column headings within the same sheet.
  width <- max(length(columns), 8L); columns <- c(columns, rep("", width - length(columns)))
  rows <- list(sheet_row("Number of families", 1L, types = "count"), sheet_row("Number of comparisons per family", nrow(x), types = "count"),
    sheet_row(if (fdr) "Target FDR Q" else "Alpha", if (fdr) cfg$fdr_q else cfg$alpha, types = "number"), sheet_row(kind = "blank"))
  correction <- switch(cfg$adjustment, sidak = "Sidak's", holm_sidak = "Holm-Sidak", bonferroni = "Bonferroni", holm = "Holm", tukey = "Tukey", BKY = "Two-stage BKY", BH = "Benjamini-Hochberg", BY = "Benjamini-Yekutieli", none = "Unadjusted")
  has_intervals <- all(c("estimate", "ci_low", "ci_high") %in% names(x))
  if (!has_intervals) columns[1:2] <- c("Mean rank difference", "Test statistic")
  rows <- c(rows, list(sheet_row(paste(correction, "multiple comparisons"), columns, "header")))
  for (i in seq_len(nrow(x))) {
    probability <- if ("adjusted_p_value" %in% names(x)) x$adjusted_p_value[i] else x$p_value[i]
    decision <- if (fdr) x$discovery[i] else if ("significant" %in% names(x)) x$significant[i] else probability < cfg$alpha
    values <- if (has_intervals) list(x$estimate[i], paste(sheet_number(x$ci_low[i], 12), "to", sheet_number(x$ci_high[i], 12))) else list(x$mean_rank_difference[i], x$statistic[i])
    values <- c(values, list(decision, if (fdr) "" else significance_summary(probability)))
    types <- c("number", if (has_intervals) "interval" else "number", "decision", "text")
    if (show_adjusted) { values <- c(values, list(probability)); types <- c(types, "p") }
    values <- c(values, list(if ("group_1" %in% names(x)) x$group_1[i] else "", if ("group_2" %in% names(x)) x$group_2[i] else "")); types <- c(types, "text", "text")
    label <- if (all(c("group_1", "group_2") %in% names(x))) paste(x$group_1[i], "vs.", x$group_2[i]) else x$contrast[i]
    rows <- c(rows, list(sheet_row(label, values, types = types)))
  }
  if (all(c("mean_1", "mean_2", "n1", "n2") %in% names(x))) {
    headers <- c(if (predicted) "Predicted (LS) Mean 1" else "Mean 1", if (predicted) "Predicted (LS) Mean 2" else "Mean 2", estimate_label, "SE of diff.", "n1", "n2", "t", "DF")
    rows <- c(rows, list(sheet_row(kind = "blank"), sheet_row("Test details", headers, "header")))
    for (i in seq_len(nrow(x))) rows <- c(rows, list(sheet_row(paste(x$group_1[i], "vs.", x$group_2[i]),
      list(x$mean_1[i], x$mean_2[i], x$estimate[i], x$se[i], x$n1[i], x$n2[i], abs(x$statistic[i]), x$df[i]), types = c(rep("number", 4), "count", "count", "number", "number"))))
  }
  if (log_scale && all(c("geometric_mean_ratio", "ratio_ci_low", "ratio_ci_high") %in% names(x))) {
    rows <- c(rows, list(sheet_row(kind = "blank"), sheet_row("Original-scale comparisons", c("Geometric-mean ratio", "Ratio CI"), "header")))
    for (i in seq_len(nrow(x))) rows <- c(rows, list(sheet_row(x$contrast[i], list(x$geometric_mean_ratio[i], paste(sheet_number(x$ratio_ci_low[i], 12), "to", sheet_number(x$ratio_ci_high[i], 12))), types = c("number", "interval"))))
  }
  note <- c(if (predicted) "Predicted means come from the model used for each contrast; n1/n2 count its available observations. Pairwise models may have different predicted means from the overall model." else if (result$method %in% c("rm_anova", "log_rm_anova")) "Matched comparisons use the same complete pairs for Mean 1, Mean 2 and n1/n2." else "Mean 1 and Mean 2 are the observed means of the analyzed groups.",
    if (has_intervals && "interval" %in% names(x)) paste("Intervals:", paste(unique(x$interval), collapse = "; ")))
  if ("status" %in% names(x)) note <- c(note, unique(x$status[x$status != "OK"]))
  result_sheet("Multiple comparisons", columns, rows, paste(note, collapse = " "), headers = rep("", length(columns)))
}
tabular_sheet <- function(result, table_name = "") {
  cfg <- result$config; t <- result$tables; mixed <- result$method %in% c("mixed", "log_mixed")
  rows <- list(sheet_row("Table Analyzed", table_name), sheet_row("Distribution assumption", if (identical(cfg$scale, "log")) "Lognormal (model uses natural logs)" else if (result$method %in% c("kruskal", "friedman")) "Nonparametric" else "Normal (Gaussian)"), sheet_row(kind = "blank"),
    sheet_row(if (mixed) "Mixed-effects model (REML)" else names(analysis_methods)[match(result$method, analysis_methods)], if (!is.null(cfg$id)) paste("Matching by", cfg$id) else "No matching", "section"),
    sheet_row("Alpha", cfg$alpha, types = "number"), sheet_row(kind = "blank"))
  if (result$method %in% c("rm_anova", "log_rm_anova", "mixed", "log_mixed")) rows <- append(rows, list(sheet_row("Assume sphericity?", identical(cfg$sphericity, "none"), types = "decision")), after = 4)
  columns <- c("P value", "P value summary", paste0("Statistically significant (P < ", cfg$alpha, ")?"), "Statistic (DFn, DFd)", "Geisser-Greenhouse epsilon")
  test <- t$test
  if (is.null(test) && !is.null(t$anova)) {
    a <- t$anova; i <- match("group", a$term)
    if (!is.na(i)) test <- data.frame(p_value = a[["Pr(>F)"]][i], statistic = a[["F value"]][i], numerator_df = a$Df[i], denominator_df = tail(a$Df, 1))
  }
  if (!is.null(test)) {
    if ("selected" %in% names(test)) test <- test[test$selected, , drop = FALSE]
    rows <- c(rows, list(sheet_row(if (mixed) "Fixed effect (condition)" else "Omnibus test", columns, "header")))
    for (i in seq_len(nrow(test))) {
      stat <- if (all(c("numerator_df", "denominator_df") %in% names(test))) paste0("F (", sheet_number(test$numerator_df[i]), ", ", sheet_number(test$denominator_df[i]), ") = ", sheet_number(test$statistic[i])) else paste0(sheet_number(test$statistic[i]), if ("df" %in% names(test)) paste0("; df = ", sheet_number(test$df[i])) else "")
      rows <- c(rows, list(sheet_row("Treatment (between columns)", list(test$p_value[i], significance_summary(test$p_value[i]), test$p_value[i] < cfg$alpha, stat, if ("epsilon" %in% names(test)) test$epsilon[i] else NA_real_), types = c("p", "text", "decision", "text", "number"), highlight = 1:3)))
    }
  }
  if (!is.null(t$variance_components)) {
    v <- t$variance_components
    rows <- c(rows, list(sheet_row(kind = "blank"), sheet_row("Random effects", c("SD", "Variance"), "header")))
    for (i in seq_len(nrow(v))) rows <- c(rows, list(sheet_row(if (i == 1) "Individual (between rows)" else "Residual", list(v$sd[i], v$variance[i]), types = "number")))
  }
  if (!is.null(t$matching)) {
    m <- t$matching[1, ]
    rows <- c(rows, list(sheet_row(kind = "blank"), sheet_row("Was the matching effective?", kind = "section"),
      sheet_row(if (mixed) "Chi-square, df" else "F, DFn, DFd", paste(c(sheet_number(m$statistic), if (mixed) sheet_number(m$df) else c(sheet_number(m$numerator_df), sheet_number(m$denominator_df))), collapse = ", ")),
      sheet_row("P value", m$p_value, types = "p"), sheet_row("P value summary", significance_summary(m$p_value)), sheet_row("Is there significant matching?", m$p_value < cfg$alpha, types = "decision")))
    if ("boundary_mixture_p_value" %in% names(m) && is.finite(m$boundary_mixture_p_value)) rows <- c(rows, list(sheet_row("Boundary-mixture P value", m$boundary_mixture_p_value, types = "p")))
  }
  if (!is.null(t$data_summary)) {
    s <- t$data_summary
    rows <- c(rows, list(sheet_row(kind = "blank"), sheet_row("Data summary", kind = "section"), sheet_row("Number of treatments (columns)", s$treatments, types = "count")))
    if (!is.na(s$subjects)) rows <- c(rows, list(sheet_row("Number of subjects (rows)", s$subjects, types = "count")))
    rows <- c(rows, list(sheet_row("Number of missing values", s$missing_values, types = "count"), sheet_row("Number of observations used", t$data_audit$used_rows, types = "count")))
    if (!is.null(t$model) && "subjects" %in% names(t$model)) rows <- c(rows, list(sheet_row("Number of subjects used", t$model$subjects, types = "count")))
  }
  result_sheet("Tabular results", columns, rows, headers = rep("", length(columns)))
}
dataframe_sheet <- function(table, title) {
  sheet <- result_sheet(title, names(table), lapply(seq_len(nrow(table)), function(i) sheet_row(values = as.list(table[i, , drop = FALSE]),
    types = vapply(names(table), function(n) if (is.numeric(table[[n]])) if (n %in% c("p_value", "adjusted_p_value", "Pr(>F)")) "p" else "number" else if (is.logical(table[[n]])) "decision" else "text", character(1)))))
  sheet$statistic_column <- FALSE
  sheet
}
cleaned_worksheet <- function(result, worksheet = NULL) {
  if (is.null(worksheet) || !identical(worksheet$layout, "wide")) return(result$tables$cleaned_data)
  d <- worksheet$data; columns <- worksheet$columns
  assert(nrow(result$tables$cleaned_data) == nrow(d)*length(columns), "Saved worksheet does not match cleaned measurements.")
  for (j in seq_along(columns)) d[[columns[j]]] <- result$tables$cleaned_data[[result$config$y]][((j-1)*nrow(d)+1):(j*nrow(d))]
  d
}
result_sheets <- function(result, table_name = "", worksheet = NULL) {
  # Arrange each analysis into useful result tabs. Worksheet metadata supplies
  # original column names and counts when the fitted data omitted measurements.
  if (!nzchar(table_name) && !is.null(worksheet$name)) table_name <- worksheet$name
  t <- visible_result_tables(result); sheets <- list(); consumed <- character()
  if (isTRUE(result$nested)) {
    sheets <- lapply(names(result$formatted_tables), function(n) dataframe_sheet(result$formatted_tables[[n]], gsub("_", " ", n)))
    names(sheets) <- names(result$formatted_tables)
    attr(sheets, "additional") <- list()
    return(sheets)
  }
  if (result$method == "outliers") {
    x <- t$outlier_summary
    sheets$summary <- result_sheet("Summary", as.character(x$group), list(sheet_row("Method", kind = "section"), sheet_row(paste0("ROUT (Q = ", result$config$rout_q*100, "%)")), sheet_row(kind = "blank"), sheet_row("Number of points", kind = "section"),
      sheet_row("Y values analyzed", x$n, types = "count"), sheet_row("Outliers", x$outliers, types = "count"), sheet_row("Status", x$status)))
    sheets$outliers <- dataframe_sheet(t$outliers, "Outliers")
    sheets$cleaned_data <- dataframe_sheet(cleaned_worksheet(result, worksheet), "Cleaned data")
    consumed <- c("outlier_summary", "outliers", "cleaned_data", "summary", "descriptive_statistics", "data_summary")
  } else {
    if (result$method == "normality") {
      sheets$tabular_results <- normality_sheet(result); consumed <- c("distribution_tests", "distribution_counts")
    } else if (result$method %in% c("anova", "welch_anova", "rm_anova", "mixed", "log_rm_anova", "log_mixed", "kruskal", "friedman")) {
      sheets$tabular_results <- tabular_sheet(result, table_name)
      consumed <- c("test", "anova", "matching", "variance_components", "sphericity", "model", "data_summary")
    } else if (!is.null(t$test) || !is.null(t$coefficients)) {
      field <- if (!is.null(t$test)) "test" else "coefficients"
      sheets$tabular_results <- dataframe_sheet(t[[field]], "Tabular results"); consumed <- field
    }
    if (!is.null(t$posthoc)) { sheets$multiple_comparisons <- comparison_sheet(result); consumed <- c(consumed, "posthoc") }
    if (!is.null(t$descriptive_statistics)) { sheets$descriptive_statistics <- descriptive_sheet(result); consumed <- c(consumed, "descriptive_statistics", "summary", "original_summary") }
    if (!is.null(t$predicted_means)) {
      x <- t$predicted_means
      sheets$predicted_means <- result_sheet("Predicted means", as.character(x$group), list(sheet_row("Number of available values", x$n, types = "count"), sheet_row(if (identical(result$config$scale, "log")) "Predicted (LS) mean (ln response)" else "Predicted (LS) mean", x$mean, types = "number"), sheet_row("Standard error", x$se, types = "number")), "Means from the overall model. Pairwise-model means can differ; see Multiple comparisons.")
      consumed <- c(consumed, "predicted_means")
    }
    if (result$method != "normality" && !is.null(t$distribution_tests)) {
      sheets$normality_distribution <- normality_sheet(result); sheets$normality_distribution$title <- "Normality & lognormality"
      consumed <- c(consumed, "distribution_tests", "distribution_counts")
    }
  }
  attr(sheets, "additional") <- t[setdiff(names(t), consumed)]
  sheets
}
result_export_tables <- function(result, table_name = "", worksheet = NULL, options = list()) {
  sheets <- result_sheets(result, table_name, worksheet)
  arranged <- lapply(sheets, sheet_frame, options = options); names(arranged) <- paste0("sheet_", names(sheets))
  c(arranged, visible_result_tables(result))
}
result_export_choices <- function(result, table_name = "", worksheet = NULL) {
  sheets <- result_sheets(result, table_name, worksheet); raw <- names(visible_result_tables(result))
  c(setNames(paste0("sheet_", names(sheets)), vapply(sheets, `[[`, character(1), "title")), setNames(raw, paste("Raw:", gsub("_", " ", raw))))
}
