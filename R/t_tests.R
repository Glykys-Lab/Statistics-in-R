# One-sample, two-sample and multiple t-tests on raw or natural-log responses.
# The same functions serve the website, saved projects and standalone R exports.
t_test_methods <- c("welch", "student", "paired_t", "one_sample", "log_welch", "log_student", "log_paired_t", "log_one_sample", "multiple_t", "log_multiple_t")
multiple_t_methods <- c("multiple_t", "log_multiple_t")
one_sample_methods <- c("one_sample", "log_one_sample")

# Map explicit lognormal method names to the underlying t-test design.
t_test_kind <- function(method, config = list()) {
  if (method %in% multiple_t_methods) return(if (is.null(config$t_test_kind)) "welch" else config$t_test_kind)
  sub("^log_", "", method)
}

# Check scalar numerical settings before calling R's statistical functions.
t_test_number <- function(value, label, positive = FALSE) {
  assert(is.numeric(value) && length(value) == 1L && is.finite(value), paste(label, "must be a finite number."))
  if (positive) assert(value > 0, paste(label, "must be positive."))
  value
}

# Keep diagnostic values scoped to the observations and contrast actually tested.
t_test_samples <- function(d, groups, paired = FALSE, label = NULL) {
  if (!nrow(d)) return(NULL)
  if (paired) {
    a <- d[as.character(d$group) == groups[1], , drop = FALSE]
    b <- d[as.character(d$group) == groups[2], , drop = FALSE]
    ids <- intersect(as.character(a$id), as.character(b$id))
    values <- a$y[match(ids, a$id)] - b$y[match(ids, b$id)]
    name <- if (is.null(label)) "Within-pair differences" else paste(label, "Within-pair differences", sep = " / ")
    return(data.frame(value = values, sample = rep(name, length(values))))
  }
  data.frame(value = d$y, sample = if (is.null(label)) as.character(d$group) else paste(label, d$group, sep = " / "))
}

# Calculate one planned test, retaining a readable unavailable row for a family.
# Paired measurements are matched by subject ID separately for each comparison.
t_test_contrast <- function(d, groups, cfg, label, confidence) {
  one <- cfg$method %in% one_sample_methods
  paired <- identical(cfg$t_test_kind, "paired_t") && !one
  a <- d[as.character(d$group) == groups[1], , drop = FALSE]
  b <- if (one) d[FALSE, , drop = FALSE] else d[as.character(d$group) == groups[2], , drop = FALSE]
  if (paired) {
    a <- a[!is.na(a$id) & !is.na(a$y), , drop = FALSE]
    b <- b[!is.na(b$id) & !is.na(b$y), , drop = FALSE]
    ids <- intersect(as.character(a$id), as.character(b$id))
    a <- a[match(ids, a$id), , drop = FALSE]; b <- b[match(ids, b$id), , drop = FALSE]
  } else {
    a <- a[!is.na(a$y), , drop = FALSE]; b <- b[!is.na(b$y), , drop = FALSE]
  }
  x <- a$y; y <- b$y
  null <- if (cfg$scale == "log") log(cfg$null_ratio) else cfg$null_value
  # Undefined or undersized tests remain in the planned correction family.
  test <- tryCatch({
    assert(length(x) >= 2 && (one || length(y) >= 2), if (paired) "At least two complete matched pairs are required." else "At least two observed measurements per tested sample are required.")
    stats::t.test(x, if (one) NULL else y, paired = paired,
      var.equal = identical(cfg$t_test_kind, "student"), mu = null,
      alternative = cfg$alternative, conf.level = confidence)
  }, error = identity)
  valid <- !inherits(test, "error")
  estimate <- if (one) mean(x) else mean(x) - mean(y)
  if (!length(x) || (!one && !length(y))) estimate <- NA_real_
  row <- data.frame(test = label, contrast = if (one) paste(groups[1], "against hypothesis") else paste(groups, collapse = if (cfg$scale == "log") " / " else " minus "),
    group_1 = groups[1], group_2 = if (one) "Hypothesized value" else groups[2],
    method = if (valid) test$method else if (one) "One Sample t-test" else cfg$t_test_kind,
    n1 = length(x), n2 = if (one) NA_integer_ else length(y), n_pairs = if (paired) length(x) else NA_integer_,
    mean_1 = if (length(x)) mean(x) else NA_real_, mean_2 = if (one) null else if (length(y)) mean(y) else NA_real_,
    estimate = estimate, se = if (valid) test$stderr else NA_real_,
    statistic = if (valid) unname(test$statistic) else NA_real_, df = if (valid) unname(test$parameter) else NA_real_,
    p_value = if (valid) test$p.value else NA_real_,
    ci_low = if (valid) test$conf.int[1] else NA_real_, ci_high = if (valid) test$conf.int[2] else NA_real_,
    null_value = if (cfg$scale == "log") cfg$null_ratio else cfg$null_value,
    alternative = cfg$alternative, confidence = cfg$confidence,
    interval = paste0(cfg$confidence * 100, "% ", if (cfg$alternative == "two.sided") "two-sided" else "one-sided", if (cfg$method %in% multiple_t_methods && cfg$adjustment %in% c("bonferroni", "sidak")) paste(" simultaneous", cfg$adjustment) else " unadjusted"),
    status = if (valid) "OK" else paste("Unavailable:", conditionMessage(test)), stringsAsFactors = FALSE)
  if (cfg$scale == "log") {
    row$original_scale_estimate <- exp(row$estimate)
    row$original_ci_low <- exp(row$ci_low); row$original_ci_high <- exp(row$ci_high)
    row$geometric_mean_1 <- exp(row$mean_1); row$geometric_mean_2 <- exp(row$mean_2)
    if (one) row$geometric_mean <- row$original_scale_estimate else row$geometric_mean_ratio <- row$original_scale_estimate
    row$ratio_ci_low <- row$original_ci_low; row$ratio_ci_high <- row$original_ci_high
  }
  list(table = row, rows = if (valid) c(a$source_row, b$source_row) else integer())
}

# Run a t-test workflow and return the same result structure as other analyses.
# Multiple tests can compare column pairs or two groups within each endpoint.
run_t_analysis <- function(data, config) {
  cfg <- modifyList(list(t_test_mode = "columns", t_test_kind = "welch", test_column = NULL,
    one_sample_group = NULL, alternative = "two.sided", null_value = 0, null_ratio = 1,
    confidence = .95, alpha = .05, adjustment = "holm", comparisons = "all", control = NULL,
    selected_pairs = character(), reverse_comparisons = FALSE, fdr_q = .05, scale = "raw"), config)
  assert(is.data.frame(data) && nrow(data) > 0, "Import a nonempty data table first.")
  assert(!anyDuplicated(names(data)), "Column names must be unique.")
  assert(cfg$method %in% t_test_methods, "Unknown t-test method.")
  cfg$t_test_kind <- t_test_kind(cfg$method, cfg)
  multiple <- cfg$method %in% multiple_t_methods; one <- cfg$method %in% one_sample_methods
  lognormal <- startsWith(cfg$method, "log_") || identical(cfg$scale, "log")
  cfg$scale <- if (lognormal) "log" else "raw"
  assert(cfg$t_test_kind %in% c("welch", "student", "paired_t", "one_sample"), "Choose a valid t-test design.")
  assert(!multiple || cfg$t_test_kind != "one_sample", "Multiple t-tests compare two samples; select Welch, equal variance or paired tests.")
  assert(cfg$alternative %in% c("two.sided", "greater", "less") && length(cfg$alternative) == 1L, "Choose a two-tailed, greater-than or less-than hypothesis.")
  for (name in c("confidence", "alpha", "fdr_q")) {
    value <- t_test_number(cfg[[name]], name)
    assert(value > 0 && value < 1, paste(name, "must lie between 0 and 1."))
  }
  if (lognormal) t_test_number(cfg$null_ratio, "Hypothesized geometric mean or ratio", TRUE) else t_test_number(cfg$null_value, "Hypothesized mean or difference")
  assert(cfg$adjustment %in% c("holm", "holm_sidak", "bonferroni", "BH", "BY", "BKY", "none", "sidak"), "Choose a valid t-test multiple-comparison correction (Tukey is not a t-test correction).")
  assert(cfg$t_test_mode %in% c("columns", "endpoints"), "Choose column comparisons or one test per endpoint.")
  paired <- identical(cfg$t_test_kind, "paired_t") && !one
  endpoints <- multiple && cfg$t_test_mode == "endpoints"
  d <- data.frame(y = numeric_column(data, cfg$y),
    group = if (is.null(cfg$group) || identical(cfg$group, "")) factor(rep("All", nrow(data))) else group_column(data, cfg$group), source_row = seq_len(nrow(data)))
  if (paired) {
    assert(length(cfg$id) == 1L && !cfg$id %in% c(cfg$y, cfg$group), "Select a subject ID column distinct from the response and group columns.")
    d$id <- group_column(data, cfg$id)
  }
  if (endpoints) {
    assert(length(cfg$test_column) == 1L && !cfg$test_column %in% c(cfg$y, cfg$group, if (paired) cfg$id), "Select a distinct marker or time-point column for the tests.")
    d$factor2 <- group_column(data, cfg$test_column)
  }
  mapped_missing <- sum(!complete.cases(d[, c("y", "group", if (paired) "id", if (endpoints) "factor2"), drop = FALSE]))
  d <- d[!is.na(d$group) & if (endpoints) !is.na(d$factor2) else TRUE, , drop = FALSE]
  groups <- levels(d$group)
  if (one) {
    if (length(cfg$one_sample_group)) {
      assert(length(cfg$one_sample_group) == 1L && cfg$one_sample_group %in% groups, "Select a valid group for the one-sample test.")
      groups <- cfg$one_sample_group
    }
    assert(length(groups) == 1L, "Select one group to test against the hypothesized value.")
    d <- d[as.character(d$group) == groups, , drop = FALSE]; d$group <- droplevels(d$group)
  } else assert(length(groups) >= 2L && (multiple && !endpoints || length(groups) == 2L), if (endpoints) "Endpoint t-tests require exactly two groups across the table." else "This t-test requires exactly two groups.")
  if (paired) {
    keys <- d[!is.na(d$id), c(if (endpoints) "factor2", "id", "group"), drop = FALSE]
    assert(!anyDuplicated(keys), "Duplicate subject-condition records found within a test. Aggregate technical replicates before analysis.")
  }
  if (lognormal) assert(all(is.na(d$y) | d$y > 0), "Lognormal t-tests require every observed measurement to be positive; zeros and negative values cannot be logged.")
  original <- d
  if (lognormal) d$y <- log(d$y)
  if (multiple && !endpoints) {
    assert(cfg$comparisons %in% c("all", "control", "selected"), "Choose all pairs, a control or selected pairs for multiple t-tests.")
    pairs <- comparison_pairs(groups, cfg)
    keys <- vapply(pairs, function(pair) paste(sort(pair), collapse = "\r"), character(1))
    assert(!anyDuplicated(keys), "Select each pair of groups only once.")
    labels <- vapply(pairs, paste, character(1), collapse = if (lognormal) " / " else " minus ")
    subsets <- rep(list(d), length(pairs))
  } else {
    if (!one && isTRUE(cfg$reverse_comparisons)) groups <- rev(groups)
    if (endpoints) {
      labels <- levels(d$factor2)
      subsets <- lapply(labels, function(label) d[as.character(d$factor2) == label, , drop = FALSE])
      pairs <- rep(list(groups), length(labels))
    } else {
      labels <- if (one) groups else paste(groups, collapse = if (lognormal) " / " else " minus ")
      subsets <- list(d); pairs <- list(groups)
    }
  }
  m <- length(pairs); assert(m > 0L, "Select at least one test.")
  confidence <- if (multiple && cfg$adjustment == "bonferroni") 1 - (1 - cfg$confidence)/m else if (multiple && cfg$adjustment == "sidak") cfg$confidence^(1/m) else cfg$confidence
  # Each comparison has its own complete pairs and error estimate.
  outcomes <- lapply(seq_len(m), function(i) t_test_contrast(subsets[[i]], pairs[[i]], cfg, labels[i], confidence))
  tests <- do.call(rbind, lapply(outcomes, `[[`, "table")); rownames(tests) <- NULL
  if (!multiple) assert(tests$status == "OK", tests$status)
  if (multiple) {
    tests$family_size <- m
    tests$adjusted_p_value <- adjust_family(tests$p_value, cfg$adjustment, n = m)
    fdr <- cfg$adjustment %in% c("BH", "BY", "BKY")
    if (fdr) {
      discoveries <- if (cfg$adjustment == "BKY") bky_rejections(tests$p_value, cfg$fdr_q, n = m) else tests$adjusted_p_value <= cfg$fdr_q
      tests$discovery <- ifelse(is.na(tests$p_value), NA, discoveries); tests$target_FDR <- cfg$fdr_q
    } else {
      tests$significant <- tests$adjusted_p_value < cfg$alpha
      tests$p_summary <- significance_summary(tests$adjusted_p_value)
    }
  } else tests <- annotate_table(tests, cfg$alpha)
  used <- unique(unlist(lapply(outcomes, `[[`, "rows")))
  analyzed <- d[d$source_row %in% used, , drop = FALSE]
  raw <- original[original$source_row %in% used, , drop = FALSE]
  if (paired) {
    analyzed <- analyzed[order(analyzed$id, analyzed$group), , drop = FALSE]
    raw <- raw[order(raw$id, raw$group), , drop = FALSE]
  }
  # Use endpoint-and-group labels for summaries and distribution checks, avoiding
  # a misleading pool of different markers or time points in the same sample.
  display_groups <- if (endpoints) factor(paste(original$factor2, original$group, sep = " / "), levels = unique(paste(original$factor2, original$group, sep = " / "))) else original$group
  log_groups <- if (endpoints) factor(paste(d$factor2, d$group, sep = " / "), levels = levels(display_groups)) else d$group
  missing_rows <- mapped_missing
  incomplete <- if (paired && !multiple) length(unique(original$id[!is.na(original$id)])) - tests$n_pairs else 0L
  tables <- list(descriptive_statistics = summary_table(original$y, display_groups, cfg$confidence),
    summary = summary_table(d$y, log_groups, cfg$confidence),
    data_audit = data.frame(input_rows = nrow(data), missing_rows = missing_rows, incomplete_subjects = incomplete,
      retained_incomplete_subjects = 0L, used_rows = length(used), excluded_rows = nrow(data) - length(used)))
  if (multiple) tables$posthoc <- tests else tables$test <- tests
  if (lognormal) {
    tables$original_summary <- tables$descriptive_statistics
    z <- tables$summary
    tables$geometric_summary <- data.frame(group = z$group, n = z$n, geometric_mean = exp(z$mean), geometric_sd = exp(z$sd), ci_low = exp(z$ci_low), ci_high = exp(z$ci_high))
  }
  # Diagnostics for paired tests concern the within-pair differences actually tested.
  diagnostics <- lapply(seq_len(m), function(i) {
    selected <- d[d$source_row %in% outcomes[[i]]$rows, , drop = FALSE]
    if (!nrow(selected)) return(NULL)
    selected <- droplevels(selected)
    if (paired) selected <- selected[order(selected$id, selected$group), , drop = FALSE]
    z <- normality_table(selected, if (paired) "paired_t" else "summary")
    if (multiple && !is.null(z)) z$sample <- paste(labels[i], z$sample, sep = " / ")
    z
  })
  tables$normality <- do.call(rbind, diagnostics)
  # QQ plots use these same test-specific samples rather than pooling endpoints.
  qq <- lapply(seq_len(m), function(i) t_test_samples(d[d$source_row %in% outcomes[[i]]$rows, , drop = FALSE],
    pairs[[i]], paired, if (multiple) labels[i] else NULL))
  normality_samples <- do.call(rbind, qq)
  observed <- !is.na(original$y)
  tables$distribution_tests <- distribution_tests(original$y[observed], display_groups[observed], cfg$alpha)
  counts <- lapply(split(original$y[observed], display_groups[observed], drop = TRUE), function(x) data.frame(n = length(x), zeros = sum(x == 0), negatives = sum(x < 0)))
  if (length(counts)) {
    counts <- do.call(rbind, counts); tables$distribution_counts <- data.frame(group = rownames(counts), counts, row.names = NULL)
  }
  notes <- c(method_note(cfg$method), paste("Alternative hypothesis:", cfg$alternative, "; comparison direction is printed in each test row."),
    if (one) "The one-sample interval estimates the sample population mean, not its difference from the hypothesized value." else "Contrasts follow group order; greater means the first group exceeds the second by more than the hypothesized difference or ratio.",
    if (lognormal) "Tests use natural logs. Exponentiated estimates and interval limits give geometric means for one-sample tests and geometric-mean ratios for two-sample tests. One-sided limits can be zero or infinity.",
    if (paired) "Subjects are matched by ID within each test. Only complete pairs for that test contribute; input measurements are preserved.",
    if (multiple) paste("One family of", m, "planned tests; correction:", cfg$adjustment, ". Unavailable tests remain in the family. Intervals are simultaneous only for Bonferroni or Sidak; other corrections leave them unadjusted."),
    if (multiple && cfg$adjustment == "sidak") "Sidak simultaneous coverage assumes independent tests; shared groups and repeated endpoints may violate that assumption. Bonferroni does not require independence.",
    if (multiple && cfg$adjustment %in% c("BH", "BY", "BKY")) "FDR discoveries use the selected Q; confidence intervals are unadjusted and are not FDR confidence intervals.",
    if (lognormal) "Data graphs retain original response units and arithmetic mean summaries; geometric-mean estimates and inferential intervals are in the t-test results.",
    "Normality diagnostics are exploratory and do not choose the test or its tail.")
  # If every family member is unavailable, keep observed data for inspection.
  if (!nrow(analyzed)) {
    analyzed <- d[complete.cases(d[, c("y", if (paired) "id"), drop = FALSE]), , drop = FALSE]
    raw <- original[original$source_row %in% analyzed$source_row, , drop = FALSE]
  }
  column_data <- data.frame(y = original$y, group = display_groups)
  list(method = cfg$method, config = cfg, tables = tables, fit = NULL, data = analyzed,
    raw_data = raw, column_data = column_data, normality_samples = normality_samples,
    notes = notes, warnings = unique(tests$status[tests$status != "OK"]))
}

# Supply a reproducible example for tests of two groups at several endpoints.
multiple_t_demo <- function() {
  endpoint <- rep(c("Marker A", "Marker B", "Marker C"), each = 12)
  group <- rep(rep(c("Control", "Treatment"), each = 6), 3)
  data.frame(endpoint = endpoint, subject = rep(rep(paste0("S", 1:6), 2), 3), group = group,
    response = c(8, 10, 9, 12, 11, 13, 11, 12, 14, 13, 15, 14,
      21, 18, 23, 20, 24, 19, 22, 20, 25, 23, 27, 21,
      4, 6, 5, 7, 8, 6, 8, 9, 7, 11, 10, 9), stringsAsFactors = FALSE)
}
