# LabStats analysis engine. This file can be sourced without starting the GUI.
# run_analysis(data, config) validates the selected design and returns numeric
# tables, fitted models, settings and audit notes. Formatting is handled elsewhere.
# Subject matching uses IDs; complete-case tests and mixed models intentionally
# handle missing measurements differently. Nested designs dispatch to nested.R.
if (file.exists("R/options.R")) source("R/options.R", local = TRUE)
if (file.exists("R/results.R")) source("R/results.R", local = TRUE)
if (file.exists("R/nested.R")) source("R/nested.R", local = TRUE)
analysis_methods <- c(
  "Descriptive statistics" = "summary", "Unpaired t test (Welch)" = "welch",
  "Unpaired t test (equal variance)" = "student", "Paired t test" = "paired_t",
  "Mann-Whitney test" = "mann_whitney", "Wilcoxon matched-pairs test" = "paired_wilcox",
  "One-way ANOVA" = "anova", "Welch one-way ANOVA" = "welch_anova",
  "Kruskal-Wallis" = "kruskal",
  "Friedman repeated-measures test" = "friedman",
  "Repeated-measures one-way ANOVA" = "rm_anova",
  "Mixed-effects model (REML)" = "mixed",
  "Automatic repeated ANOVA / mixed effects" = "rm_auto",
  "Automatic lognormal repeated ANOVA / mixed effects" = "log_rm_auto",
  "Lognormal repeated-measures ANOVA" = "log_rm_anova",
  "Lognormal mixed-effects model (REML)" = "log_mixed",
  "Nested mixed model: repeated cells" = "nested_repeated",
  "Nested mixed model: between groups" = "nested_between",
  "Normality and lognormality tests" = "normality",
  "Identify outliers (ROUT)" = "outliers",
  "Two-way ANOVA (balanced, independent)" = "two_way",
  "Pearson correlation" = "pearson", "Spearman correlation" = "spearman",
  "Linear regression" = "linear", "Dose-response (4-parameter logistic)" = "dose",
  "Contingency: chi-square" = "chisq", "Contingency: Fisher exact" = "fisher",
  "Kaplan-Meier + log-rank" = "survival"
)

method_note <- function(method) switch(method,
  nested_repeated = "Repeated cells nested within slices and mice. Optional slice x condition effects and unequal cell noise; Satterthwaite inference from the supplied nested-model script.",
  nested_between = "Different groups of cells nested within slices and mice. Separate mouse and slice random effects; optional per-condition slice and residual SDs, with Satterthwaite inference.",
  rm_auto = "Use repeated-measures ANOVA when all subjects have every selected condition; otherwise use REML mixed effects with the available observations.",
  log_rm_auto = "Choose repeated ANOVA or REML mixed effects based on completeness, using natural-log responses and geometric-mean comparisons.",
  summary = "Mean, SD, SEM, median, quartiles, and a t-based confidence interval for each group mean. Each row should represent an independent experimental unit.",
  welch = "Two independent groups. Welch's t test allows unequal variances. The estimate is the first group minus the second in table order.",
  student = "Two independent groups with a scientifically justified equal-variance assumption. Consider Welch's test otherwise.",
  paired_t = "Two conditions on matched subjects. Pairs are matched by subject ID, never by row order. Incomplete pairs are excluded and reported. The t assumption applies to within-pair differences.",
  mann_whitney = "Two independent groups; tests a distributional shift. A median-shift interpretation requires similarly shaped distributions. Uses a normal approximation with continuity correction.",
  paired_wilcox = "Two matched conditions. The signed-rank test assumes symmetric within-pair differences. Uses a normal approximation with continuity correction.",
  anova = "Independent groups, approximately normal residuals, common variance. Follow-ups can compare all pairs, a control, selected pairs, or linear trend. Tukey supplies all-pairs comparisons. Sidak and Bonferroni supply adjusted intervals; other p-value corrections leave intervals unadjusted.",
  welch_anova = "Independent groups, unequal variances allowed. Pairwise Welch tests use the selected p-value adjustment. No equal-variance pooling.",
  kruskal = "Independent groups. Omnibus Kruskal-Wallis with optional Dunn rank comparisons or pairwise Wilcoxon rank-sum tests, using the selected comparison family and correction.",
  friedman = "Matched subjects measured once in every condition. Only complete subjects are used; duplicate subject-condition records are rejected. Optional follow-ups use tie-corrected Dunn rank comparisons or paired Wilcoxon tests.",
  rm_anova = "One within-subject factor, complete subjects only. Geisser-Greenhouse correction is available. Comparisons use paired differences when sphericity is not assumed, or pooled repeated-measures error when it is assumed. All-pairs, control or selected pairs form one correction family; linear trend requires sphericity assumed.",
  mixed = "REML condition effect with compound symmetry and common residual variance; available observations are retained. Choose removal of a nonpositive subject component, a signed covariance component, or a nonnegative random intercept. GG epsilon uses complete subject covariance or an observed-data EM estimate with missing measurements. Missingness must be ignorable. Exact Prism equivalence is unverified.",
  log_rm_anova = "Repeated-measures ANOVA on natural-log responses, with GG correction and comparisons reported as geometric-mean ratios. All observed responses must be positive. Complete subjects only.",
  log_mixed = "REML repeated mixed model on natural-log responses. Includes GG correction, matching diagnostic and geometric-mean ratio comparisons. All observed responses must be positive. Available observations are retained.",
  normality = "Test each column for Gaussian and lognormal distributions using Shapiro-Wilk and D'Agostino-Pearson. Lognormal tests use natural logarithms and are unavailable if any observed value is nonpositive. A non-rejection does not prove the assumed distribution.",
  outliers = "Column ROUT: adaptive Lorentzian constant fit, robust residual scale and the Motulsky-Brown FDR procedure. Q applies separately to each column. Flagged values appear in a separate sheet; the original table is preserved. Cleaned data keep the same rows with flagged responses blank.",
  two_way = "Independent, fully crossed, balanced design only. Includes both main effects and their interaction (ordinary ANOVA). Interpret main effects cautiously when an interaction is present. No post-hoc tests are supplied.",
  pearson = "Linear association between two numeric variables. The test and interval assume independent observations and the usual bivariate-normal model.",
  spearman = "Monotonic association between two numeric variables. Uses an asymptotic p-value; no confidence interval is calculated.",
  linear = "Ordinary least squares with intercept. Inspect residual diagnostics. The shaded band is a confidence interval for the mean response, not a prediction interval.",
  dose = "Unweighted four-parameter logistic curve with free plateaus, fitted to one curve at a time. Supply positive concentrations on their original scale. EC50 is the concentration halfway between fitted plateaus. Intervals use a local Wald approximation; no relative-potency or curve-comparison test is supplied.",
  chisq = "Use one row per category combination, with a nonnegative integer count. Counts are summed within cells. Pearson chi-square without Yates correction; inspect expected counts.",
  fisher = "Use one row per category combination and integer counts. Exact Fisher test; an odds ratio and confidence interval are returned only for 2-by-2 tables.",
  survival = "One independent subject per row. Time must be nonnegative; event must be numeric 1 = event, 0 = right-censored. Includes Kaplan-Meier estimates and an unadjusted log-rank comparison. No covariate adjustment.")

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
numeric_column <- function(data, name) {
  # Reject text and non-finite measurements rather than silently coercing them.
  assert(length(name) == 1L && name %in% names(data), "Select a valid numeric column.")
  assert(is.numeric(data[[name]]), paste0("Column '", name, "' must be numeric. Check for text, units, or decimal separators."))
  assert(!any(is.infinite(data[[name]])), paste0("Column '", name, "' contains infinite values."))
  data[[name]]
}
group_column <- function(data, name) {
  # Preserve an explicit factor order or first appearance, never alphabetical order.
  assert(length(name) == 1L && name %in% names(data), "Select a valid group/category column.")
  x <- trimws(as.character(data[[name]])); x[x == ""] <- NA_character_
  order <- if (is.factor(data[[name]])) trimws(levels(data[[name]])) else unique(x)
  factor(x, levels = unique(order[!is.na(order) & nzchar(order)]))
}
summary_table <- function(y, group, confidence = .95) {
  do.call(rbind, lapply(split(y, group, drop = TRUE), function(x) {
    missing <- sum(is.na(x)); x <- x[!is.na(x)]
    n <- length(x); se <- if (n > 1) sd(x)/sqrt(n) else NA_real_
    margin <- if (n > 1) qt((1 + confidence)/2, n - 1) * se else NA_real_
    data.frame(n = n, missing = missing, minimum = if (n) min(x) else NA_real_, maximum = if (n) max(x) else NA_real_, mean = if (n) mean(x) else NA_real_, sd = sd(x), sem = se, median = median(x),
               q1 = if (n) unname(quantile(x, .25)) else NA_real_, q3 = if (n) unname(quantile(x, .75)) else NA_real_,
               ci_low = if (n) mean(x) - margin else NA_real_, ci_high = if (n) mean(x) + margin else NA_real_)
  })) |> (function(x) data.frame(group = rownames(x), x, row.names = NULL))()
}
htest_table <- function(test) {
  data.frame(method = test$method, statistic = if (length(test$statistic)) unname(test$statistic[1]) else NA_real_,
    df = if (length(test$parameter)) unname(test$parameter[1]) else NA_real_, p_value = test$p.value,
    estimate = if (length(test$estimate) == 1) unname(test$estimate) else NA_real_,
    ci_low = if (length(test$conf.int)) test$conf.int[1] else NA_real_,
    ci_high = if (length(test$conf.int)) test$conf.int[2] else NA_real_)
}
pairwise_table <- function(object) {
  x <- as.data.frame(as.table(object$p.value), stringsAsFactors = FALSE)
  names(x) <- c("group_1", "group_2", "adjusted_p_value")
  x[!is.na(x$adjusted_p_value), , drop = FALSE]
}

# Design choices narrow the available analyses without a preliminary p-value gate.
design_methods <- function(table_type = "column", family = "many", matched = FALSE,
                           distribution = "normal", variance = "unequal", repeated = "complete") {
  if (table_type == "xy") return(analysis_methods[analysis_methods %in% c("linear", "pearson", "spearman")])
  if (table_type == "nested") return(analysis_methods[analysis_methods %in% nested_methods])
  if (table_type == "other") return(analysis_methods[!analysis_methods %in% nested_methods])
  if (family %in% c("normality", "outliers")) return(analysis_methods[analysis_methods == family])
  if (distribution == "lognormal" && matched && family == "many") return(analysis_methods[analysis_methods == switch(repeated, mixed = "log_mixed", auto = "log_rm_auto", "log_rm_anova")])
  method <- if (family == "summary") "summary" else if (family == "two") {
    if (matched) { if (distribution == "normal") "paired_t" else "paired_wilcox" }
    else if (distribution == "rank") "mann_whitney" else if (variance == "equal") "student" else "welch"
  } else if (matched) {
    if (distribution == "rank") "friedman" else switch(repeated, mixed = "mixed", auto = "rm_auto", "rm_anova")
  } else if (distribution == "rank") "kruskal" else if (variance == "equal") "anova" else "welch_anova"
  analysis_methods[analysis_methods == method]
}

normality_table <- function(d, method, fit = NULL) {
  samples <- if (method %in% c("paired_t", "paired_wilcox")) {
    v <- split(d$y, d$group); setNames(list(v[[1]] - v[[2]]), "Within-pair differences")
  } else if (method %in% c("linear", "dose")) list("Model residuals" = residuals(fit))
  else if (method %in% c("pearson", "spearman")) list(X = d$x, Y = d$y)
  else if (method %in% c("rm_anova", "mixed", "friedman")) {
    # Within-subject contrasts, rather than pooled correlated residuals.
    subjects <- split(d, d$id)
    pairs <- combn(levels(d$group), 2, simplify = FALSE)
    setNames(lapply(pairs, function(p) unlist(lapply(subjects, function(s) {
      a <- s$y[s$group == p[1]]; b <- s$y[s$group == p[2]]
      if (length(a) == 1 && length(b) == 1) a - b else NULL
    }))), vapply(pairs, paste, character(1), collapse = " minus "))
  } else if ("group" %in% names(d)) split(d$y, d$group) else return(NULL)
  do.call(rbind, lapply(names(samples), function(name) {
    x <- samples[[name]]; n <- length(x)
    status <- if (n < 3) "Requires at least 3 values" else if (n > 5000) "Maximum 5000 values; no subsampling" else if (length(unique(x)) < 2) "Constant values" else "Exploratory; assess QQ plots and study design"
    test <- if (n >= 3 && n <= 5000 && length(unique(x)) > 1) shapiro.test(x) else NULL
    data.frame(sample = name, n = n, W = if (is.null(test)) NA_real_ else unname(test$statistic),
      p_value = if (is.null(test)) NA_real_ else test$p.value, note = status)
  }))
}

repeated_anova <- function(d, correction = "gg") {
  # Work on matched subject rows so the covariance used for GG correction has
  # the same condition ordering as the repeated-measures model.
  fit <- lm(y ~ id + group, data = d)
  a <- anova(fit); k <- nlevels(d$group)
  wide <- matrix(d$y, ncol = k, byrow = TRUE)
  basis <- qr.Q(qr(contr.helmert(k)))
  s <- cov(wide %*% basis)
  assert(sum(s^2) > 0 && a["Residuals", "Mean Sq"] > 0, "Repeated-measures error variance is zero; the F test is undefined.")
  epsilon <- max(1/(k - 1), min(1, sum(diag(s))^2 / ((k - 1) * sum(s^2))))
  f <- a["group", "F value"]; df1 <- k - 1; df2 <- a["Residuals", "Df"]
  tests <- data.frame(correction = c("Sphericity assumed", "Geisser-Greenhouse"),
    statistic = f, numerator_df = c(df1, epsilon * df1), denominator_df = c(df2, epsilon * df2),
    epsilon = c(1, epsilon), p_value = c(pf(f, df1, df2, lower.tail = FALSE), pf(f, epsilon * df1, epsilon * df2, lower.tail = FALSE)),
    selected = c(correction == "none", correction == "gg"))
  list(fit = fit, tests = tests)
}

# D'Agostino's skew transformation and Anscombe-Glynn kurtosis transformation.
# References: D'Agostino et al. (1990), Anscombe & Glynn (1983).
dagostino_pearson <- function(x) {
  n <- length(x)
  if (n < 8 || length(unique(x)) < 2) return(c(K2 = NA_real_, p_value = NA_real_, skew_z = NA_real_, kurtosis_z = NA_real_))
  # Standardization avoids overflow and preserves location/scale invariance.
  z <- (x - mean(x)) / sd(x)
  m2 <- mean(z^2); skew <- mean(z^3)/m2^1.5; kurt <- mean(z^4)/m2^2
  y <- skew * sqrt((n + 1) * (n + 3)/(6 * (n - 2)))
  beta <- 3 * (n^2 + 27*n - 70) * (n + 1) * (n + 3)/((n - 2) * (n + 5) * (n + 7) * (n + 9))
  w2 <- sqrt(2 * (beta - 1)) - 1
  skew_z <- asinh(y/sqrt(2/(w2 - 1))) / sqrt(.5 * log(w2))
  expectation <- 3 * (n - 1)/(n + 1)
  variance <- 24*n*(n - 2)*(n - 3)/((n + 1)^2*(n + 3)*(n + 5))
  standardized <- (kurt - expectation)/sqrt(variance)
  root_beta <- 6*(n^2 - 5*n + 2)/((n + 7)*(n + 9)) * sqrt(6*(n + 3)*(n + 5)/(n*(n - 2)*(n - 3)))
  a <- 6 + 8/root_beta * (2/root_beta + sqrt(1 + 4/root_beta^2))
  denom <- 1 + standardized * sqrt(2/(a - 4))
  kurt_z <- if (denom == 0) NA_real_ else ((1 - 2/(9*a)) - sign(denom)*((1 - 2/a)/abs(denom))^(1/3))/sqrt(2/(9*a))
  k2 <- skew_z^2 + kurt_z^2
  c(K2 = k2, p_value = pchisq(k2, 2, lower.tail = FALSE), skew_z = skew_z, kurtosis_z = kurt_z)
}

distribution_tests <- function(y, group, alpha = .05) {
  do.call(rbind, lapply(names(split(y, group, drop = TRUE)), function(name) {
    x <- y[as.character(group) == name]; x <- x[!is.na(x)]
    do.call(rbind, lapply(c("Gaussian", "Lognormal"), function(distribution) {
      valid <- distribution == "Gaussian" || all(x > 0)
      values <- if (valid && distribution == "Lognormal") log(x) else x
      n <- length(x); constant <- length(unique(values)) < 2
      shapiro <- if (valid && !constant && n >= 3 && n <= 5000) shapiro.test(values) else NULL
      dp <- if (valid && !constant && n >= 8) dagostino_pearson(values) else c(K2 = NA_real_, p_value = NA_real_)
      p <- c(if (is.null(shapiro)) NA_real_ else shapiro$p.value, dp["p_value"])
      status <- if (!valid) rep("Unavailable: nonpositive observations", 2) else if (constant) rep("Unavailable: constant or empty sample", 2) else c(
        if (n < 3 || n > 5000) "Shapiro-Wilk requires 3 to 5000 observations" else "Non-rejection does not establish normality",
        if (n < 8) "D'Agostino-Pearson requires at least 8 observations" else if (n < 20) "Small-sample approximation (n < 20)" else "Non-rejection does not establish normality")
      data.frame(group = name, distribution = distribution, test = c("Shapiro-Wilk", "D'Agostino-Pearson"), n = n,
        statistic = c(if (is.null(shapiro)) NA_real_ else unname(shapiro$statistic), dp["K2"]), p_value = unname(p),
        alpha = alpha, not_rejected = ifelse(is.na(p), NA, p >= alpha), note = status, row.names = NULL)
    }))
  }))
}

adjust_p <- function(p, adjustment) {
  if (adjustment == "sidak") -expm1(length(p) * log1p(-p)) else p.adjust(p, method = adjustment)
}

significance_summary <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p <= .0001, "****", ifelse(p <= .001, "***", ifelse(p <= .01, "**", ifelse(p <= .05, "*", "ns")))))
}

annotate_table <- function(table, alpha = .05) {
  columns <- intersect(c("adjusted_p_value", "p adj", "p_value", "Pr(>F)", "Pr(>|t|)"), names(table))
  if (!length(columns)) return(table)
  p <- table[[columns[1]]]
  table$p_summary <- significance_summary(p); table$significant <- ifelse(is.na(p), NA, p < alpha)
  table
}

mean_comparisons <- function(d, cfg, kind = "paired", fit = NULL) {
  groups <- levels(d$group); pairs <- comparison_pairs(groups, cfg)
  if (!length(pairs)) return(NULL)
  m <- length(pairs)
  tail <- (1 - cfg$confidence)/2
  adjusted_tail <- if (cfg$adjustment == "sidak") -expm1(log(cfg$confidence)/m)/2 else if (cfg$adjustment == "bonferroni") tail/m else tail
  output <- lapply(pairs, function(p) {
    a <- d[d$group == p[1], , drop = FALSE]; b <- d[d$group == p[2], , drop = FALSE]
    status <- "OK"; estimate <- se <- df <- mean_1 <- mean_2 <- NA_real_; n_pairs <- NA_integer_
    n1 <- nrow(a); n2 <- nrow(b); mean_type <- "Observed"
    if (kind == "paired") {
      ids <- intersect(as.character(a$id), as.character(b$id))
      values_1 <- a$y[match(ids, a$id)]; values_2 <- b$y[match(ids, b$id)]
      delta <- values_1 - values_2
      mean_1 <- mean(values_1); mean_2 <- mean(values_2); n1 <- n2 <- length(ids)
      n_pairs <- length(delta); estimate <- mean(delta); se <- sd(delta)/sqrt(n_pairs); df <- n_pairs - 1
    } else if (kind %in% c("welch", "ordinary", "rm_pooled")) {
      mean_1 <- mean(a$y); mean_2 <- mean(b$y)
      estimate <- mean(a$y) - mean(b$y)
      if (kind == "welch") {
        v1 <- var(a$y)/nrow(a); v2 <- var(b$y)/nrow(b)
        se <- sqrt(v1 + v2); df <- (v1 + v2)^2/(v1^2/(nrow(a)-1) + v2^2/(nrow(b)-1))
      } else {
        se <- sqrt(sigma(fit)^2 * (1/nrow(a) + 1/nrow(b))); df <- df.residual(fit)
        if (kind == "rm_pooled") n_pairs <- nrow(a)
      }
    } else {
      mean_type <- "Predicted (LS)"
      pair_d <- droplevels(d[d$group %in% p, , drop = FALSE])
      pair_fit <- fit
      if (kind == "mixed_pairs") {
        coverage <- table(pair_d$id, pair_d$group); n_pairs <- sum(rowSums(coverage) == 2)
        if (n_pairs < 2 || nrow(pair_d) <= nlevels(pair_d$id) + 1) status <- "Unavailable: insufficient overlap for this pair"
        else pair_fit <- tryCatch(fit_mixed_options(pair_d, cfg), error = function(e) NULL)
        if (is.null(pair_fit)) status <- "Unavailable: pairwise model did not converge"
      }
      if (status == "OK") {
        lev <- levels(pair_fit$data$group)
        design <- model.matrix(~group, data.frame(group = factor(p, levels = lev)))
        contrast <- design[1, ] - design[2, ]
        predicted <- as.numeric(design %*% mixed_beta(pair_fit)); mean_1 <- predicted[1]; mean_2 <- predicted[2]
        estimate <- as.numeric(contrast %*% mixed_beta(pair_fit)); se <- sqrt(as.numeric(contrast %*% vcov(pair_fit) %*% contrast))
        df <- min(mixed_coefficients(pair_fit)[contrast != 0, "DF"])
      }
    }
    if (!is.finite(se) || se <= 0 || !is.finite(df) || df <= 0) status <- if (status == "OK") "Unavailable: no estimable comparison error" else status
    test <- if (status == "OK") estimate/se else NA_real_
    probability <- if (status == "OK") 2*pt(-abs(test), df) else NA_real_
    margin <- if (status == "OK") qt(1 - adjusted_tail, df)*se else NA_real_
    data.frame(group_1 = p[1], group_2 = p[2], contrast = paste(p, collapse = " minus "), estimate = estimate,
      mean_1 = mean_1, mean_2 = mean_2, mean_type = mean_type, n1 = n1, n2 = n2,
      se = se, df = df, n_pairs = n_pairs, statistic = test, p_value = probability, ci_low = estimate - margin, ci_high = estimate + margin,
      interval = if (cfg$adjustment %in% c("sidak", "bonferroni")) paste(cfg$confidence*100, "% simultaneous", cfg$adjustment) else paste(cfg$confidence*100, "% unadjusted"), status = status)
  })
  result <- do.call(rbind, output)
  # Preserve family size including unavailable contrasts rather than shrinking it.
  result$adjusted_p_value <- adjust_family(result$p_value, cfg$adjustment, n = m)
  if (cfg$scale == "log") {
    result$geometric_mean_ratio <- exp(result$estimate); result$ratio_ci_low <- exp(result$ci_low); result$ratio_ci_high <- exp(result$ci_high)
  }
  result
}

subject_matrix <- function(d) {
  subjects <- levels(d$id); groups <- levels(d$group)
  x <- matrix(NA_real_, length(subjects), length(groups), dimnames = list(subjects, groups))
  x[cbind(match(d$id, subjects), match(d$group, groups))] <- d$y
  x
}

gg_epsilon <- function(covariance) {
  k <- nrow(covariance)
  if (k == 2) return(1)
  center <- diag(k) - matrix(1/k, k, k); s <- center %*% covariance %*% center
  assert(all(is.finite(s)) && sum(s^2) > 0, "GG epsilon cannot be estimated from zero or undefined covariance.")
  max(1/(k - 1), min(1, sum(diag(s))^2/((k - 1)*sum(s^2))))
}

# Observed-data multivariate Gaussian EM. Missing outcomes are integrated out;
# conditional covariance contributes to E[YY']; values are never imputed for tests.
em_covariance <- function(x, maxiter = 2000L, tolerance = 1e-7) {
  x <- x[rowSums(!is.na(x)) > 0, , drop = FALSE]; n <- nrow(x); k <- ncol(x)
  assert(n >= 3 && all(colSums(!is.na(x)) >= 2), "GG covariance requires at least three subjects and two observations per condition.")
  overlap <- crossprod(!is.na(x))
  assert(all(overlap >= 2), "GG correction with missing data requires at least two shared subjects for every condition pair. Choose sphericity assumed or collect overlapping observations.")
  mu <- colMeans(x, na.rm = TRUE); s <- diag(apply(x, 2, var, na.rm = TRUE), k)
  assert(all(diag(s) > 0), "GG covariance requires variation in every condition.")
  patterns <- split(seq_len(n), apply(is.na(x), 1, paste, collapse = ""))
  converged <- FALSE
  for (iteration in seq_len(maxiter)) {
    expected <- x; uncertainty <- matrix(0, k, k)
    for (rows in patterns) {
      missing <- which(is.na(x[rows[1], ])); observed <- setdiff(seq_len(k), missing)
      if (!length(missing)) next
      gain <- s[missing, observed, drop = FALSE] %*% solve(s[observed, observed, drop = FALSE])
      expected[rows, missing] <- sweep(sweep(x[rows, observed, drop = FALSE], 2, mu[observed]) %*% t(gain), 2, mu[missing], "+")
      uncertainty[missing, missing] <- uncertainty[missing, missing, drop = FALSE] + length(rows) * (s[missing, missing, drop = FALSE] - gain %*% s[observed, missing, drop = FALSE])
    }
    next_mu <- colMeans(expected); centered <- sweep(expected, 2, next_mu)
    next_s <- (crossprod(centered) + uncertainty)/n
    change <- max(abs(next_s - s))/(max(abs(s)) + .Machine$double.eps) + max(abs(next_mu - mu))/(max(abs(mu)) + 1)
    s <- (next_s + t(next_s))/2; mu <- next_mu
    assert(min(eigen(s, symmetric = TRUE, only.values = TRUE)$values) > .Machine$double.eps * max(diag(s)), "Unstructured covariance is singular; GG correction cannot be estimated reliably for these missing-data patterns.")
    if (change < tolerance) { converged <- TRUE; break }
  }
  assert(converged, "Missing-data GG covariance estimation did not converge; choose sphericity assumed or review data coverage.")
  list(covariance = s*n/(n - 1), iterations = iteration)
}

mixed_epsilon <- function(d) {
  x <- subject_matrix(d)
  if (ncol(x) == 2) return(list(epsilon = 1, source = "Two conditions: epsilon = 1", iterations = 0))
  if (!anyNA(x)) return(list(epsilon = gg_epsilon(cov(x)), source = "Observed complete-subject covariance", iterations = 0))
  em <- em_covariance(x)
  list(epsilon = gg_epsilon(em$covariance), source = "Observed-data Gaussian EM unstructured covariance (not Prism-validated)", iterations = em$iterations)
}

rout_column <- function(x, q = .01, tolerance = 1e-9, maxiter = 1000L) {
  # Identify unusually large deviations using a robust center/scale and FDR rule;
  # the caller decides whether to flag observations or create cleaned data.
  n <- length(x); assert(all(is.finite(x)), "ROUT requires finite observed values.")
  assert(length(q) == 1 && is.finite(q) && q > 0 && q < 1, "ROUT Q must lie between 0 and 1.")
  location <- median(x); scale <- NA_real_; status <- "OK"; converged <- FALSE
  if (n < 3) status <- "Unavailable: at least three observations required"
  else if (length(unique(x)) < 2) status <- "Unavailable: constant sample"
  else {
    # Constant-model specialization of the adaptive Lorentzian derivative/Hessian
    # step in Motulsky & Brown 2006; recompute RSDR at each iteration (K = 1).
    for (iteration in seq_len(maxiter)) {
      residual <- x - location
      scale <- unname(quantile(abs(residual), .6827, type = 7))*n/(n - 1)
      if (scale <= .Machine$double.eps*max(1, max(abs(x)))) { status <- "Unavailable: zero robust residual scale"; break }
      weight <- 1/(1 + (residual/scale)^2)
      next_location <- sum(weight*x)/sum(weight)
      if (abs(next_location - location) <= tolerance*max(scale, abs(location), 1)) { location <- next_location; converged <- TRUE; break }
      location <- next_location
    }
    if (!converged && status == "OK") status <- "Unavailable: robust fit did not converge"
  }
  residual <- x - location
  if (status == "OK") scale <- unname(quantile(abs(residual), .6827, type = 7))*n/(n - 1)
  p <- if (status == "OK") 2*pt(-abs(residual)/scale, df = n - 1) else rep(NA_real_, n)
  order <- order(abs(residual)); threshold <- rep(NA_real_, n); flag <- rep(FALSE, n)
  # Published scan starts at integer(0.70*N), stopping at the first rejection.
  candidates <- if (n >= 3) seq.int(max(1L, as.integer(.70*n)), n) else integer()
  if (status == "OK") for (i in candidates) {
    threshold[order[i]] <- q*(n - i + 1)/n
    if (p[order[i]] < threshold[order[i]]) { flag[order[i:n]] <- TRUE; break }
  }
  list(location = location, scale = scale, residual = residual, p_value = p, threshold = threshold, outlier = flag, status = status)
}

run_analysis <- function(data, config) {
  # Route hierarchical data before ordinary-column validation: nested models need
  # Mouse/Slice/Cell columns rather than a single group and subject mapping.
  if (config$method %in% nested_methods) return(run_nested_analysis(data, config))
  defaults <- list(method = "summary", y = NULL, group = NULL, x = NULL, id = NULL,
    factor2 = NULL, event = NULL, count = NULL, confidence = .95, adjustment = "holm", direction = "decreasing", sphericity = "gg",
    alpha = .05, comparisons = "all", control = NULL, scale = "raw", rout_q = .01,
    selected_pairs = character(), reverse_comparisons = FALSE, fdr_q = .05, rank_test = "wilcox",
    random_effect = "keep", glm_initial = TRUE, report_adjusted_p = TRUE, report_descriptive = TRUE,
    report_goodness = FALSE, report_effects = FALSE, report_aicc = FALSE,
    residual_normality = FALSE, residual_variance = FALSE, diagnostic_graphs = character())
  if (identical(config$method, "anova")) defaults$adjustment <- "tukey"
  cfg <- modifyList(defaults, config)
  cfg$selected_pairs <- as.character(unlist(cfg$selected_pairs))
  cfg$diagnostic_graphs <- as.character(unlist(cfg$diagnostic_graphs))
  if (cfg$method %in% c("rm_auto", "log_rm_auto")) {
    group <- group_column(data, cfg$group); id <- group_column(data, cfg$id); y <- numeric_column(data, cfg$y)
    coverage <- table(id, group)
    incomplete <- anyNA(y) || anyNA(id) || anyNA(group) || any(coverage != 1)
    requested <- cfg$method
    cfg$method <- if (requested == "log_rm_auto") if (incomplete) "log_mixed" else "log_rm_anova" else if (incomplete) "mixed" else "rm_anova"
    result <- run_analysis(data, cfg)
    result$config$requested_method <- requested
    result$notes <- c(paste("Automatic selection:", names(analysis_methods)[match(result$method, analysis_methods)],
      if (incomplete) "because selected measurements are incomplete." else "because every subject has every selected condition."), result$notes)
    return(result)
  }
  if (cfg$method %in% c("log_rm_anova", "log_mixed")) {
    original <- cfg$method
    cfg$method <- if (original == "log_rm_anova") "rm_anova" else "mixed"; cfg$scale <- "log"
    result <- run_analysis(data, cfg); result$method <- original; result$config$method <- original
    result$notes <- c(method_note(original), result$notes)
    return(result)
  }
  assert(is.data.frame(data) && nrow(data) > 0, "Import a nonempty data table first.")
  assert(!anyDuplicated(names(data)), "Column names must be unique.")
  assert(cfg$method %in% analysis_methods, "Unknown analysis method.")
  assert(is.numeric(cfg$confidence) && length(cfg$confidence) == 1 && is.finite(cfg$confidence) && cfg$confidence > 0 && cfg$confidence < 1, "Confidence level must lie between 0 and 1.")
  assert(cfg$adjustment %in% c("holm", "holm_sidak", "bonferroni", "BH", "BY", "BKY", "none", "sidak", "tukey"), "Unknown multiple-comparison adjustment.")
  assert(cfg$adjustment != "tukey" || cfg$method == "anova", "Tukey comparisons require ordinary one-way ANOVA.")
  assert(cfg$sphericity %in% c("gg", "none"), "Unknown sphericity setting.")
  assert(length(cfg$alpha) == 1 && is.finite(cfg$alpha) && cfg$alpha > 0 && cfg$alpha < 1, "Alpha must lie between 0 and 1.")
  assert(cfg$comparisons %in% c("all", "control", "none", "selected", "trend"), "Choose a valid comparison family.")
  assert(cfg$random_effect %in% c("keep", "remove", "full"), "Choose a valid random-effect policy.")
  assert(cfg$rank_test %in% c("dunn", "wilcox"), "Choose Dunn or pairwise Wilcoxon rank comparisons.")
  assert(length(cfg$fdr_q) == 1 && is.finite(cfg$fdr_q) && cfg$fdr_q > 0 && cfg$fdr_q < 1, "Target FDR must lie between 0 and 1.")
  assert(cfg$comparisons != "trend" || cfg$method %in% c("anova", "rm_anova", "mixed"), "Linear trend requires parametric one-way analysis.")
  assert(cfg$comparisons != "trend" || cfg$method == "anova" || cfg$sphericity == "none", "Linear trend requires sphericity assumed for repeated data.")
  assert(all(cfg$diagnostic_graphs %in% c("residual", "homoscedasticity", "residual_qq", "residual_heatmap", "comparison", "ranks", "differences")), "Choose valid diagnostic graphs.")
  assert(cfg$scale %in% c("raw", "log"), "Unknown response scale.")
  method <- cfg$method
  notes <- method_note(method); warnings <- character(); tables <- list(); fit <- NULL
  calc <- function() {
    if (method %in% c("chisq", "fisher")) {
      assert(cfg$group != cfg$factor2, "Select two different category columns.")
      d <- data.frame(group = group_column(data, cfg$group), category = group_column(data, cfg$factor2), count = numeric_column(data, cfg$count))
    } else if (method == "survival") {
      d <- data.frame(y = numeric_column(data, cfg$y), event = numeric_column(data, cfg$event), group = group_column(data, cfg$group))
    } else if (method %in% c("pearson", "spearman", "linear", "dose")) {
      assert(cfg$x != cfg$y, "Select different X and Y columns.")
      d <- data.frame(x = numeric_column(data, cfg$x), y = numeric_column(data, cfg$y))
    } else {
      d <- data.frame(y = numeric_column(data, cfg$y), group = if (is.null(cfg$group) || identical(cfg$group, "")) factor(rep("All", nrow(data))) else group_column(data, cfg$group))
      if (method == "outliers") d$source_row <- seq_len(nrow(data))
      if (method %in% c("paired_t", "paired_wilcox", "friedman", "rm_anova", "mixed")) {
        assert(cfg$id != cfg$group && cfg$id != cfg$y, "Select a subject ID column distinct from the response and condition columns.")
        d$id <- group_column(data, cfg$id)
      }
      if (method == "two_way") {
        assert(cfg$group != cfg$factor2, "Select two different factors.")
        d$factor2 <- group_column(data, cfg$factor2)
      }
    }
    if ("count" %in% names(d)) assert(all(is.na(d$count) | (d$count >= 0 & d$count == floor(d$count))), "Counts must be nonnegative integers.")
    if (method == "survival") {
      assert(all(is.na(d$event) | d$event %in% c(0, 1)), "Event must be 0 (censored) or 1 (event).")
      assert(all(is.na(d$y) | d$y >= 0), "Survival times must be nonnegative.")
    }
    # Check duplicate subject-condition records even when one outcome is missing.
    if ("id" %in% names(d)) {
      keys <- d[complete.cases(d[, c("id", "group")]), c("id", "group")]
      assert(!anyDuplicated(keys), "Duplicate subject-condition records found. Aggregate technical replicates before analysis.")
    }
    original_groups <- if ("group" %in% names(d)) levels(d$group) else NULL
    if ("y" %in% names(d) && "group" %in% names(d) && method != "survival") {
      tables$descriptive_statistics <<- summary_table(d$y, d$group, cfg$confidence)
      subjects <- if ("id" %in% names(d)) unique(as.character(d$id[!is.na(d$id)])) else character()
      if (length(subjects) && !anyNA(d$id)) tables$descriptive_statistics$missing <<- length(subjects) - tables$descriptive_statistics$n
      missing_values <- sum(is.na(d$y) & !is.na(d$group))
      if (length(subjects)) missing_values <- length(subjects)*length(original_groups) - sum(complete.cases(d[, c("y", "group", "id")]))
      tables$data_summary <<- data.frame(treatments = length(original_groups), subjects = if (length(subjects)) length(subjects) else NA_integer_, missing_values = missing_values)
    }
    column_values <- if ("y" %in% names(d) && "group" %in% names(d) && method != "survival") d[!is.na(d$y) & !is.na(d$group), c("y", "group"), drop = FALSE] else NULL
    missing_rows <- sum(!complete.cases(d)); d <- droplevels(d[complete.cases(d), , drop = FALSE])
    assert(nrow(d) > 0, "No complete observations remain in the selected columns.")
    incomplete_subjects <- 0L; retained_incomplete_subjects <- 0L
    if ("id" %in% names(d)) {
      observed_groups <- original_groups
      all_ids <- unique(as.character(group_column(data, cfg$id))); all_ids <- all_ids[!is.na(all_ids)]
      coverage <- table(d$id, factor(d$group, levels = observed_groups))
      keep <- rownames(coverage)[rowSums(coverage == 1) == length(observed_groups)]
      incomplete_subjects <- length(setdiff(all_ids, keep))
      if (method == "mixed") {
        retained_incomplete_subjects <- sum(!rownames(coverage) %in% keep)
        assert(identical(levels(d$group), original_groups), "Every condition needs observed data for the mixed-effects model.")
        assert(nlevels(d$id) >= 3 && sum(rowSums(coverage) >= 2) >= 2, "Mixed-effects analysis requires at least three subjects and repeated observations for at least two subjects.")
        # Conditions must overlap across subjects to identify a within-subject effect.
        connected <- 1L
        repeat {
          overlap <- crossprod(coverage > 0) > 0
          expanded <- union(connected, which(colSums(overlap[connected, , drop = FALSE]) > 0))
          if (length(expanded) == length(connected)) break
          connected <- expanded
        }
        assert(length(connected) == length(original_groups), "Conditions are disconnected across subjects; a repeated mixed effect is not identifiable.")
      } else {
        d <- droplevels(d[as.character(d$id) %in% keep, , drop = FALSE])
        assert(nlevels(d$id) >= 2, "At least two complete matched subjects are required.")
      }
      d <- d[order(d$id, d$group), , drop = FALSE]
    }
    tables$data_audit <<- data.frame(input_rows = nrow(data), missing_rows = missing_rows,
      incomplete_subjects = incomplete_subjects, retained_incomplete_subjects = retained_incomplete_subjects,
      used_rows = nrow(d), excluded_rows = nrow(data) - nrow(d))
    if (missing_rows) notes <<- c(notes, paste(missing_rows, "rows excluded for missing values in selected columns."))
    if (incomplete_subjects) notes <<- c(notes, if (method == "mixed") paste(retained_incomplete_subjects, "incomplete subjects retained;", incomplete_subjects - retained_incomplete_subjects, "subjects with no usable observations excluded.") else paste(incomplete_subjects, "incomplete matched subjects excluded."))
    raw_d <- d
    if (cfg$scale == "log") {
      assert(method %in% c("rm_anova", "mixed"), "Log response models are available for repeated ANOVA and mixed effects.")
      assert(all(column_values$y > 0), "Lognormal analysis requires every observed response to be positive; zeros and negatives cannot be logged.")
      d$y <- log(d$y)
      tables$original_summary <<- summary_table(raw_d$y, raw_d$group, cfg$confidence)
      notes <<- c(notes, "Model and diagnostic differences use natural-log responses. Comparisons exponentiate log differences to geometric-mean ratios. The raw-data graph uses original response units.")
    }
    group_methods <- c("welch", "student", "paired_t", "mann_whitney", "paired_wilcox", "anova", "welch_anova", "kruskal", "friedman", "rm_anova", "mixed", "two_way")
    if (method %in% group_methods) {
      assert(nlevels(d$group) >= 2, "At least two groups are required.")
      assert(all(table(d$group) >= 2), "At least two complete observations per group are required.")
    }
    if ("y" %in% names(d) && "group" %in% names(d) && method != "survival") tables$summary <<- summary_table(d$y, d$group, cfg$confidence)
    if (cfg$scale == "log") {
      s <- tables$summary
      tables$geometric_summary <<- data.frame(group = s$group, n = s$n, geometric_mean = exp(s$mean), geometric_sd = exp(s$sd), ci_low = exp(s$ci_low), ci_high = exp(s$ci_high))
    }
    if (method == "outliers") {
      assert(length(cfg$rout_q) == 1 && is.finite(cfg$rout_q) && cfg$rout_q > 0 && cfg$rout_q < 1, "ROUT Q must lie between 0 and 1.")
      detections <- lapply(split(d, d$group, drop = TRUE), function(z) {
        r <- rout_column(z$y, cfg$rout_q)
        data.frame(group = as.character(z$group), source_row = z$source_row, value = z$y,
          robust_mean = r$location, RSDR = r$scale, residual = r$residual, p_value = r$p_value,
          threshold = r$threshold, outlier = r$outlier, status = r$status)
      })
      detections <- do.call(rbind, detections); rownames(detections) <- NULL
      tables$outlier_details <<- detections
      tables$outliers <<- detections[detections$outlier, , drop = FALSE]
      tables$outlier_summary <<- do.call(rbind, lapply(split(detections, detections$group), function(z) data.frame(group = z$group[1], method = "ROUT constant fit", Q = cfg$rout_q,
        n = nrow(z), outliers = sum(z$outlier), robust_mean = z$robust_mean[1], RSDR = z$RSDR[1], status = z$status[1])))
      cleaned <- data; cleaned[detections$source_row[detections$outlier], cfg$y] <- NA_real_
      tables$cleaned_data <<- cleaned
      notes <<- c(notes, "Outliers are identified, not removed from the original input. The cleaned-data sheet blanks flagged responses while preserving rows and subject matching. Q is an FDR target per column, not a family correction across all columns. Constant or too-small samples cannot be assessed.")
    }
    if (method %in% c("welch", "student", "paired_t", "mann_whitney", "paired_wilcox")) {
      assert(nlevels(d$group) == 2, "This test requires exactly two groups.")
      values <- split(d$y, d$group); is_paired <- method %in% c("paired_t", "paired_wilcox")
      if (method %in% c("welch", "student", "paired_t")) {
        test <- t.test(values[[1]], values[[2]], paired = is_paired, var.equal = method == "student", conf.level = cfg$confidence)
        tables$test <<- htest_table(test)
        tables$test$estimate <<- mean(values[[1]]) - mean(values[[2]])
      } else tables$test <<- htest_table(wilcox.test(values[[1]], values[[2]], paired = is_paired, exact = FALSE, correct = TRUE))
      tables$test$contrast <<- paste(names(values), collapse = " minus ")
    }
    if (method == "anova") {
      fit <<- aov(y ~ group, data = d)
      a <- as.data.frame(summary(fit)[[1]]); tables$anova <<- data.frame(term = trimws(rownames(a)), a, row.names = NULL, check.names = FALSE)
      if (cfg$adjustment == "tukey" && cfg$comparisons != "none") {
        assert(cfg$comparisons == "all", "Tukey supplies all pairwise comparisons. Select Sidak or another correction for control comparisons.")
        p <- as.data.frame(TukeyHSD(fit, conf.level = cfg$confidence)$group)
        tables$posthoc <<- data.frame(contrast = rownames(p), p, row.names = NULL, check.names = FALSE)
        if (isTRUE(cfg$reverse_comparisons)) {
          tables$posthoc$diff <<- -tables$posthoc$diff
          low <- tables$posthoc$lwr; tables$posthoc$lwr <<- -tables$posthoc$upr; tables$posthoc$upr <<- -low
          pair_names <- combn(levels(d$group), 2, simplify = FALSE)
          tables$posthoc$contrast <<- vapply(pair_names, paste, character(1), collapse = " minus ")
        }
        pairs <- combn(levels(d$group), 2, simplify = FALSE)
        if (!isTRUE(cfg$reverse_comparisons)) pairs <- lapply(pairs, rev)
        g1 <- vapply(pairs, `[`, character(1), 1); g2 <- vapply(pairs, `[`, character(1), 2)
        means <- tapply(d$y, d$group, mean); counts <- table(d$group)
        tables$posthoc$group_1 <<- g1; tables$posthoc$group_2 <<- g2
        tables$posthoc$mean_1 <<- unname(means[g1]); tables$posthoc$mean_2 <<- unname(means[g2]); tables$posthoc$mean_type <<- "Observed"
        tables$posthoc$n1 <<- as.numeric(counts[g1]); tables$posthoc$n2 <<- as.numeric(counts[g2])
        tables$posthoc$se <<- sqrt(sigma(fit)^2*(1/as.numeric(counts[g1]) + 1/as.numeric(counts[g2])))
        tables$posthoc$df <<- df.residual(fit); tables$posthoc$statistic <<- tables$posthoc$diff/tables$posthoc$se
      } else tables$posthoc <<- mean_comparisons(d, cfg, kind = "ordinary", fit = fit)
    }
    if (method == "welch_anova") {
      tables$test <<- htest_table(oneway.test(y ~ group, data = d, var.equal = FALSE))
      tables$test$denominator_df <<- unname(oneway.test(y ~ group, data = d)$parameter[2])
      tables$posthoc <<- mean_comparisons(d, cfg, kind = "welch")
    }
    if (method == "kruskal") {
      tables$test <<- htest_table(kruskal.test(y ~ group, data = d))
      tables$posthoc <<- rank_comparisons(d, cfg)
    }
    if (method == "friedman") {
      tables$test <<- htest_table(friedman.test(y ~ group | id, data = d))
      tables$posthoc <<- rank_comparisons(d, cfg, matched = TRUE)
    }
    if (method == "rm_anova") {
      model <- repeated_anova(d, cfg$sphericity); fit <<- model$fit
      tables$test <<- model$tests
      tables$posthoc <<- mean_comparisons(d, cfg, kind = if (cfg$sphericity == "gg") "paired" else "rm_pooled", fit = fit)
      a <- anova(fit)
      tables$matching <<- data.frame(method = "Subject fixed-effect F test", statistic = a["id", "F value"], numerator_df = a["id", "Df"], denominator_df = a["Residuals", "Df"], p_value = a["id", "Pr(>F)"])
      notes <<- c(notes, paste("Selected omnibus test:", if (cfg$sphericity == "gg") "Geisser-Greenhouse corrected." else "Sphericity assumed."))
    }
    if (method == "mixed") {
      assert(requireNamespace("nlme", quietly = TRUE), "Install the nlme package using setup.R.")
      assert(nrow(d) > nlevels(d$id) + nlevels(d$group) - 1, "Not enough observations to estimate subject and residual variation.")
      fit <<- fit_mixed_options(d, cfg)
      a <- mixed_omnibus(fit); info <- attr(fit, "labstats")
      tables$test <<- data.frame(method = if (info$removed) "REML independent model after removing subject effect" else "REML compound-symmetry model: condition effect", statistic = unname(a["statistic"]),
        numerator_df = unname(a["df1"]), denominator_df = unname(a["df2"]), p_value = pf(a["statistic"], a["df1"], a["df2"], lower.tail = FALSE))
      uncorrected <- tables$test
      gg <- if (cfg$sphericity == "gg") mixed_epsilon(d) else list(epsilon = 1, source = "Sphericity assumed", iterations = 0)
      tables$sphericity <<- data.frame(epsilon = gg$epsilon, estimator = gg$source, iterations = gg$iterations)
      tables$test$uncorrected_p_value <<- uncorrected$p_value
      tables$test$numerator_df <<- uncorrected$numerator_df * gg$epsilon
      tables$test$denominator_df <<- uncorrected$denominator_df * gg$epsilon
      tables$test$epsilon <<- gg$epsilon
      tables$test$correction <<- if (cfg$sphericity == "gg") "Geisser-Greenhouse" else "Sphericity assumed"
      tables$test$p_value <<- pf(uncorrected$statistic, uncorrected$numerator_df * gg$epsilon, uncorrected$denominator_df * gg$epsilon, lower.tail = FALSE)
      tables$posthoc <<- mean_comparisons(d, cfg, kind = if (cfg$sphericity == "gg") "mixed_pairs" else "mixed_pooled", fit = fit)
      independent <- nlme::gls(y ~ group, data = d, method = "REML", na.action = na.fail)
      likelihood_ratio <- max(0, 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(independent))))
      tables$matching <<- data.frame(method = "REML likelihood ratio: subject intercept vs independent residuals", statistic = likelihood_ratio,
        df = 1, p_value = pchisq(likelihood_ratio, 1, lower.tail = FALSE), boundary_mixture_p_value = if (cfg$random_effect == "full") NA_real_ else if (likelihood_ratio == 0) 1 else .5 * pchisq(likelihood_ratio, 1, lower.tail = FALSE))
      notes <<- c(notes, if (cfg$random_effect == "full") "Matching diagnostic uses a conventional chi-square(1) likelihood-ratio test for signed within-subject covariance. A boundary-mixture p-value does not apply to the signed model." else "Matching diagnostic reports the conventional chi-square(1) likelihood-ratio p-value and the boundary-adjusted 50:50 chi-square(0)/chi-square(1) p-value for a nonnegative subject component. Fixed effects and input observations are identical in both REML fits.",
        paste("GG epsilon estimator:", gg$source, ". Missing-data covariance estimation does not fill or change the input table. Exact Prism equivalence has not been validated."))
      coefficients <- mixed_coefficients(fit); critical <- qt((1 + cfg$confidence)/2, coefficients[, "DF"])
      tables$coefficients <<- data.frame(term = rownames(coefficients), estimate = coefficients[, "Value"], se = coefficients[, "Std.Error"],
        df = coefficients[, "DF"], p_value = coefficients[, "p-value"], ci_low = coefficients[, "Value"] - critical * coefficients[, "Std.Error"],
        ci_high = coefficients[, "Value"] + critical * coefficients[, "Std.Error"], row.names = NULL)
      design <- model.matrix(~group, data.frame(group = factor(levels(d$group), levels = levels(d$group))))
      tables$predicted_means <<- data.frame(group = levels(d$group), n = as.numeric(table(d$group)),
        mean = as.numeric(design %*% mixed_beta(fit)), se = sqrt(diag(design %*% vcov(fit) %*% t(design))))
      v <- c(info$subject_variance, info$residual_variance)
      tables$variance_components <<- data.frame(component = c("Subject covariance component", "Residual"), variance = v, sd = sqrt(pmax(v, 0)),
        note = c(if (v[1] < 0) "Signed covariance component, not a negative physical variance; SD not defined" else if (info$removed) "Removed and model refitted" else "Subject random intercept", "Common residual variance"))
      if (v[1] < 0) tables$variance_components$sd[1] <<- NA_real_
      tables$model <<- data.frame(observations = nrow(d), subjects = nlevels(d$id), estimation = "REML", covariance = if (info$removed) "Independent residuals" else "Compound symmetry", AIC = AIC(fit),
        random_effect_policy = cfg$random_effect, GLM_initial_values = cfg$glm_initial, starts_tried = info$starts_tried, selected_start = info$selected_start)
      if (v[1] < 1e-6*v[2]) warnings <<- c(warnings, if (info$removed) "Subject component was nonpositive; removed and the simpler model refitted." else "The subject covariance component is near zero or negative; review the model summary.")
      if (cfg$random_effect == "full") notes <<- c(notes, "Full model uses REML compound symmetry, allowing a signed subject covariance component while keeping the observation covariance positive definite. Denominator df use the repeated-measures containment convention, not Satterthwaite. Negative components are not physical variances.")
    }
    if (cfg$comparisons == "trend") tables$trend <<- trend_table(d, cfg, fit, method)
    if (method %in% c("anova", "rm_anova", "mixed")) tables <<- c(tables, extra_model_tables(d, fit, cfg, method))
    if (method == "two_way") {
      cells <- table(d$group, d$factor2)
      assert(nlevels(d$factor2) >= 2 && all(cells >= 2) && length(unique(as.vector(cells))) == 1, "Two-way ANOVA requires at least two replicates per cell and equal cell sizes in a complete factorial design.")
      fit <<- aov(y ~ group * factor2, data = d)
      a <- as.data.frame(summary(fit)[[1]]); tables$anova <<- data.frame(term = trimws(rownames(a)), a, row.names = NULL, check.names = FALSE)
    }
    if (method %in% c("pearson", "spearman", "linear")) {
      assert(nrow(d) >= 3 && length(unique(d$x)) >= 2 && length(unique(d$y)) >= 2, "At least three complete rows and variation in X and Y are required.")
      if (method != "linear") tables$test <<- htest_table(cor.test(d$x, d$y, method = method, exact = FALSE, conf.level = cfg$confidence))
      else {
        fit <<- lm(y ~ x, data = d); s <- summary(fit); ci <- confint(fit, level = cfg$confidence)
        tables$coefficients <<- data.frame(term = rownames(s$coefficients), s$coefficients, ci_low = ci[, 1], ci_high = ci[, 2], row.names = NULL, check.names = FALSE)
        tables$model <<- data.frame(n = nrow(d), r_squared = s$r.squared, adjusted_r_squared = s$adj.r.squared, residual_sd = s$sigma)
      }
    }
    if (method == "dose") {
      assert(all(d$x > 0), "Dose-response fitting requires positive concentrations. Exclude zero-dose controls or use a model that explicitly supports them.")
      assert(nrow(d) >= 6 && length(unique(d$x)) >= 5 && diff(range(d$y)) > 0, "Use at least six observations, five distinct positive concentrations, and a varying response.")
      assert(cfg$direction %in% c("increasing", "decreasing"), "Select an increasing or decreasing curve.")
      sign_hill <- if (cfg$direction == "decreasing") 1 else -1
      candidates <- lapply(c(.25, .5, .75), function(q) tryCatch(nls(
        y ~ bottom + exp(logSpan) / (1 + exp(hill * (log(x) - logEC50))), data = d,
        start = list(bottom = min(d$y), logSpan = log(diff(range(d$y))), logEC50 = as.numeric(quantile(log(d$x), q)), hill = sign_hill),
        algorithm = "port", lower = c(-Inf, -30, log(min(d$x)) - 10, if (sign_hill > 0) .01 else -10),
        upper = c(Inf, 30, log(max(d$x)) + 10, if (sign_hill > 0) 10 else -.01),
        control = nls.control(maxiter = 500)), error = function(e) NULL))
      candidates <- Filter(function(f) !is.null(f) && isTRUE(f$convInfo$isConv), candidates)
      assert(length(candidates) > 0, "The 4PL model did not converge. Check curve direction and whether doses span both plateaus.")
      fit <<- candidates[[which.min(vapply(candidates, deviance, numeric(1)))]]
      b <- coef(fit); v <- vcov(fit); critical <- qt((1 + cfg$confidence)/2, df.residual(fit))
      se <- sqrt(v["logEC50", "logEC50"])
      tables$parameters <<- data.frame(parameter = c("Bottom", "Top", "EC50", "Hill slope"), estimate = c(b["bottom"], b["bottom"] + exp(b["logSpan"]), exp(b["logEC50"]), b["hill"]), row.names = NULL)
      tables$ec50_interval <<- data.frame(EC50 = exp(b["logEC50"]), ci_low = exp(b["logEC50"] - critical * se), ci_high = exp(b["logEC50"] + critical * se), confidence = cfg$confidence)
      tables$model <<- data.frame(n = nrow(d), residual_sd = sigma(fit), residual_df = df.residual(fit))
      if (exp(b["logEC50"]) < min(d$x) || exp(b["logEC50"]) > max(d$x)) warnings <<- c(warnings, "Fitted EC50 lies outside measured doses: extrapolated estimate.")
      if (kappa(v) > 1e8) warnings <<- c(warnings, "The curve parameters are poorly identified; collect data spanning both plateaus before relying on EC50.")
      if (abs(b["hill"]) < .011 || abs(b["hill"]) > 9.99) warnings <<- c(warnings, "Hill slope reached a fitting bound. The curve may not be supported by these data.")
    }
    if (method %in% c("chisq", "fisher")) {
      tab <- xtabs(count ~ group + category, data = d)
      assert(nrow(tab) >= 2 && ncol(tab) >= 2 && all(rowSums(tab) > 0) && all(colSums(tab) > 0), "A table of at least 2 by 2 with nonzero row and column totals is required.")
      test <- if (method == "chisq") chisq.test(tab, correct = FALSE) else fisher.test(tab, conf.level = cfg$confidence)
      tables$test <<- htest_table(test)
      tables$counts <<- data.frame(group = rownames(tab), as.data.frame.matrix(tab), row.names = NULL, check.names = FALSE)
      if (method == "chisq") {
        tables$expected <<- data.frame(group = rownames(tab), as.data.frame.matrix(test$expected), row.names = NULL, check.names = FALSE)
        if (any(test$expected < 5)) warnings <<- c(warnings, "Some expected counts are below 5. The chi-square approximation may be unreliable; consider Fisher's exact test.")
      }
    }
    if (method == "survival") {
      assert(requireNamespace("survival", quietly = TRUE), "Install the survival package.")
      assert(sum(d$event) > 0, "At least one event is required.")
      fit <<- survival::survfit(survival::Surv(y, event) ~ group, data = d, conf.int = cfg$confidence)
      if (nlevels(d$group) >= 2) {
        test <- survival::survdiff(survival::Surv(y, event) ~ group, data = d)
        tables$test <<- data.frame(method = "Log-rank", statistic = test$chisq, df = length(test$n) - 1, p_value = pchisq(test$chisq, length(test$n) - 1, lower.tail = FALSE))
      }
      s <- summary(fit)
      tables$survival <<- data.frame(time = s$time, n_risk = s$n.risk, n_event = s$n.event, survival = s$surv, ci_low = s$lower, ci_high = s$upper, group = if (is.null(s$strata)) rep(levels(d$group)[1], length(s$time)) else sub("^group=", "", as.character(s$strata)))
    }
    if (!method %in% c("chisq", "fisher", "survival", "outliers")) {
      tables$normality <<- normality_table(d, method, fit)
      notes <<- c(notes, "Shapiro-Wilk diagnostics are exploratory and unadjusted. A large p-value does not establish normality; inspect distributions/QQ plots and the experimental design. The diagnostic does not select or change the analysis.")
    }
    if (!is.null(column_values)) {
      tables$distribution_tests <<- distribution_tests(column_values$y, column_values$group, cfg$alpha)
      tables$distribution_counts <<- do.call(rbind, lapply(split(column_values$y, column_values$group), function(x) data.frame(n = length(x), zeros = sum(x == 0), negatives = sum(x < 0))))
      tables$distribution_counts <<- data.frame(group = rownames(tables$distribution_counts), tables$distribution_counts, row.names = NULL)
    }
    if (method %in% c("anova", "welch_anova", "kruskal", "rm_anova", "mixed") && cfg$adjustment != "tukey") notes <<- c(notes, paste("Pairwise p-value adjustment:", cfg$adjustment, "; comparison family:", cfg$comparisons, ".", if (method == "kruskal") "Rank comparisons supply adjusted p-values without confidence intervals." else "Intervals use Sidak or Bonferroni simultaneous limits; Holm/BH corrections retain clearly labelled unadjusted intervals."))
    if (method %in% c("anova", "welch_anova", "kruskal", "rm_anova", "mixed") && cfg$adjustment == "sidak") notes <<- c(notes, "Sidak familywise interpretation assumes independent comparisons. Shared-control and paired contrasts are correlated, so familywise coverage is approximate; Bonferroni does not require comparison independence.")
    residuals_for_tests <- residual_frame(list(method = method, data = d, fit = fit))
    if (isTRUE(cfg$residual_normality)) {
      assert(!is.null(residuals_for_tests), "Residual normality is unavailable for this analysis.")
      tables$residual_normality <<- residual_normality_tests(residuals_for_tests$residual)
    }
    if (isTRUE(cfg$residual_variance)) {
      assert(!method %in% c("rm_anova", "mixed", "friedman", "paired_t", "paired_wilcox"), "Independent-group variance tests are unavailable for matched data.")
      assert(!is.null(residuals_for_tests), "Residual variance diagnostics are unavailable for this analysis.")
      tables$residual_variance <<- residual_variance_tests(residuals_for_tests)
    }
    tables <- lapply(tables, annotate_table, alpha = cfg$alpha)
    if (!is.null(tables$posthoc) && cfg$adjustment %in% c("BH", "BY", "BKY")) {
      tables$posthoc$discovery <- if (cfg$adjustment == "BKY") bky_rejections(tables$posthoc$p_value, cfg$fdr_q) else tables$posthoc$adjusted_p_value <= cfg$fdr_q
      tables$posthoc$target_FDR <- cfg$fdr_q
      tables$posthoc$significant <- NULL; tables$posthoc$p_summary <- NULL
      notes <- c(notes, "FDR comparisons report discoveries at the selected Q. Any comparison intervals are unadjusted, not simultaneous FDR confidence intervals.")
    }
    list(method = method, config = cfg, tables = tables, fit = fit, data = d, raw_data = raw_d, column_data = column_values, notes = notes)
  }
  result <- withCallingHandlers(calc(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  result$warnings <- unique(warnings)
  result
}

analysis_report <- function(result, options = list(), worksheet = NULL) {
  settings <- modifyList(list(output_digits = 4L, p_style = "prism", p_digits = 6L), modifyList(result$config, options))
  sheets <- result_sheets(result, worksheet = worksheet)
  c("LabStats analysis report", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    names(analysis_methods)[match(result$method, analysis_methods)], "", result$notes,
    if (length(result$warnings)) c("", "WARNINGS", result$warnings),
    unlist(lapply(sheets, function(sheet) c("", toupper(sheet$title), capture.output(print(sheet_frame(sheet, settings), row.names = FALSE)), sheet$note))),
    unlist(lapply(names(attr(sheets, "additional")), function(name) {
      table <- attr(sheets, "additional")[[name]]
      c("", toupper(name), capture.output(print(format_result_table(table, settings$output_digits, settings$p_style, settings$p_digits), row.names = FALSE)))
    })))
}
