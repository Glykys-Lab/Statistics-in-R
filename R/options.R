# Shared calculations for the one-way analysis options. Embedded in R exports.
comparison_pairs <- function(groups, cfg) {
  if (cfg$comparisons %in% c("none", "trend")) return(list())
  if (cfg$comparisons == "control") {
    assert(length(cfg$control) == 1 && cfg$control %in% groups, "Select a valid control group.")
    pairs <- lapply(setdiff(groups, cfg$control), function(g) c(cfg$control, g))
  } else if (cfg$comparisons == "selected") {
    indices <- as.character(unlist(cfg$selected_pairs))
    assert(length(indices) > 0, "Select at least one pair of columns.")
    pairs <- lapply(indices, function(index) {
      pair <- vapply(strsplit(index, "::", fixed = TRUE)[[1]], utils::URLdecode, character(1))
      assert(length(pair) == 2 && all(pair %in% groups) && pair[1] != pair[2], "A selected pair is no longer valid for this table or has no observed data. Select pairs again.")
      pair
    })
    assert(!anyDuplicated(indices), "Select each comparison only once.")
  } else pairs <- combn(groups, 2, simplify = FALSE)
  if (isTRUE(cfg$reverse_comparisons)) pairs <- lapply(pairs, rev)
  pairs
}

bky_rejections <- function(p, q, n = length(p)) {
  finite <- which(is.finite(p)); values <- sort(p[finite]); ranks <- seq_along(values)
  first <- which(values <= ranks * q/(1 + q)/n)
  r1 <- if (length(first)) max(first) else 0L
  r <- if (r1 == 0L || r1 == n) r1 else {
    second <- which(values <= ranks * q/(1 + q)/(n - r1))
    if (length(second)) max(second) else 0L
  }
  cutoff <- if (r) values[r] else -Inf
  !is.na(p) & p <= cutoff
}

adjust_family <- function(p, adjustment, n = length(p)) {
  if (adjustment == "sidak") return(-expm1(n * log1p(-p)))
  if (adjustment == "holm_sidak") {
    ii <- order(p, na.last = NA); out <- rep(NA_real_, length(p))
    out[ii] <- cummax(-expm1((n - seq_along(ii) + 1) * log1p(-p[ii])))
    return(out)
  }
  if (adjustment != "BKY") return(p.adjust(p, method = adjustment, n = n))
  # Adjusted q = smallest target FDR at which the two-stage procedure rejects.
  vapply(seq_along(p), function(i) {
    if (!is.finite(p[i])) return(NA_real_)
    low <- 0; high <- 1
    for (j in seq_len(50)) {
      q <- (low + high)/2
      if (bky_rejections(p, q, n)[i]) high <- q else low <- q
    }
    high
  }, numeric(1))
}

fit_mixed_options <- function(d, cfg, method = "REML", formula = y ~ group) {
  policy <- cfg$random_effect
  if (is.null(policy)) policy <- "keep"
  glm_start <- !identical(cfg$glm_initial, FALSE)
  attempts <- list(); labels <- character(); attempted <- 0L
  add_fit <- function(fit, label) {
    attempted <<- attempted + 1L
    value <- tryCatch(fit(), error = function(e) NULL)
    if (!is.null(value) && is.finite(as.numeric(logLik(value)))) {
      attempts[[length(attempts) + 1L]] <<- value; labels <<- c(labels, label)
    }
  }
  a <- tryCatch(anova(lm(y ~ id + group, data = d)), error = function(e) NULL)
  residual <- if (!is.null(a)) a["Residuals", "Mean Sq"] else var(d$y)
  tau <- if (!is.null(a)) (a["id", "Mean Sq"] - residual)/mean(table(d$id)) else .1
  if (!is.finite(tau) || !is.finite(residual) || residual <= 0) { tau <- .1; residual <- 1 }
  starts <- .1
  if (glm_start) starts <- c(tau/residual, starts)
  if (policy == "keep") {
    for (s in unique(pmax(starts, 1e-5))) add_fit(function() nlme::lme(formula,
      random = list(id = nlme::pdDiag(matrix(s, 1, 1), form = ~1)), data = d,
      method = method, na.action = na.fail, control = nlme::lmeControl(msMaxIter = 200)), paste("variance ratio", signif(s, 4)))
  } else {
    lower <- -1/(max(table(d$id)) - 1) + 1e-5
    rhos <- unique(pmin(.99, pmax(lower, starts/(1 + starts))))
    for (rho in rhos) add_fit(function() nlme::gls(formula, data = d, method = method,
      correlation = nlme::corCompSymm(value = rho, form = ~1 | id),
      na.action = na.fail, control = nlme::glsControl(msMaxIter = 200)), paste("correlation", signif(rho, 4)))
  }
  assert(length(attempts) > 0, "The mixed model did not converge with the chosen initial values.")
  best <- which.max(vapply(attempts, function(f) as.numeric(logLik(f)), numeric(1)))
  fit <- attempts[[best]]; removed <- FALSE
  if (inherits(fit, "lme")) {
    v <- nlme::VarCorr(fit); tau <- as.numeric(v[1, "Variance"]); residual <- as.numeric(v[2, "Variance"])
  } else {
    rho <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
    tau <- sigma(fit)^2 * rho; residual <- sigma(fit)^2 * (1 - rho)
    if (policy == "remove" && tau <= 1e-7 * residual) {
      fit <- nlme::gls(formula, data = d, method = method, na.action = na.fail)
      tau <- 0; residual <- sigma(fit)^2; removed <- TRUE
    }
  }
  fit$data <- d
  p <- length(coef(if (inherits(fit, "lme")) lm(formula, d) else fit))
  attr(fit, "labstats") <- list(policy = policy, removed = removed, subject_variance = tau,
    residual_variance = residual, df = if (removed) nrow(d) - p else nrow(d) - nlevels(d$id) - (p - 1),
    glm_initial = glm_start, starts_tried = attempted, starts_converged = length(attempts), selected_start = labels[best])
  fit
}

mixed_coefficients <- function(fit) {
  if (inherits(fit, "lme")) return(summary(fit)$tTable)
  info <- attr(fit, "labstats"); beta <- coef(fit); se <- sqrt(diag(vcov(fit)))
  df <- rep(info$df, length(beta))
  if (!info$removed) df[names(beta) == "(Intercept)"] <- nlevels(fit$data$id) - 1
  cbind(Value = beta, Std.Error = se, DF = df, `t-value` = beta/se, `p-value` = 2*pt(-abs(beta/se), df))
}

mixed_beta <- function(fit) if (inherits(fit, "lme")) nlme::fixef(fit) else coef(fit)

mixed_omnibus <- function(fit) {
  if (inherits(fit, "lme")) {
    a <- anova(fit, type = "marginal")
    return(c(statistic = a["group", "F-value"], df1 = a["group", "numDF"], df2 = a["group", "denDF"]))
  }
  beta <- coef(fit); ii <- which(names(beta) != "(Intercept)"); v <- vcov(fit)[ii, ii, drop = FALSE]
  c(statistic = as.numeric(crossprod(beta[ii], solve(v, beta[ii])))/length(ii), df1 = length(ii), df2 = attr(fit, "labstats")$df)
}

rank_comparisons <- function(d, cfg, matched = FALSE) {
  pairs <- comparison_pairs(levels(d$group), cfg)
  if (!length(pairs)) return(NULL)
  scores <- if (matched) ave(d$y, d$id, FUN = rank) else rank(d$y)
  means <- tapply(scores, d$group, mean); n <- table(d$group); N <- nrow(d); k <- nlevels(d$group)
  ties <- if (matched) sum(vapply(split(d$y, d$id), function(x) sum(table(x)^3 - table(x)), numeric(1))) else sum(table(d$y)^3 - table(d$y))
  variance <- if (matched) k*(k + 1)/(6*nlevels(d$id)) * (1 - ties/(nlevels(d$id)*(k^3 - k))) else N*(N + 1)/12 - ties/(12*(N - 1))
  result <- do.call(rbind, lapply(pairs, function(p) {
    a <- d[d$group == p[1], ]; b <- d[d$group == p[2], ]
    se <- if (matched) sqrt(variance) else sqrt(variance*(1/n[p[1]] + 1/n[p[2]]))
    z <- (means[p[1]] - means[p[2]])/se
    if (cfg$rank_test == "wilcox") {
      if (matched) b <- b[match(a$id, b$id), ]
      t <- wilcox.test(a$y, b$y, paired = matched, exact = FALSE, correct = TRUE)
      z <- unname(t$statistic); probability <- t$p.value
    } else probability <- if (is.finite(z)) 2*pnorm(-abs(z)) else NA_real_
    data.frame(group_1 = p[1], group_2 = p[2], contrast = paste(p, collapse = " minus "),
      mean_rank_difference = unname(means[p[1]] - means[p[2]]), statistic = unname(z), p_value = probability,
      test = if (cfg$rank_test == "wilcox") "Pairwise Wilcoxon (normal approximation)" else "Dunn rank comparison with tie correction")
  }))
  result$adjusted_p_value <- adjust_family(result$p_value, cfg$adjustment)
  result
}

trend_table <- function(d, cfg, fit, method) {
  assert(method %in% c("anova", "rm_anova", "mixed"), "Linear trend requires ordinary or repeated parametric ANOVA/mixed analysis.")
  assert(method == "anova" || cfg$sphericity == "none", "Linear trend requires sphericity assumed for repeated data; it is unavailable with Geisser-Greenhouse correction.")
  k <- nlevels(d$group); n <- as.numeric(table(d$group)); x <- seq_len(k)
  if (method == "anova") {
    center <- sum(n*x)/sum(n); weights <- n*(x - center)/sum(n*(x - center)^2)
  } else weights <- (x - mean(x))/sum((x - mean(x))^2)
  design <- model.matrix(~group, data.frame(group = factor(levels(d$group), levels = levels(d$group))))
  if (method == "mixed") {
    contrast <- as.numeric(weights %*% design); estimate <- sum(contrast*mixed_beta(fit))
    se <- sqrt(as.numeric(contrast %*% vcov(fit) %*% contrast)); df <- min(mixed_coefficients(fit)[contrast != 0, "DF"])
  } else {
    means <- tapply(d$y, d$group, mean); estimate <- sum(weights*means)
    se <- sqrt(sigma(fit)^2*sum(weights^2/n)); df <- df.residual(fit)
  }
  critical <- qt((1 + cfg$confidence)/2, df)
  data.frame(test = "Linear trend across equally spaced column positions", slope = estimate, se = se,
    statistic = (estimate/se)^2, numerator_df = 1, denominator_df = df,
    p_value = pf((estimate/se)^2, 1, df, lower.tail = FALSE), ci_low = estimate - critical*se,
    ci_high = estimate + critical*se, column_sequence = paste(levels(d$group), collapse = " -> "))
}

residual_frame <- function(result) {
  d <- result$data; fit <- result$fit; method <- result$method
  if (method %in% c("survival", "chisq", "fisher", "pearson", "spearman", "outliers")) return(NULL)
  if (!is.null(fit)) {
    predicted <- as.numeric(fitted(fit)); r <- as.numeric(residuals(fit))
    if (inherits(fit, "gls") && !is.null(attr(fit, "labstats"))) {
      info <- attr(fit, "labstats")
      if (!info$removed && info$subject_variance > 0) {
        corrections <- tapply(r, d$id, sum)*info$subject_variance/(info$residual_variance + as.numeric(table(d$id))*info$subject_variance)
        predicted <- predicted + corrections[as.character(d$id)]; r <- d$y - predicted
      }
    }
  } else if ("group" %in% names(d) && method %in% c("welch_anova", "kruskal", "friedman", "summary", "normality")) {
    fun <- if (method %in% c("kruskal", "friedman")) median else mean
    predicted <- ave(d$y, d$group, FUN = fun)
    if (method == "friedman") predicted <- predicted + ave(d$y, d$id, FUN = median) - median(d$y)
    r <- d$y - predicted
  } else return(NULL)
  frame <- data.frame(fitted = as.numeric(predicted), residual = as.numeric(r), observation = seq_len(nrow(d)))
  if ("group" %in% names(d)) frame$group <- d$group else frame$group <- factor("All")
  frame$subject <- if ("id" %in% names(d)) as.character(d$id) else as.character(ave(seq_len(nrow(d)), frame$group, FUN = seq_along))
  frame
}

residual_normality_tests <- function(x) {
  assert(requireNamespace("nortest", quietly = TRUE), "Install nortest with source('setup.R') for residual normality tests.")
  checks <- list(
    `D'Agostino-Pearson` = function() { p <- dagostino_pearson(x); assert(is.finite(p["p_value"]), "Requires at least 8 varying residuals"); c(p["K2"], p["p_value"]) },
    `Anderson-Darling` = function() { t <- nortest::ad.test(x); c(t$statistic, t$p.value) },
    `Shapiro-Wilk` = function() { t <- shapiro.test(x); c(t$statistic, t$p.value) },
    `Kolmogorov-Smirnov (Lilliefors correction)` = function() { t <- nortest::lillie.test(x); c(t$statistic, t$p.value) })
  do.call(rbind, lapply(names(checks), function(name) {
    status <- "Exploratory: fitted residuals are not independent observations"
    t <- tryCatch({ assert(sd(x) > 0, "Residuals are constant"); checks[[name]]() }, error = function(e) { status <<- paste("Unavailable:", conditionMessage(e)); c(NA_real_, NA_real_) })
    data.frame(test = name, n = length(x), statistic = unname(t[1]), p_value = unname(t[2]), note = status)
  }))
}

residual_variance_tests <- function(frame) {
  do.call(rbind, lapply(c("Brown-Forsythe", "Bartlett"), function(name) {
    status <- "Exploratory test of group variance equality"
    t <- tryCatch({
      assert(nlevels(frame$group) >= 2 && all(table(frame$group) >= 2), "At least two residuals in each group required")
      if (name == "Bartlett") { b <- bartlett.test(residual ~ group, frame); c(unname(b$statistic), b$p.value) }
      else {
        z <- abs(frame$residual - ave(frame$residual, frame$group, FUN = median))
        a <- anova(lm(z ~ frame$group)); c(a[1, "F value"], a[1, "Pr(>F)"])
      }
    }, error = function(e) { status <<- paste("Unavailable:", conditionMessage(e)); c(NA_real_, NA_real_) })
    data.frame(test = name, statistic = unname(t[1]), p_value = unname(t[2]), note = status)
  }))
}

extra_model_tables <- function(d, fit, cfg, method) {
  tables <- list()
  if (is.null(fit)) return(tables)
  if (isTRUE(cfg$report_goodness)) {
    if (method == "mixed") {
      info <- attr(fit, "labstats"); fixed <- as.numeric(model.matrix(~group, d) %*% mixed_beta(fit))
      vf <- var(fixed); total <- vf + info$subject_variance + info$residual_variance
      tables$goodness_of_fit <- data.frame(metric = c("Marginal R squared", "Conditional R squared", "Log likelihood", "AIC"),
        value = c(if (info$subject_variance >= 0) vf/total else NA_real_, if (info$subject_variance >= 0) (vf + info$subject_variance)/total else NA_real_, as.numeric(logLik(fit)), AIC(fit)),
        note = c(rep(if (info$subject_variance < 0) "Variance-based R squared unavailable for a signed covariance component" else "Variance decomposition for Gaussian random-intercept model", 2), rep("REML fit", 2)))
    } else tables$goodness_of_fit <- data.frame(R_squared = 1 - sum(residuals(fit)^2)/sum((d$y - mean(d$y))^2), residual_SD = sigma(fit), residual_df = df.residual(fit))
  }
  if (isTRUE(cfg$report_effects)) {
    if (method == "mixed") {
      info <- attr(fit, "labstats"); b <- mixed_beta(fit)[-1]
      tables$effect_sizes <- data.frame(contrast = names(b), standardized_difference = unname(b)/sqrt(info$residual_variance),
        denominator = "REML residual SD; treatment difference vs first column (model scale)")
    } else {
      a <- anova(fit); ss <- a["group", "Sum Sq"]; error <- a["Residuals", "Sum Sq"]
      tables$effect_sizes <- data.frame(effect = "Condition", partial_eta_squared = ss/(ss + error), eta_squared = ss/sum((d$y - mean(d$y))^2))
    }
  }
  if (isTRUE(cfg$report_aicc)) {
    # Different fixed-effect models must be compared with ML, not REML.
    formulas <- if (method == "rm_anova") list(`Condition + subject` = y ~ id + group, `Subject only` = y ~ id) else list(`Condition means` = y ~ group, `Common mean` = y ~ 1)
    fits <- lapply(formulas, function(f) if (method == "mixed") fit_mixed_options(d, cfg, method = "ML", formula = f) else lm(f, d))
    tables$model_comparison <- do.call(rbind, lapply(seq_along(fits), function(i) {
      ll <- logLik(fits[[i]]); K <- attr(ll, "df"); n <- nrow(d)
      data.frame(model = names(fits)[i], estimation = "ML", observations = n, parameters = K,
        AIC = AIC(fits[[i]]), AICc = if (n > K + 1) AIC(fits[[i]]) + 2*K*(K + 1)/(n - K - 1) else NA_real_)
    }))
    aa <- tables$model_comparison$AICc
    if (all(is.finite(aa))) {
      tables$model_comparison$delta_AICc <- aa - min(aa)
      w <- exp(-.5*(aa - min(aa))); tables$model_comparison$Akaike_weight <- w/sum(w)
    }
  }
  tables
}

visible_result_tables <- function(result, options = result$config) {
  tables <- result$tables
  if (identical(options$report_descriptive, FALSE)) tables[c("summary", "original_summary", "geometric_summary", "descriptive_statistics")] <- NULL
  if (identical(options$report_adjusted_p, FALSE) && !is.null(tables$posthoc)) {
    tables$posthoc$adjusted_p_value <- NULL; tables$posthoc[["p adj"]] <- NULL
  }
  tables
}

format_result_table <- function(table, digits = 4L, p_style = "prism", p_digits = 6L) {
  assert(length(digits) == 1 && is.finite(digits) && digits %in% 2:12, "Choose 2 to 12 significant digits.")
  assert(length(p_digits) == 1 && is.finite(p_digits) && p_digits %in% 1:15, "Choose 1 to 15 p-value digits.")
  assert(p_style %in% c("prism", "decimal", "scientific", "stars"), "Choose a valid p-value style.")
  pcolumns <- intersect(names(table), c("p_value", "adjusted_p_value", "p adj", "Pr(>F)", "Pr(>|t|)", "uncorrected_p_value", "boundary_mixture_p_value"))
  for (name in pcolumns) {
    p <- table[[name]]
    table[[name]] <- if (p_style == "stars") significance_summary(p) else if (p_style == "scientific") ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = p_digits - 1L)) else if (p_style == "decimal") ifelse(is.na(p), NA_character_, formatC(p, format = "f", digits = p_digits)) else ifelse(is.na(p), NA_character_, paste(format.pval(p, digits = 4, eps = .0001), significance_summary(p)))
  }
  for (name in setdiff(names(table), pcolumns)) if (is.numeric(table[[name]])) table[[name]] <- vapply(table[[name]], function(x) if (is.na(x)) NA_character_ else format(signif(x, digits), digits = digits, trim = TRUE), character(1))
  table
}
