# Graph construction for analyzed results: data plots, QQ/residual diagnostics,
# comparisons, ranks, matched differences and saved-graph layouts.
# make_graph() is the shared dispatcher for display and export. Nested analyses
# supply their own model-aware plots; factor levels preserve experimental order.
make_plot <- function(result, style = "points", error = "ci", title = "", x_label = "", y_label = "",
                      point_size = 2.6, symbol = 16, palette_name = "color", connect = TRUE,
                      varied_symbols = FALSE, legend = "bottom", x_angle = 0, show_stars = FALSE) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  d <- if (identical(result$config$scale, "log") && !is.null(result$raw_data)) result$raw_data else result$data
  method <- result$method; cfg <- result$config
  colors <- c("#007F83", "#E8914A", "#576CC0", "#C35B83", "#649D49", "#A477B2", "#A78744", "#738495")
  if (palette_name == "mono") colors <- c("#202020", "#666666", "#999999", "#BBBBBB")
  if (palette_name == "prism") colors <- c("#202020", "#2478B5", "#F28E2B", "#BD35DD", "#35A66D", "#D94855")
  palette <- grDevices::colorRampPalette(colors)
  if (method == "survival") {
    f <- result$fit
    labels <- if (is.null(f$strata)) rep(levels(d$group)[1], length(f$time)) else rep(sub("^group=", "", names(f$strata)), f$strata)
    curve <- data.frame(x = f$time, y = f$surv, group = labels, censored = f$n.censor)
    curve <- rbind(data.frame(x = 0, y = 1, group = unique(labels), censored = 0), curve)
    curve$group <- factor(curve$group, levels = levels(d$group))
    curve <- curve[order(curve$group, curve$x), ]
    p <- ggplot2::ggplot(curve, ggplot2::aes(x, y, color = group)) + ggplot2::geom_step(linewidth = .9) +
      ggplot2::geom_point(data = curve[curve$censored > 0, ], shape = 3, size = 2.8) +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    xl <- cfg$y; yl <- "Survival probability"
  } else if (method %in% c("chisq", "fisher")) {
    p <- ggplot2::ggplot(d, ggplot2::aes(group, count, fill = category)) + ggplot2::geom_col(position = "dodge", width = .7)
    xl <- cfg$group; yl <- "Count"
  } else if (method %in% c("pearson", "spearman", "linear", "dose")) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point(color = colors[1], shape = as.numeric(symbol), size = point_size, alpha = .85)
    if (method == "linear") p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, level = cfg$confidence, color = colors[1], fill = "#B7DADB", linewidth = .8)
    if (method == "dose") {
      grid <- data.frame(x = exp(seq(log(min(d$x)), log(max(d$x)), length.out = 250)))
      grid$y <- as.numeric(predict(result$fit, newdata = grid))
      p <- p + ggplot2::geom_line(data = grid, color = colors[1], linewidth = .9) + ggplot2::scale_x_log10()
    }
    xl <- cfg$x; yl <- cfg$y
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(group, y, color = group))
    if (style == "box") p <- p + ggplot2::geom_boxplot(width = .48, outlier.shape = NA, alpha = .13, ggplot2::aes(fill = group))
    if (style == "violin") p <- p + ggplot2::geom_violin(alpha = .12, width = .8, trim = FALSE, ggplot2::aes(fill = group))
    if ("id" %in% names(d) && isTRUE(connect)) p <- p + ggplot2::geom_line(ggplot2::aes(group = id), color = "#BECACC", alpha = .65)
    shape_mapping <- if (isTRUE(varied_symbols)) ggplot2::aes(shape = group) else NULL
    position <- ggplot2::position_jitter(width = if ("id" %in% names(d)) 0 else .09, height = 0, seed = 42)
    p <- p + if (isTRUE(varied_symbols)) ggplot2::geom_point(shape_mapping, position = position, size = point_size, alpha = .82) else ggplot2::geom_point(position = position, size = point_size, alpha = .82, shape = as.numeric(symbol))
    if (style == "points") {
      if ("factor2" %in% names(d)) {
        sums <- do.call(rbind, lapply(split(d, d$factor2), function(z) {
          s <- summary_table(z$y, z$group, cfg$confidence); s$factor2 <- as.character(z$factor2[1]); s
        }))
      } else sums <- if (method %in% t_test_methods) summary_table(d$y, d$group, cfg$confidence) else if (identical(cfg$scale, "log")) result$tables$original_summary else result$tables$summary
      sums$low <- if (error == "sd") sums$mean - sums$sd else if (error == "sem") sums$mean - sums$sem else sums$ci_low
      sums$high <- if (error == "sd") sums$mean + sums$sd else if (error == "sem") sums$mean + sums$sem else sums$ci_high
      if (error != "none") p <- p + ggplot2::geom_errorbar(data = sums, ggplot2::aes(x = group, ymin = low, ymax = high), inherit.aes = FALSE, width = .15, color = "#23353D", linewidth = .7)
      p <- p + ggplot2::geom_point(data = sums, ggplot2::aes(group, mean), inherit.aes = FALSE, shape = 95, size = 8, color = "#23353D")
    }
    if ("factor2" %in% names(d)) p <- p + ggplot2::facet_wrap(~factor2)
    xl <- if (is.null(cfg$group) || cfg$group == "") "Group" else cfg$group; yl <- cfg$y
  }
  n_colors <- max(8, if ("group" %in% names(d)) nlevels(d$group) else 1, if ("category" %in% names(d)) nlevels(d$category) else 1)
  if (isTRUE(varied_symbols) && "group" %in% names(d)) p <- p + ggplot2::scale_shape_manual(values = rep(c(16, 15, 17, 18, 1, 0, 2, 5), length.out = nlevels(d$group)))
  if (isTRUE(show_stars)) {
    tests <- result$tables$test
    if (is.null(tests) && !is.null(result$tables$anova)) {
      tests <- result$tables$anova
      tests$p_value <- tests[["Pr(>F)"]]
      tests <- tests[!is.na(tests$p_value), , drop = FALSE]
    }
    if (is.null(tests) && method == "linear") tests <- result$tables$coefficients[result$tables$coefficients$term != "(Intercept)", , drop = FALSE]
    if (!is.null(tests) && "p_value" %in% names(tests)) {
      if ("selected" %in% names(tests)) tests <- tests[tests$selected, , drop = FALSE]
      p <- p + ggplot2::labs(subtitle = paste("Test p =", format.pval(tests$p_value[1], digits = 3), significance_summary(tests$p_value[1])))
    }
  }
  p + ggplot2::scale_color_manual(values = palette(n_colors)) + ggplot2::scale_fill_manual(values = palette(n_colors)) +
    ggplot2::labs(title = if (nzchar(title)) title else names(analysis_methods)[match(method, analysis_methods)],
      x = if (nzchar(x_label)) x_label else xl, y = if (nzchar(y_label)) y_label else yl, color = NULL, fill = NULL,
      caption = if (method == "survival") "+ indicates censoring; confidence intervals and numbers at risk are in the results table." else if ("group" %in% names(d) && !method %in% c("chisq", "fisher") && style == "points") if (error == "none") "Points = observations; horizontal marks = means" else paste("Points = observations; bars =", switch(error, ci = paste0(cfg$confidence * 100, "% CI of mean"), sd = "SD", sem = "SEM")) else NULL) +
    ggplot2::theme_classic(base_size = 13) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 18)),
      axis.title = ggplot2::element_text(color = "#324650"), legend.position = legend, axis.text.x = ggplot2::element_text(angle = x_angle, hjust = if (x_angle == 0) .5 else 1), plot.caption = ggplot2::element_text(hjust = 0, color = "#64757C", size = 10),
      plot.margin = ggplot2::margin(18, 20, 18, 18))
}

diagnostic_plot <- function(result) {
  if (isTRUE(result$nested)) return(result$graphs$residual)
  d <- residual_frame(result)
  if (is.null(d)) return(NULL)
  ggplot2::ggplot(d, ggplot2::aes(fitted, residual)) + ggplot2::geom_hline(yintercept = 0, color = "#A0ADAF", linetype = 2) +
    ggplot2::geom_point(color = "#007F83", size = 2.5) + ggplot2::theme_classic(base_size = 13) + ggplot2::labs(x = "Fitted response", y = "Residual", title = "Residuals vs fitted values")
}

residual_graph <- function(result, kind) {
  # Use residuals from the fitted analysis, not differences recomputed from a plot.
  d <- residual_frame(result)
  assert(!is.null(d), "Residual graphs are unavailable for this analysis.")
  if (kind == "residual") return(diagnostic_plot(result))
  if (kind == "homoscedasticity") {
    p <- ggplot2::ggplot(d, ggplot2::aes(fitted, abs(residual), color = group)) + ggplot2::geom_point(size = 2.5) +
      ggplot2::labs(x = "Fitted response", y = "Absolute residual", title = "Homoscedasticity plot", color = NULL)
  } else if (kind == "residual_qq") {
    x <- sort(d$residual)
    assert(length(x) >= 2 && sd(x) > 0, "Residual QQ plots require varying residuals.")
    z <- data.frame(actual = x, predicted = mean(x) + sd(x)*qnorm(ppoints(length(x))))
    p <- ggplot2::ggplot(z, ggplot2::aes(actual, predicted)) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#69788a") +
      ggplot2::geom_point(color = "#2467ae", size = 2.5) + ggplot2::labs(x = "Actual residual", y = "Predicted normal residual", title = "Residual QQ plot")
  } else {
    subjects <- unique(d$subject)
    d$subject <- factor(d$subject, levels = rev(subjects))
    # Limit labels, not tiles, so large tables remain readable at app height.
    row_labels <- if (length(subjects) <= 10L) subjects else subjects[unique(round(seq(1, length(subjects), length.out = 10L)))]
    p <- ggplot2::ggplot(d, ggplot2::aes(group, subject, fill = residual)) + ggplot2::geom_tile(color = "white", linewidth = .2) +
      ggplot2::scale_y_discrete(breaks = row_labels) +
      ggplot2::scale_fill_gradient2(low = "#2467ae", mid = "white", high = "#b94040", midpoint = 0) +
      ggplot2::labs(x = "Condition", y = if ("id" %in% names(result$data)) "Subject" else "Observation within group", title = "Residual heatmap", fill = "Residual")
  }
  p + ggplot2::theme_classic(base_size = 13)
}

rank_graph <- function(result) {
  assert(result$method %in% c("kruskal", "friedman", "mann_whitney", "paired_wilcox"), "Rank graphs require a nonparametric analysis.")
  d <- result$data
  d$rank <- if ("id" %in% names(d)) ave(d$y, d$id, FUN = rank) else rank(d$y)
  ggplot2::ggplot(d, ggplot2::aes(group, rank, color = group)) + ggplot2::geom_point(position = ggplot2::position_jitter(width = .07, height = 0, seed = 42), alpha = .75, size = 2.5) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 8, color = "#25364a") + ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = "Condition", y = "Rank", color = NULL, title = "Ranks by condition", caption = if ("id" %in% names(d)) "Ranks calculated within each subject; horizontal marks show mean ranks." else "Ranks calculated across all selected observations; horizontal marks show mean ranks.")
}

difference_graph <- function(result) {
  d <- result$data; assert("id" %in% names(d), "Difference graphs require matched subject IDs.")
  cfg <- result$config
  if (cfg$comparisons %in% c("none", "trend")) cfg$comparisons <- "all"
  pairs <- comparison_pairs(levels(d$group), cfg)
  values <- lapply(pairs, function(p) {
    a <- d[d$group == p[1], ]; b <- d[d$group == p[2], ]; ids <- intersect(as.character(a$id), as.character(b$id))
    data.frame(subject = ids, difference = a$y[match(ids, a$id)] - b$y[match(ids, b$id)], comparison = rep(paste(p, collapse = " minus "), length(ids)))
  })
  z <- do.call(rbind, values); assert(!is.null(z) && nrow(z) > 0, "No overlapping subject measurements for the difference graph.")
  z$comparison <- factor(z$comparison, levels = unique(z$comparison))
  ggplot2::ggplot(z, ggplot2::aes(comparison, difference, color = comparison)) + ggplot2::geom_hline(yintercept = 0, color = "#69788a", linetype = 2) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .07, height = 0, seed = 42), size = 2.5, alpha = .8) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 8, color = "#25364a") + ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = if (cfg$scale == "log") "Within-subject log difference" else "Within-subject difference", title = "Matched differences", caption = "Each point is a subject with both measurements; horizontal marks show means.")
}

normality_plot <- function(result) {
  if (isTRUE(result$nested)) return(result$graphs$assumption_qq)
  d <- result$data; m <- result$method
  if (m %in% c("chisq", "fisher", "survival")) return(NULL)
  if (!is.null(result$normality_samples)) samples <- result$normality_samples
  else if (m %in% c("paired_t", "log_paired_t", "paired_wilcox")) {
    v <- split(d$y, d$group); samples <- data.frame(value = v[[1]] - v[[2]], sample = "Within-pair differences")
  } else if (m %in% c("linear", "dose")) samples <- data.frame(value = residuals(result$fit), sample = "Model residuals")
  else if (m %in% c("pearson", "spearman")) samples <- rbind(data.frame(value = d$x, sample = "X"), data.frame(value = d$y, sample = "Y"))
  else if (m %in% c("rm_anova", "mixed", "friedman", "log_rm_anova", "log_mixed")) {
    subjects <- split(d, d$id)
    samples <- do.call(rbind, lapply(combn(levels(d$group), 2, simplify = FALSE), function(p) {
      v <- unlist(lapply(subjects, function(s) {
        a <- s$y[s$group == p[1]]; b <- s$y[s$group == p[2]]
        if (length(a) == 1 && length(b) == 1) a - b else NULL
      }))
      data.frame(value = v, sample = rep(paste(p, collapse = " minus "), length(v)))
    }))
  } else samples <- data.frame(value = d$y, sample = d$group)
  if (is.null(samples) || !nrow(samples)) return(NULL)
  if (!is.factor(samples$sample)) samples$sample <- factor(samples$sample, levels = unique(samples$sample))
  ggplot2::ggplot(samples, ggplot2::aes(sample = value)) + ggplot2::stat_qq(color = "#007F83") +
    ggplot2::stat_qq_line(color = "#829399") + ggplot2::facet_wrap(~sample, scales = "free") +
    ggplot2::theme_classic(base_size = 12) + ggplot2::labs(x = "Normal theoretical quantile", y = "Observed quantile", title = "Normal QQ plots")
}

distribution_qq_plot <- function(result, lognormal = FALSE, title = "", palette_name = "prism", legend = "right", point_size = 2.6, symbol = 16, varied_symbols = TRUE, x_angle = 0, x_label = "", y_label = "", ...) {
  # Column-distribution QQ plots use observed measurements. Lognormal plots
  # require positive values and are distinct from model-residual diagnostics.
  assert(!result$method %in% c("survival", "chisq", "fisher"), "Column distribution QQ plots require uncensored numeric responses, rather than event-time or count-table analyses.")
  d <- result$column_data
  if (is.null(d)) {
    d <- result$raw_data
    if (is.null(d) || !all(c("group", "y") %in% names(d))) {
      d <- result$data; d <- data.frame(y = d$y, group = factor("Y"))
    }
  }
  skipped <- character()
  points <- lapply(split(d$y, d$group, drop = TRUE), function(x) {
    x <- sort(x[is.finite(x)])
    if (length(x) < 2 || length(unique(x)) < 2 || (lognormal && any(x <= 0))) return(NULL)
    model_values <- if (lognormal) log(x) else x
    prediction <- mean(model_values) + sd(model_values)*qnorm(ppoints(length(x)))
    data.frame(actual = x, predicted = if (lognormal) exp(prediction) else prediction)
  })
  valid <- !vapply(points, is.null, logical(1)); skipped <- names(points)[!valid]
  assert(any(valid), if (lognormal) "No varying, strictly positive column is available for the lognormal QQ plot." else "No varying column has enough values for a QQ plot.")
  p <- do.call(rbind, lapply(names(points)[valid], function(g) { z <- points[[g]]; z$group <- g; z }))
  p$group <- factor(p$group, levels = names(points)[valid])
  colors <- if (palette_name == "mono") c("#202020", "#999999") else if (palette_name == "prism") c("#202020", "#2478B5", "#F28E2B", "#BD35DD") else c("#007F83", "#E8914A", "#576CC0", "#C35B83")
  plot <- ggplot2::ggplot(p, ggplot2::aes(actual, predicted, color = group)) + ggplot2::geom_abline(slope = 1, intercept = 0, color = "#D45A5A", linetype = 3)
  plot <- plot + if (varied_symbols) ggplot2::geom_point(ggplot2::aes(shape = group), size = point_size, alpha = .7) else ggplot2::geom_point(shape = as.numeric(symbol), size = point_size, alpha = .7)
  plot + ggplot2::scale_color_manual(values = grDevices::colorRampPalette(colors)(length(unique(p$group)))) +
    ggplot2::scale_shape_manual(values = rep(c(16, 15, 17, 18, 1, 0, 2, 5), length.out = length(unique(p$group)))) +
    ggplot2::theme_classic(base_size = 13) + ggplot2::theme(legend.position = legend, axis.text.x = ggplot2::element_text(angle = x_angle, hjust = if (x_angle == 0) .5 else 1)) +
    ggplot2::labs(title = if (nzchar(title)) title else if (lognormal) "Lognormal QQ plot" else "Normal QQ plot", x = if (nzchar(x_label)) x_label else "Actual", y = if (nzchar(y_label)) y_label else "Predicted", color = NULL, shape = NULL,
      caption = if (length(skipped)) paste("Unavailable columns:", paste(skipped, collapse = ", ")) else NULL)
}

comparison_plot <- function(result, title = "", x_label = "", y_label = "", point_size = 2.6, symbol = 16, ...) {
  table <- result$tables$posthoc
  assert(!is.null(table), "This analysis has no multiple-comparison table.")
  if (all(c("diff", "lwr", "upr") %in% names(table))) {
    table$estimate <- table$diff; table$ci_low <- table$lwr; table$ci_high <- table$upr
  }
  assert(all(c("estimate", "ci_low", "ci_high") %in% names(table)), "This comparison procedure does not supply confidence intervals.")
  ratio <- "geometric_mean_ratio" %in% names(table)
  if (ratio) { table$estimate <- table$geometric_mean_ratio; table$ci_low <- table$ratio_ci_low; table$ci_high <- table$ratio_ci_high }
  table <- table[is.finite(table$estimate) & !is.na(table$ci_low) & !is.na(table$ci_high), , drop = FALSE]
  assert(nrow(table) > 0, "There are no estimable comparison intervals to plot.")
  if (result$method %in% multiple_t_methods && identical(result$config$t_test_mode, "endpoints")) table$contrast <- paste(table$test, table$contrast, sep = ": ")
  table$open_low <- !is.finite(table$ci_low) | (ratio & table$ci_low <= 0)
  table$open_high <- !is.finite(table$ci_high)
  # Draw arrows for unbounded one-sided limits rather than discarding the tests.
  values <- c(table$estimate, table$ci_low[!table$open_low], table$ci_high[!table$open_high])
  if (ratio) values <- log10(values[values > 0])
  extent <- range(values[is.finite(values)]); pad <- max(diff(extent) * .15, .1)
  low <- extent[1] - pad; high <- extent[2] + pad
  table$ci_low[table$open_low] <- if (ratio) 10^low else low
  table$ci_high[table$open_high] <- if (ratio) 10^high else high
  if (!"contrast" %in% names(table)) table$contrast <- paste(table$group_1, table$group_2, sep = " minus ")
  table$contrast <- factor(table$contrast, levels = rev(unique(table$contrast)))
  label <- if ("interval" %in% names(table)) unique(table$interval)[1] else paste(result$config$confidence*100, "% simultaneous Tukey")
  reference <- if (result$method %in% multiple_t_methods) if (ratio) result$config$null_ratio else result$config$null_value else if (ratio) 1 else 0
  p <- ggplot2::ggplot(table, ggplot2::aes(estimate, contrast)) + ggplot2::geom_vline(xintercept = reference, linetype = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = ci_low, xend = ci_high, yend = contrast)) +
    ggplot2::geom_errorbar(data = table[!table$open_low & !table$open_high, , drop = FALSE], ggplot2::aes(xmin = ci_low, xmax = ci_high), width = .15, orientation = "y") +
    ggplot2::geom_segment(data = table[table$open_low, , drop = FALSE], ggplot2::aes(xend = ci_low, yend = contrast), arrow = grid::arrow(length = grid::unit(.08, "inches"))) +
    ggplot2::geom_segment(data = table[table$open_high, , drop = FALSE], ggplot2::aes(xend = ci_high, yend = contrast), arrow = grid::arrow(length = grid::unit(.08, "inches"))) +
    ggplot2::geom_point(size = point_size, shape = as.numeric(symbol)) +
    ggplot2::theme_classic(base_size = 13) + ggplot2::labs(title = if (nzchar(title)) title else paste(label, "confidence intervals"),
      caption = if (any(table$open_low | table$open_high)) "Arrows mark unbounded confidence limits; the displayed endpoints are not finite bounds." else NULL,
      x = if (nzchar(x_label)) x_label else if (ratio) "Ratio of geometric means" else "Difference between group means", y = if (nzchar(y_label)) y_label else NULL)
  if (ratio) p <- p + ggplot2::scale_x_log10()
  p
}

graph_kinds <- c("Data graph" = "raw", "Normal QQ (columns)" = "normal_qq", "Lognormal QQ (columns)" = "lognormal_qq",
  "Publication SuperPlot (linear mixed analysis)" = "publication", "Residual QQ by condition (linear mixed analysis)" = "condition_residual_qq",
  "Multiple-comparison intervals" = "comparison", "Residual diagnostic" = "residual", "Assumption QQ (differences / residuals)" = "assumption_qq",
  "Homoscedasticity" = "homoscedasticity", "Residual QQ" = "residual_qq", "Residual heatmap" = "residual_heatmap", "Ranks" = "ranks", "Matched differences" = "differences")

make_graph <- function(result, kind = "raw", options = list()) {
  assert(kind %in% graph_kinds, "Unknown graph type.")
  if (isTRUE(result$nested) && !kind %in% c("normal_qq", "lognormal_qq")) {
    plot <- result$graphs[[kind]]
    assert(!is.null(plot), "Enable this linear mixed analysis figure in Analysis options before running the analysis.")
    return(plot)
  }
  assert(!kind %in% c("publication", "condition_residual_qq"), "This graph requires a linear mixed analysis.")
  if (kind == "raw") return(do.call(make_plot, c(list(result = result), options)))
  if (kind == "comparison") return(do.call(comparison_plot, c(list(result = result), options)))
  if (kind %in% c("normal_qq", "lognormal_qq")) return(do.call(distribution_qq_plot, c(list(result = result, lognormal = kind == "lognormal_qq"), options)))
  plot <- if (kind %in% c("residual", "homoscedasticity", "residual_qq", "residual_heatmap")) residual_graph(result, kind) else if (kind == "ranks") rank_graph(result) else if (kind == "differences") difference_graph(result) else normality_plot(result)
  assert(!is.null(plot), "This analysis has no graph of the selected type.")
  if (!is.null(options$title) && nzchar(options$title)) plot <- plot + ggplot2::labs(title = options$title)
  plot
}

draw_graph_layout <- function(plots, columns = 2L, title = "") {
  # Draw ggplots and publication grid grobs through their respective APIs so
  # fixed physical dimensions in nested publication panels remain meaningful.
  assert(length(plots) > 0, "Select at least one saved graph for this layout.")
  columns <- max(1L, min(6L, as.integer(columns))); rows <- ceiling(length(plots)/columns)
  grid::grid.newpage()
  heights <- if (nzchar(title)) grid::unit.c(grid::unit(.4, "inches"), grid::unit(rep(1, rows), "null")) else grid::unit(rep(1, rows), "null")
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(rows + as.integer(nzchar(title)), columns, heights = heights)))
  on.exit(grid::popViewport())
  if (nzchar(title)) grid::grid.text(title, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(columns)), gp = grid::gpar(fontface = "bold", fontsize = 16))
  offset <- as.integer(nzchar(title))
  for (i in seq_along(plots)) {
    vp <- grid::viewport(layout.pos.row = ceiling(i/columns) + offset, layout.pos.col = (i - 1) %% columns + 1)
    if (inherits(plots[[i]], "grob")) {
      grid::pushViewport(vp); grid::grid.draw(plots[[i]]); grid::popViewport()
    } else print(plots[[i]], newpage = FALSE, vp = vp)
  }
  invisible(NULL)
}

overlay_graphs <- function(entries, title = "") {
  assert(length(entries) >= 2, "Select at least two data graphs for an overlay.")
  assert(!any(vapply(entries, function(e) isTRUE(e$result$nested), logical(1))), "Use separate panels for linear mixed analysis SuperPlots so model means and mouse / slice identities are preserved.")
  assert(all(vapply(entries, function(e) e$kind == "raw", logical(1))), "Overlays combine data graphs; use a panel layout for QQ, residual, or comparison graphs.")
  xy <- vapply(entries, function(e) e$result$method %in% c("linear", "pearson", "spearman"), logical(1))
  column <- vapply(entries, function(e) all(c("group", "y") %in% names(e$result$data)) && !e$result$method %in% c("survival", "chisq", "fisher", "two_way"), logical(1))
  assert(all(xy) || all(column), "Overlay graphs need compatible XY or Column axes. Use a panel layout for different analysis types.")
  parts <- lapply(entries, function(e) {
    d <- if (identical(e$result$config$scale, "log")) e$result$raw_data else e$result$data
    d$series <- e$name
    if (all(xy)) d[, c("x", "y", "series")] else {
      d$subject <- if ("id" %in% names(d) && !identical(e$options$connect, FALSE)) paste(e$name, d$id, sep = ":") else NA_character_
      d[, c("group", "y", "series", "subject")]
    }
  })
  d <- do.call(rbind, parts)
  if (all(xy)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y, color = series, shape = series)) + ggplot2::geom_point(size = 2.4, alpha = .7)
    if (all(vapply(entries, function(e) e$result$method == "linear", logical(1)))) p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE)
    xl <- entries[[1]]$result$config$x
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(group, y, color = series, shape = series)) + ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .08, dodge.width = .45, seed = 42), size = 2.3, alpha = .65)
    if (any(!is.na(d$subject))) p <- p + ggplot2::geom_line(data = d[!is.na(d$subject), ], ggplot2::aes(group = subject), alpha = .3, linewidth = .4, position = ggplot2::position_dodge(width = .45))
    sums <- do.call(rbind, lapply(split(d, d$series), function(z) { s <- summary_table(z$y, z$group); s$series <- z$series[1]; s }))
    p <- p + ggplot2::geom_errorbar(data = sums, ggplot2::aes(group, ymin = ci_low, ymax = ci_high, color = series, group = series), inherit.aes = FALSE, position = ggplot2::position_dodge(.45), width = .12) +
      ggplot2::geom_point(data = sums, ggplot2::aes(group, mean, color = series, group = series), inherit.aes = FALSE, position = ggplot2::position_dodge(.45), shape = 95, size = 7)
    xl <- "Condition"
  }
  p + ggplot2::theme_classic(base_size = 13) + ggplot2::theme(legend.position = "bottom") + ggplot2::labs(title = title, x = xl, y = entries[[1]]$result$config$y, color = NULL, shape = NULL)
}
