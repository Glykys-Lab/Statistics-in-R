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
      } else sums <- if (identical(cfg$scale, "log")) result$tables$original_summary else result$tables$summary
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
  if (is.null(result$fit) || !result$method %in% c("anova", "two_way", "linear", "dose", "rm_anova", "mixed", "log_rm_anova", "log_mixed")) return(NULL)
  d <- data.frame(fitted = fitted(result$fit), residual = residuals(result$fit))
  ggplot2::ggplot(d, ggplot2::aes(fitted, residual)) + ggplot2::geom_hline(yintercept = 0, color = "#A0ADAF", linetype = 2) +
    ggplot2::geom_point(color = "#007F83", size = 2.5) + ggplot2::theme_classic(base_size = 13) + ggplot2::labs(x = "Fitted response", y = "Residual", title = "Residuals vs fitted values")
}

normality_plot <- function(result) {
  d <- result$data; m <- result$method
  if (m %in% c("chisq", "fisher", "survival")) return(NULL)
  if (m %in% c("paired_t", "paired_wilcox")) {
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
  ggplot2::ggplot(samples, ggplot2::aes(sample = value)) + ggplot2::stat_qq(color = "#007F83") +
    ggplot2::stat_qq_line(color = "#829399") + ggplot2::facet_wrap(~sample, scales = "free") +
    ggplot2::theme_classic(base_size = 12) + ggplot2::labs(x = "Normal theoretical quantile", y = "Observed quantile", title = "Normal QQ plots")
}

distribution_qq_plot <- function(result, lognormal = FALSE, title = "", palette_name = "prism", legend = "right", point_size = 2.6, symbol = 16, varied_symbols = TRUE, x_angle = 0, x_label = "", y_label = "", ...) {
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
  table <- table[is.finite(table$estimate) & is.finite(table$ci_low) & is.finite(table$ci_high), , drop = FALSE]
  assert(nrow(table) > 0, "There are no estimable comparison intervals to plot.")
  if (!"contrast" %in% names(table)) table$contrast <- paste(table$group_1, table$group_2, sep = " minus ")
  table$contrast <- factor(table$contrast, levels = rev(unique(table$contrast)))
  label <- if ("interval" %in% names(table)) unique(table$interval)[1] else paste(result$config$confidence*100, "% simultaneous Tukey")
  p <- ggplot2::ggplot(table, ggplot2::aes(estimate, contrast)) + ggplot2::geom_vline(xintercept = if (ratio) 1 else 0, linetype = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = .15, orientation = "y") + ggplot2::geom_point(size = point_size, shape = as.numeric(symbol)) +
    ggplot2::theme_classic(base_size = 13) + ggplot2::labs(title = if (nzchar(title)) title else paste(label, "confidence intervals"),
      x = if (nzchar(x_label)) x_label else if (ratio) "Ratio of geometric means" else "Difference between group means", y = if (nzchar(y_label)) y_label else NULL)
  if (ratio) p <- p + ggplot2::scale_x_log10()
  p
}

graph_kinds <- c("Data graph" = "raw", "Normal QQ (columns)" = "normal_qq", "Lognormal QQ (columns)" = "lognormal_qq",
  "Multiple-comparison intervals" = "comparison", "Residual diagnostic" = "residual", "Assumption QQ (differences / residuals)" = "assumption_qq")

make_graph <- function(result, kind = "raw", options = list()) {
  assert(kind %in% graph_kinds, "Unknown graph type.")
  if (kind == "raw") return(do.call(make_plot, c(list(result = result), options)))
  if (kind == "comparison") return(do.call(comparison_plot, c(list(result = result), options)))
  if (kind %in% c("normal_qq", "lognormal_qq")) return(do.call(distribution_qq_plot, c(list(result = result, lognormal = kind == "lognormal_qq"), options)))
  plot <- if (kind == "residual") diagnostic_plot(result) else normality_plot(result)
  assert(!is.null(plot), "This analysis has no graph of the selected type.")
  if (!is.null(options$title) && nzchar(options$title)) plot <- plot + ggplot2::labs(title = options$title)
  plot
}

draw_graph_layout <- function(plots, columns = 2L, title = "") {
  assert(length(plots) > 0, "Select at least one saved graph for this layout.")
  columns <- max(1L, min(6L, as.integer(columns))); rows <- ceiling(length(plots)/columns)
  grid::grid.newpage()
  heights <- if (nzchar(title)) grid::unit.c(grid::unit(.4, "inches"), grid::unit(rep(1, rows), "null")) else grid::unit(rep(1, rows), "null")
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(rows + as.integer(nzchar(title)), columns, heights = heights)))
  on.exit(grid::popViewport())
  if (nzchar(title)) grid::grid.text(title, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(columns)), gp = grid::gpar(fontface = "bold", fontsize = 16))
  offset <- as.integer(nzchar(title))
  for (i in seq_along(plots)) print(plots[[i]], newpage = FALSE, vp = grid::viewport(layout.pos.row = ceiling(i/columns) + offset, layout.pos.col = (i - 1) %% columns + 1))
  invisible(NULL)
}

overlay_graphs <- function(entries, title = "") {
  assert(length(entries) >= 2, "Select at least two data graphs for an overlay.")
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
