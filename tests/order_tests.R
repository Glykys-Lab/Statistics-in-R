# Order regressions sourced by run_tests.R using its check/near helpers.
# Deliberately reverse selector order and leave a first measurement blank to verify
# that tables, contrasts, graph axes and legends still follow the worksheet sequence.
# Exercise actual ggplot scales with the experimental sequence from the worksheet.
steps <- c("Baseline2 - Baseline", "aCSF1 - Baseline", "aCSF2 - Baseline")
worksheet <- data.frame(Subject = paste0("S", 1:8), check.names = FALSE)
worksheet[[steps[1]]] <- c(NA, 2, 4, 3, 6, 5, 9, 7)
worksheet[[steps[2]]] <- c(3, 5, 4, 8, 7, 6, 10, 9)
worksheet[[steps[3]]] <- c(5, 8, 7, 9, 12, 10, 14, 11)
observations <- wide_to_long(worksheet, rev(steps), "Subject")
check(identical(levels(observations$group), steps), "Worksheet order overrides measurement selector click order")
order_result <- run_analysis(observations, list(method = "welch_anova", y = "response", group = "group"))
check(identical(order_result$tables$summary$group, steps), "Missing first measurement preserves the experimental group sequence")
for (style in c("points", "box", "violin")) {
  built <- ggplot2::ggplot_build(make_plot(order_result, style = style, varied_symbols = TRUE))
  check(identical(built$layout$panel_params[[1]]$x$get_labels(), steps) &&
    identical(as.character(built$plot$scales$get_scales("colour")$get_breaks()), steps) &&
    identical(as.character(built$plot$scales$get_scales("shape")$get_breaks()), steps),
    paste(style, "graph axis and legends follow worksheet order"))
}
for (kind in c("normal_qq", "lognormal_qq")) {
  built <- ggplot2::ggplot_build(make_graph(order_result, kind))
  check(identical(as.character(built$plot$scales$get_scales("colour")$get_breaks()), steps), paste(kind, "legend follows worksheet order"))
}
long_order <- observations
long_order$group <- as.character(long_order$group)
check(identical(levels(group_column(long_order, "group")), steps), "Long-format categories follow first appearance")
shuffled_order <- observations[c(9:24, 1:8), ]
check(identical(levels(group_column(shuffled_order, "group")), steps), "Explicit factor sequence survives row rearrangement")
before_after <- wide_to_long(column_demo("paired"), c("Before", "After"), "Subject")
paired_order <- run_analysis(before_after, list(method = "paired_t", y = "response", group = "group", id = "subject"))
check(paired_order$tables$test$contrast == "Before minus After" &&
  near(paired_order$tables$test$estimate, mean(column_demo("paired")$Before - column_demo("paired")$After)),
  "Paired contrast follows table order and preserves subject matching")
