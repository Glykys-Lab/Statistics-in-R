# Settings and browser preferences; no experimental data is stored in preferences.
one_way_methods <- c("anova", "welch_anova", "kruskal", "friedman", "rm_anova", "mixed", "rm_auto", "log_rm_anova", "log_mixed", "log_rm_auto")
repeated_methods <- c("rm_anova", "mixed", "rm_auto", "log_rm_anova", "log_mixed", "log_rm_auto")
input_or <- function(id, fallback) { value <- input[[id]]; if (is.null(value)) fallback else value }
pending_preferences <- reactiveVal(list())
advanced_config <- reactive({
  m <- input$method; req(m); z <- list()
  if (m %in% one_way_methods) {
    z <- list(reverse_comparisons = isTRUE(input$reverse_comparisons), report_adjusted_p = input_or("report_adjusted_p", TRUE), report_descriptive = input_or("report_descriptive", TRUE))
    if (identical(input$comparisons, "selected")) z$selected_pairs <- as.character(input_or("selected_pairs", character()))
    if (input_or("adjustment", "sidak") %in% c("BKY", "BH", "BY") || identical(input$error_control, "fdr")) z$fdr_q <- input_or("fdr_q", .05)
    if (m %in% c("kruskal", "friedman")) z$rank_test <- input_or("rank_test", "wilcox")
    if (m %in% c("anova", repeated_methods)) for (id in c("report_goodness", "report_effects", "report_aicc")) z[[id]] <- isTRUE(input[[id]])
    if (m %in% c("mixed", "rm_auto", "log_mixed", "log_rm_auto")) { z$random_effect <- input_or("random_effect", "keep"); z$glm_initial <- input_or("glm_initial", TRUE) }
    z$residual_normality <- isTRUE(input$residual_normality)
    if (!m %in% c(repeated_methods, "friedman")) z$residual_variance <- isTRUE(input$residual_variance)
    graphs <- c(residual = "graph_residual", homoscedasticity = "graph_homoscedasticity", residual_qq = "graph_residual_qq", residual_heatmap = "graph_residual_heatmap")
    z$diagnostic_graphs <- names(graphs)[vapply(graphs, function(id) isTRUE(input[[id]]), logical(1))]
    if (!input_or("comparisons", "all") %in% c("none", "trend") && !m %in% c("kruskal", "friedman") && isTRUE(input_or("graph_comparison", TRUE))) z$diagnostic_graphs <- c(z$diagnostic_graphs, "comparison")
    if (m %in% c("kruskal", "friedman") && isTRUE(input$graph_ranks)) z$diagnostic_graphs <- c(z$diagnostic_graphs, "ranks")
    if (m %in% c(repeated_methods, "friedman") && isTRUE(input$graph_differences)) z$diagnostic_graphs <- c(z$diagnostic_graphs, "differences")
  }
  z
})
output$comparison_family_ui <- renderUI({
  choices <- c("None: no follow-up comparisons" = "none", "Every pair of columns" = "all", "Each column versus a control" = "control", "Selected pairs of columns" = "selected")
  can_trend <- identical(input$method, "anova") || (input$method %in% repeated_methods && identical(input$sphericity, "none"))
  if (can_trend) choices <- c(choices, "Linear trend across column positions" = "trend")
  previous <- pending_preferences()$comparisons
  if (is.null(previous)) previous <- isolate(input$comparisons)
  if (is.null(previous) || !previous %in% choices) previous <- "all"
  tagList(radioButtons("comparisons", "Follow-up tests", choices, selected = previous), if (!can_trend) p(class = "hint", "Linear trend requires ordinary ANOVA, or repeated parametric analysis with sphericity assumed."))
})
output$selected_pairs_ui <- renderUI({
  d <- data(); req(input$group %in% names(d)); groups <- levels(group_column(d, input$group)); req(length(groups) >= 2)
  pairs <- combn(groups, 2, simplify = FALSE)
  values <- vapply(pairs, function(p) paste(vapply(p, utils::URLencode, character(1), reserved = TRUE), collapse = "::"), character(1))
  selectInput("selected_pairs", "Pairs to compare", setNames(values, vapply(pairs, paste, character(1), collapse = " versus ")), selected = intersect(isolate(input$selected_pairs), values), multiple = TRUE)
})
output$analysis_plan <- renderText({ m <- input$method; req(m); paste(names(analysis_methods)[match(m, analysis_methods)], "\n", method_note(m), "\nFollow-up comparisons:", input_or("comparisons", "all"), "\nCorrection:", input_or("adjustment", "sidak")) })
observeEvent(input$show_options, updateTabsetPanel(session, "tabs", selected = "options"))
observeEvent(input$alpha, { req(is.finite(input$alpha), input$alpha > 0, input$alpha < 1); value <- 1 - input$alpha; choices <- sort(unique(c(.90, .95, .99, value))); updateSelectInput(session, "confidence", choices = setNames(choices, paste0(format(100*choices, trim = TRUE), "%")), selected = value) }, ignoreInit = TRUE)
observeEvent(input$confidence, { value <- suppressWarnings(as.numeric(input$confidence)); req(length(value) == 1, is.finite(value), value > 0, value < 1); updateNumericInput(session, "alpha", value = 1 - value) }, ignoreInit = TRUE)
display_options <- reactive(list(output_digits = as.integer(input_or("output_digits", 4L)), p_style = input_or("p_style", "prism"), p_digits = as.integer(input_or("p_digits", 6L))))
display_tables <- reactive(visible_result_tables(result()))
preference_choices <- list(matching = c("independent", "matched"), distribution = c("normal", "lognormal", "rank"), variance = c("unequal", "equal"), repeated = c("complete", "mixed", "auto"), sphericity = c("gg", "none"), random_effect = c("remove", "full", "keep"), comparisons = c("none", "all", "control", "selected", "trend"), error_control = c("fwer", "fdr", "none"), adjustment = c("sidak", "holm_sidak", "bonferroni", "holm", "tukey", "BKY", "BH", "BY", "none"), rank_test = c("dunn", "wilcox"), p_style = c("prism", "decimal", "scientific", "stars"))
preference_numbers <- list(alpha = c(.001, .2), confidence = c(.8, .999), fdr_q = c(.001, .5), output_digits = c(2, 12), p_digits = c(1, 15))
preference_checks <- c("glm_initial", "reverse_comparisons", "report_adjusted_p", "report_descriptive", "report_goodness", "report_effects", "report_aicc", "graph_comparison", "graph_ranks", "graph_differences", "graph_residual", "graph_homoscedasticity", "graph_residual_qq", "graph_residual_heatmap", "residual_normality", "residual_variance")
preferences <- reactiveVal(list())
restore_preferences <- function(values) {
  req(is.list(values), length(values) <= 50); accepted <- list()
  for (id in intersect(names(values), c(names(preference_choices), names(preference_numbers), preference_checks))) {
    value <- values[[id]]
    if (length(value) != 1 || is.list(value) || is.na(value)) next
    if (id %in% names(preference_choices) && !value %in% preference_choices[[id]]) next
    if (id %in% names(preference_numbers) && (!is.numeric(value) || !is.finite(value) || value < preference_numbers[[id]][1] || value > preference_numbers[[id]][2])) next
    if (id %in% preference_checks && !is.logical(value)) next
    accepted[[id]] <- value
    if (id %in% preference_checks) updateCheckboxInput(session, id, value = value)
    else if (id %in% c("matching", "distribution", "repeated", "random_effect", "comparisons", "error_control")) updateRadioButtons(session, id, selected = value)
    else if (id %in% setdiff(names(preference_numbers), "confidence")) updateNumericInput(session, id, value = value)
    else updateSelectInput(session, id, selected = value)
  }
  preferences(accepted)
  pending_preferences(accepted[intersect(names(accepted), c("distribution", "comparisons", "adjustment"))])
}
observeEvent(list(input$distribution, input$comparisons, input$adjustment), {
  pending <- isolate(pending_preferences())
  for (id in names(pending)) if (identical(input[[id]], pending[[id]])) pending[[id]] <- NULL
  pending_preferences(pending)
})
observeEvent(input$restore_analysis_defaults, restore_preferences(input$restore_analysis_defaults))
observeEvent(input$defaults_storage_error, showNotification("This browser could not save analysis defaults. Check its storage settings.", type = "error"))
observeEvent(input$save_analysis_defaults, {
  ids <- c(names(preference_choices), names(preference_numbers), preference_checks); values <- setNames(lapply(ids, function(id) input[[id]]), ids); values <- values[!vapply(values, is.null, logical(1))]
  if (!is.null(values$confidence)) values$confidence <- as.numeric(values$confidence)
  session$sendCustomMessage("analysis-defaults", list(values = values)); showNotification("Analysis defaults saved in this browser.")
})
observeEvent(input$reset_analysis_defaults, {
  values <- list(matching = "independent", distribution = "normal", variance = "unequal", repeated = "auto", sphericity = "gg", random_effect = "remove", glm_initial = TRUE, comparisons = "all", error_control = "fwer", adjustment = "sidak", alpha = .05, confidence = .95, fdr_q = .05, reverse_comparisons = FALSE, report_adjusted_p = TRUE, report_descriptive = TRUE, rank_test = "dunn", output_digits = 4, p_style = "prism", p_digits = 6)
  for (id in setdiff(preference_checks, names(values))) values[[id]] <- identical(id, "graph_comparison")
  restore_preferences(values); session$sendCustomMessage("analysis-defaults", list(clear = TRUE)); showNotification("Analysis defaults reset.")
})
