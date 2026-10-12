# Build the ordinary one-way/repeated-measures settings tabs and correction menus.
# Values are translated into configuration by options_server.inc; keep input IDs
# stable because saved browser defaults refer to them. Nested models use nested_ui.R.
comparison_correction_choices <- function(method = NULL, comparisons = "all", approach = "fwer") {
  if (identical(approach, "fdr")) return(c("Two-stage Benjamini-Krieger-Yekutieli" = "BKY", "Benjamini-Hochberg" = "BH", "Benjamini-Yekutieli" = "BY"))
  if (identical(approach, "none")) return(c("No correction" = "none"))
  choices <- c("Sidak (adjusted p-values and intervals)" = "sidak", "Holm-Sidak (p-values; unadjusted intervals)" = "holm_sidak", "Holm (p-values; unadjusted intervals)" = "holm", "Bonferroni (p-values and intervals)" = "bonferroni")
  if (identical(method, "anova") && identical(comparisons, "all")) choices <- c(choices, "Tukey (all pairs only)" = "tukey")
  choices
}

analysis_options_ui <- function() tagList(
  h3("One-way analysis options"),
  tabsetPanel(id = "analysis_option_tabs",
    tabPanel("Experimental Design", value = "design", div(class = "options-tab-body",
      radioButtons("matching", "Experimental design", c("No matching / independent observations" = "independent", "Each row represents matched / repeated measurements" = "matched")),
      uiOutput("distribution_ui"),
      conditionalPanel("input.matching == 'independent' && input.distribution == 'normal'", selectInput("variance", "Variances", c("Unequal / Welch" = "unequal", "Assume equal" = "equal"))),
      conditionalPanel("input.matching == 'matched' && input.distribution != 'rank'", selectInput("sphericity", "Assume sphericity (equal variability of differences)?", c("No: Geisser-Greenhouse correction" = "gg", "Yes: no correction" = "none"))),
      div(class = "method-note", textOutput("analysis_plan")))),
    tabPanel("Repeated Measures", value = "repeated", div(class = "options-tab-body",
      conditionalPanel("input.matching == 'matched' && input.distribution != 'rank'",
        radioButtons("repeated", "Analyze using which method?", c("Repeated-measures ANOVA: complete subjects only" = "complete", "Mixed-effects model (REML): retain available measurements" = "mixed", "It depends: ANOVA for complete data, mixed effects when values are missing" = "auto"), selected = "auto"),
        radioButtons("random_effect", "When the subject component is zero or negative", c("Remove the subject effect and refit the simpler model" = "remove", "Fit the full model, allowing a signed covariance component" = "full", "Keep a nonnegative random intercept (legacy behavior)" = "keep")),
        checkboxInput("glm_initial", "Use GLM-based initial values as well as a small starting value", TRUE),
        p(class = "hint", "The full-model option permits negative within-subject covariance while requiring a positive-definite observation covariance. A negative covariance component is not a physical variance.")),
      conditionalPanel("input.matching != 'matched' || input.distribution == 'rank'", p("These fitting choices apply to Gaussian or lognormal matched data.")))),
    tabPanel("Multiple Comparisons", value = "comparisons", div(class = "options-tab-body",
      uiOutput("comparison_family_ui"),
      conditionalPanel("input.comparisons == 'control'", uiOutput("control_ui")),
      conditionalPanel("input.comparisons == 'selected'", uiOutput("selected_pairs_ui")),
      conditionalPanel("input.comparisons == 'trend'", p(class = "hint", "Tests trend across equally spaced column positions in worksheet order. Requires sphericity assumed for repeated data.")),
      conditionalPanel("['kruskal','friedman'].includes(input.method)", selectInput("rank_test", "Rank comparison procedure", c("Dunn comparisons with tie correction" = "dunn", "Pairwise Wilcoxon (normal approximation)" = "wilcox"))),
      p(class = "hint", "Choose the correction and reporting options on the Options tab."))),
    tabPanel("Options", value = "options", div(class = "options-tab-body",
      conditionalPanel("input.comparisons != 'none' && input.comparisons != 'trend'",
        radioButtons("error_control", "Multiple-comparison approach", c("Control the family-wise error rate" = "fwer", "Control the false discovery rate (FDR)" = "fdr", "Do not correct: each comparison stands alone" = "none")),
        selectInput("adjustment", "Multiple-comparison correction", comparison_correction_choices(), selected = "sidak"),
        conditionalPanel("input.error_control == 'fdr'", numericInput("fdr_q", "Target false discovery rate Q (0.05 = 5%)", .05, min = .001, max = .5, step = .01)),
        checkboxInput("reverse_comparisons", "Swap comparison direction (A minus B / B minus A)", FALSE),
        checkboxInput("report_adjusted_p", "Report adjusted p-values for each comparison", TRUE)),
      fluidRow(column(6, numericInput("alpha", "Family-wise significance threshold (alpha)", .05, min = .001, max = .2, step = .01)),
        column(6, selectInput("confidence", "Confidence level", c("95%" = .95, "90%" = .90, "99%" = .99)))),
      h4("Graphing"),
      conditionalPanel("input.comparisons != 'none' && input.comparisons != 'trend' && !['kruskal','friedman'].includes(input.method)", checkboxInput("graph_comparison", "Create a graph of comparison confidence intervals", TRUE)),
      conditionalPanel("['kruskal','friedman'].includes(input.method)", checkboxInput("graph_ranks", "Create a graph of ranks", FALSE)),
      conditionalPanel("input.matching == 'matched'", checkboxInput("graph_differences", "Create graphs of within-subject differences", FALSE),
        conditionalPanel("input.comparisons == 'none' || input.comparisons == 'trend'", p(class = "hint", "Difference graphs show all condition pairs when there are no pairwise follow-up tests."))),
      h4("Additional results"),
      checkboxInput("report_descriptive", "Report descriptive statistics for each data set", TRUE),
      conditionalPanel("['anova','rm_anova','mixed','rm_auto','log_rm_anova','log_mixed','log_rm_auto'].includes(input.method)",
        checkboxInput("report_goodness", "Report goodness of fit", FALSE),
        checkboxInput("report_aicc", "Compare condition-means and common-mean models using AICc", FALSE),
        checkboxInput("report_effects", "Report effect sizes", FALSE)),
      h4("Output"), numericInput("output_digits", "Significant digits (except p-values)", 4, min = 2, max = 12),
      selectInput("p_style", "P-value style", c("Prism style: value and significance stars" = "prism", "Decimal" = "decimal", "Scientific notation" = "scientific", "Significance stars only" = "stars")),
      numericInput("p_digits", "P-value decimal places / significant digits", 6, min = 1, max = 15))),
    tabPanel("Residuals", value = "residuals", div(class = "options-tab-body",
      h4("Graphs to create"), checkboxInput("graph_residual", "Residuals versus fitted values", FALSE),
      checkboxInput("graph_homoscedasticity", "Absolute residuals versus fitted values (homoscedasticity)", FALSE),
      checkboxInput("graph_residual_qq", "Residual QQ plot", FALSE), checkboxInput("graph_residual_heatmap", "Residual heatmap", FALSE),
      h4("Diagnostics for residuals"),
      conditionalPanel("!['rm_anova','mixed','rm_auto','log_rm_anova','log_mixed','log_rm_auto','friedman','paired_t','paired_wilcox'].includes(input.method)", checkboxInput("residual_variance", "Test group variance equality: Brown-Forsythe and Bartlett", FALSE)),
      conditionalPanel("['rm_anova','mixed','rm_auto','log_rm_anova','log_mixed','log_rm_auto','friedman','paired_t','paired_wilcox'].includes(input.method)", p(class = "hint", "Independent-group variance tests are unavailable for matched data.")),
      checkboxInput("residual_normality", "Test residual normality: D'Agostino-Pearson, Anderson-Darling, Shapiro-Wilk and corrected Kolmogorov-Smirnov", FALSE),
      p(class = "hint", "Residual tests are exploratory because residuals share fitted parameters. The KS test uses a Lilliefors correction for estimated normal-distribution parameters.")))),
  div(class = "download-row", actionButton("save_analysis_defaults", "Save these choices as my defaults"), actionButton("reset_analysis_defaults", "Reset defaults")),
  p(class = "hint", "Defaults are saved in this browser for future sessions. Column mappings, controls and selected pairs stay specific to the table.")
)
