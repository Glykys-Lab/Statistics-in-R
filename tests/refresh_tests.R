if (!exists("app_env")) {
  if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
  app_env <- new.env(parent = globalenv())
  sys.source("app.R", envir = app_env)
}
shiny::testServer(app_env$server, {
  session$setInputs(table_type = "column", source = "demo", demo = "groups", layout = "wide",
    wide_columns = c("Control", "Drug A", "Drug B"), wide_id = "Subject", family = "many",
    matching = "independent", distribution = "normal", variance = "equal", method = "anova",
    y = "response", group = "group", decimal = ".", comparisons = "all", error_control = "fwer", adjustment = "holm")
  stopifnot("tukey" %in% correction_update()$choices, correction_update()$selected == "holm")
  session$setInputs(comparisons = "control")
  stopifnot(!"tukey" %in% correction_update()$choices, correction_update()$selected == "holm")
  session$setInputs(error_control = "fdr")
  stopifnot(identical(unname(correction_update()$choices), c("BKY", "BH", "BY")), correction_update()$selected == "BKY")
  session$setInputs(adjustment = "BH", error_control = "none")
  stopifnot(identical(unname(correction_update()$choices), "none"), correction_update()$selected == "none")
  session$setInputs(adjustment = "none", error_control = "fwer", method = "mixed")
  stopifnot(!"tukey" %in% correction_update()$choices, correction_update()$selected == "sidak")
  session$setInputs(adjustment = "sidak")
  # Restore the last sent choice after the user changed the input themselves.
  session$setInputs(adjustment = "holm")
  restore_preferences(list(adjustment = "sidak"))
  session$flushReact()
  stopifnot(correction_update()$selected == "sidak", pending_preferences()$adjustment == "sidak")
  session$setInputs(adjustment = "sidak")
  stopifnot(is.null(pending_preferences()$adjustment))
  cat("PASS: Correction choices update across methods and comparison approaches; valid choices and restored defaults are preserved\n")
})
