sheet_column_letters <- function(n) {
  vapply(seq_len(n), function(i) {
    label <- ""
    while (i > 0) { label <- paste0(LETTERS[(i - 1L) %% 26L + 1L], label); i <- (i - 1L) %/% 26L }
    label
  }, character(1))
}
result_sheet_ui <- function(sheet, options = list(), method = "") {
  frame <- sheet_frame(sheet, options); n <- length(sheet$columns)
  div(class = "prism-sheet-wrap", tags$table(class = "prism-result-sheet",
    tags$thead(tags$tr(tags$th(class = "sheet-corner", rowspan = 2, colspan = 2, div(class = "sheet-analysis-title", method), sheet$title), lapply(sheet_column_letters(n), tags$th, class = "sheet-letter")),
      tags$tr(lapply(sheet$headers, function(label) tags$th(class = "sheet-dataset", label)))),
    tags$tbody(lapply(seq_along(sheet$rows), function(i) {
      row <- sheet$rows[[i]]
      tags$tr(class = paste0("sheet-row-", row$kind), tags$th(class = "sheet-row-number", scope = "row", i), tags$th(class = "sheet-statistic", scope = "row", frame[i, 1]),
        lapply(seq_len(n), function(j) tags$td(class = if (j %in% row$highlight) "sheet-highlight" else NULL, frame[i, j + 1])))
    }))), if (nzchar(sheet$note)) p(class = "hint sheet-note", sheet$note))
}
