# Projects store data and settings as JSON, never serialized executable R objects.
# Manage named tables, immutable analyses, saved graphs, layouts and project notes.
# read_project() validates references and reruns saved analyses from their inputs;
# pack_frame/unpack_frame preserve column types and missing cells across JSON.
new_project <- function(name = "Untitled project") list(version = 1L, name = name, info = "", sequence = 0L,
  tables = list(), analyses = list(), graphs = list(), layouts = list())

project_add <- function(project, section, record) {
  # Names may be edited, but stable IDs come from a project-wide sequence so
  # graph and layout references do not depend on the visible sheet name.
  existing <- vapply(project[[section]], function(e) e$name, character(1))
  record$name <- tail(make.unique(c(existing, record$name), sep = " #"), 1)
  project$sequence <- project$sequence + 1L
  id <- paste0(section, "_", project$sequence)
  record$id <- id; project[[section]][[id]] <- record
  list(project = project, id = id)
}

pack_frame <- function(d) list(columns = lapply(names(d), function(n) list(name = n,
  type = if (is.integer(d[[n]])) "integer" else if (is.numeric(d[[n]])) "numeric" else if (is.logical(d[[n]])) "logical" else "character",
  values = as.list(d[[n]]))))

unpack_frame <- function(x) {
  assert(is.list(x$columns) && length(x$columns) > 0 && length(x$columns) <= 1000, "Invalid project table columns.")
  columns <- lapply(x$columns, function(col) {
    assert(col$type %in% c("integer", "numeric", "logical", "character") && is.list(col$values), "Invalid project column.")
    assert(length(col$values) <= 1000000, "Project table exceeds the row limit.")
    convert <- switch(col$type, integer = as.integer, numeric = as.numeric, logical = as.logical, character = as.character)
    vapply(col$values, function(v) { if (is.null(v)) return(switch(col$type, integer = NA_integer_, numeric = NA_real_, logical = NA, character = NA_character_)); assert(length(v) == 1 && !is.list(v), "Invalid project cell."); convert(v) },
      switch(col$type, integer = integer(1), numeric = numeric(1), logical = logical(1), character = character(1)))
  })
  names(columns) <- vapply(x$columns, function(c) c$name, character(1))
  assert(length(unique(lengths(columns))) == 1, "Project columns have different row counts.")
  validate_import(as.data.frame(columns, check.names = FALSE, stringsAsFactors = FALSE))
}

project_json <- function(project) {
  payload <- project
  payload$tables <- lapply(payload$tables, function(t) { t$data <- pack_frame(t$data); t$nested_config <- pack_nested_config(t$nested_config); t })
  payload$analyses <- lapply(payload$analyses, function(a) {
    a$result <- NULL; a$input_data <- pack_frame(a$input_data); a$time <- as.character(a$time)
    a$config <- pack_nested_config(a$config)
    if (!is.null(a$worksheet)) {
      a$worksheet$data <- pack_frame(a$worksheet$data)
      a$worksheet$nested_config <- pack_nested_config(a$worksheet$nested_config)
    }
    a
  })
  jsonlite::toJSON(payload, auto_unbox = TRUE, na = "null", null = "null", digits = NA, pretty = TRUE)
}

read_project <- function(text) {
  # Reconstruct and validate a complete project before the server replaces its
  # current workspace. Saved fits are recalculated, not trusted serialized objects.
  p <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  assert(identical(p$version, 1L) && length(p$name) == 1 && is.character(p$name), "Unsupported LabStats project file.")
  assert(is.numeric(p$sequence) && length(p$sequence) == 1 && p$sequence >= 0, "Invalid project sequence.")
  for (section in c("tables", "analyses", "graphs", "layouts")) {
    assert(is.list(p[[section]]) && length(p[[section]]) <= 2000, "Invalid or oversized project section.")
    for (id in names(p[[section]])) {
      e <- p[[section]][[id]]
      assert(identical(e$id, id) && length(e$name) == 1 && is.character(e$name), "Invalid project entry.")
      suffix <- suppressWarnings(as.integer(sub(".*_", "", id)))
      assert(grepl(paste0("^", section, "_[0-9]+$"), id) && !is.na(suffix) && suffix <= p$sequence, "Invalid project entry ID.")
    }
  }
  p$tables <- lapply(p$tables, function(t) {
    t$data <- unpack_frame(t$data)
    assert(t$type %in% c("column", "xy", "other", "nested") && t$layout %in% c("wide", "long"), "Invalid table type.")
    t$columns <- as.character(unlist(t$columns)); t
  })
  p$analyses <- lapply(p$analyses, function(a) {
    assert(a$table_id %in% names(p$tables) && a$config$method %in% analysis_methods, "Invalid analysis reference.")
    a$input_data <- unpack_frame(a$input_data)
    if (!is.null(a$worksheet)) {
      a$worksheet$data <- unpack_frame(a$worksheet$data)
      a$worksheet$columns <- as.character(unlist(a$worksheet$columns))
    }
    a$result <- run_analysis(a$input_data, a$config); a
  })
  allowed_options <- names(formals(make_plot)); allowed_options <- setdiff(allowed_options, "result")
  p$graphs <- lapply(p$graphs, function(g) {
    assert(g$analysis_id %in% names(p$analyses) && g$kind %in% graph_kinds, "Invalid graph reference.")
    assert(is.list(g$options) && all(names(g$options) %in% allowed_options), "Invalid graph options.")
    for (v in g$options) assert(is.null(v) || (length(v) == 1 && !is.list(v)), "Invalid graph option value.")
    make_graph(p$analyses[[g$analysis_id]]$result, g$kind, g$options); g
  })
  p$layouts <- lapply(p$layouts, function(l) {
    l$graph_ids <- as.character(unlist(l$graph_ids))
    assert(length(l$graph_ids) > 0 && all(l$graph_ids %in% names(p$graphs)) && l$mode %in% c("grid", "overlay"), "Invalid layout reference.")
    assert(length(l$columns) == 1 && l$columns %in% 1:6, "Invalid layout column count."); l
  })
  assert(is.character(p$info) && length(p$info) == 1, "Invalid project information.")
  p
}

normalize_wide <- function(data, columns, baseline, operation = "ratio") {
  # Transform each row against its own baseline; leave undefined divisions or
  # logarithms blank and return an audit so the derived table explains exclusions.
  assert(operation %in% c("ratio", "percent", "change", "percent_change", "log"), "Unknown transformation.")
  assert(length(columns) > 0 && all(columns %in% names(data)), "Select measurement columns to transform.")
  for (n in columns) numeric_column(data, n)
  result <- data; invalid <- integer(length(columns)); names(invalid) <- columns
  if (operation != "log") {
    assert(baseline %in% columns, "Select a measurement column as the row baseline.")
    base <- numeric_column(data, baseline)
  }
  for (n in columns) {
    x <- data[[n]]
    bad <- if (operation == "log") !is.na(x) & x <= 0 else if (operation == "change") rep(FALSE, length(x)) else !is.na(x) & !is.na(base) & base == 0
    invalid[n] <- sum(bad)
    transformed <- switch(operation, ratio = x/base, percent = 100*x/base, change = x-base,
      percent_change = 100*(x-base)/base, log = suppressWarnings(log(x)))
    transformed[bad | !is.finite(transformed)] <- NA_real_; result[[n]] <- transformed
  }
  list(data = result, audit = data.frame(column = columns, invalid_values_blanked = unname(invalid)),
    note = if (operation == "log") "Natural logarithm; nonpositive responses are blank." else paste(operation, "relative to", baseline, "in each subject row; undefined divisions are blank."))
}

cleaned_wide <- function(data, result, columns, id_column = "") {
  assert(result$method == "outliers", "Run outlier identification first.")
  details <- result$tables$outliers
  for (i in seq_len(nrow(details))) {
    row <- (details$source_row[i] - 1L) %% nrow(data) + 1L
    col <- as.character(details$group[i])
    assert(!is.na(row) && col %in% columns, "The outlier snapshot does not match this worksheet.")
    data[row, col] <- NA_real_
  }
  data
}
