# Bundled repeated-cell reference, executed by nested.R rather than sourced by Shiny.
# The .inc extension prevents automatic execution of its top-level model/file code.
# LabStats injects settings, temporary paths, a normality helper and the selected
# control; existing section comments explain the reference calculations and figures.
# ==============================================================================
# Nested mixed model for the R CSVs written by the Python analysis pipeline
#
# Input layout (one row per cell, written by analysis_common.write_r_csv):
#   Mouse | Slice | Cell | <condition 1> | <condition 2> | ...
# Blank cells = missing. Works for any analysis (Area, Area_Change, Z, ...).
#
# Every cell measurement is used. The model has three layers of variability:
#   Mouse / Slice / Cell     : random intercepts (cells nested in slices in mice)
#   Slice x Condition        : each slice may respond differently to each
#                              condition (see slice_condition in SETTINGS)
#   Residual (cell level)    : cell-to-cell noise within a slice and condition,
#                              shared or per condition (see model_engine)
#
# Outputs: console tables + plots, and an Obsidian note in results_dir
# ==============================================================================

# 1. Packages ------------------------------------------------------------------
# Loaded quietly: the "masked from" and welcome messages are expected and
# harmless (lmerTest's lmer is meant to replace lme4's, and adds p-values).
# Real problems still show: errors, and warnings such as convergence.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)   # p-values via Satterthwaite approximation
  library(emmeans)    # estimated means and pairwise comparisons
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
})
# nlme is called as nlme::lme() rather than attached with library(), so it
# can't mask anything from lme4. It only needs to be installed.
stopifnot("Package 'nlme' is not installed" = requireNamespace("nlme", quietly = TRUE))

# Use the same df method as summary(model) so p-values agree between outputs.
emm_options(lmer.df = "satterthwaite")

# 2. SETTINGS -- the only lines you should need to edit ------------------------
r_csv_dir <- "H:/Obsidian Vault/Glykys Labratory/Lab Notebook/Analysis/R CSVs"
csv_name  <- "z_changes_ROIs_HIC-NaCl-Adult.csv"

# POP-UPS: TRUE = when run in RStudio, the script asks for the key choices
# (log scale, compare to 0, which comparisons, outliers) in pop-up windows,
# using the values below as the defaults, and lets you pick the CSV or a folder
# if a path below does not exist. FALSE = use the values below as they are.
interactive_prompts <- TRUE
# TRUE: each run starts with a file window opened at r_csv_dir -- pick the CSV
# to analyse there (cancel = use csv_name above). Needs interactive_prompts = TRUE.
pick_csv_each_run   <- TRUE

# DECISION 1 -- data scale
# TRUE : log10-transform the values. Use when the data are right-skewed or the
#        variance grows with the mean (check the residual plots first).
#        Requires all values > 0; the script stops if it finds zeros or negatives.
# FALSE: analyze the raw values.
log_transform <- FALSE

# OUTLIERS -- ROUT (Motulsky & Brown 2006, BMC Bioinformatics 7:123), the method
# Prism uses. Applied within each condition.
#   "flag"  : find outliers, mark them in the notebook plot and the note, keep them
#   "remove": find outliers and remove them from the analysis (you confirm first)
#   "none"  : no outlier test
# Decide this before looking at the results, and report it in the Methods.
outlier_action <- "flag"
outlier_unit   <- "cell"       # "cell", "slice" (flagged only, never removed automatically) or "both"
outlier_Q      <- 0.01         # ROUT's Q: the false discovery rate (Prism's default is 1%)
outlier_level  <- "residuals"  # "residuals": after removing each mouse's and slice's offset, so a
#   cell is judged against its own slice (respects the nesting)
# "raw": the raw values in each condition (like Prism's columns)

# DECISION 3 -- slice x condition (decision 2, the Mouse/Slice/Cell nesting, is
# always on and set automatically)
# Slice x condition random effect: lets each slice respond to each condition
# differently. Without it, every cell is treated as independent evidence about
# the drug effect even though cells in one slice share that slice's response.
#   "shared"        : one slice-to-slice SD for all conditions.
#   "per_condition" : each condition gets its own slice-to-slice SD (use when
#                     slices vary much more in one condition, e.g. washout).
#   "auto"          : always includes the term; uses "per_condition" only if a
#                     likelihood-ratio test gives p < variance_alpha.
#   "none"          : leave it out (the old model -- only if slices respond alike).
slice_condition <- "auto"

# DECISION 4 (and 5) -- cell-level noise; this also picks the software,
# because only nlme can give each condition its own cell-level SD.
# Cell-level (residual) variance:
#   "lmer" : one cell-to-cell SD shared by all conditions (lme4 + lmerTest).
#   "nlme" : each condition gets its own cell-to-cell SD (nlme::lme + varIdent).
#   "auto" : uses "nlme" only if a likelihood-ratio test gives p < variance_alpha.
model_engine <- "auto"

# Threshold for both "auto" choices above. Pick the settings BEFORE looking at
# the condition comparisons, and don't switch afterwards for nicer p-values.
variance_alpha <- 0.05

# TRUE: print (and add to the note) the variance tables and tests for all
# options, even the ones not used.
check_variance_homogeneity <- TRUE

# FALSE: each condition is compared to the FIRST condition column (reference).
# TRUE : each condition's mean is tested against 0 (e.g. change-from-baseline).
#        This does not compare conditions with each other -- see the
#        "Pairwise comparisons" table for that.
compare_to_zero <- FALSE

# DECISION 6 -- degrees of freedom (Satterthwaite, always from lmer)
#   "per_comparison": each comparison gets its own df -- the field standard
#                     (lmerTest/emmeans, SAS, Prism mixed models all do this).
#   "smallest"      : every comparison uses the smallest df -- more cautious.
df_method <- "per_comparison"

# DECISION 7 -- comparisons and correction
# Which comparisons BETWEEN conditions to test (and correct for):
#   "all_pairs"   : every condition vs every other        (Tukey correction)
#   "vs_first"    : every condition vs the first column   (Dunnett correction)
#   "consecutive" : each condition vs the one before it   (Holm correction)
#   "custom"      : only the pairs listed in custom_pairs (Holm correction)
# Fewer comparisons = lighter correction = more power for each one. Decide which
# comparisons your question needs BEFORE looking at the results.
comparison_type <- "all_pairs"

# Only used when comparison_type = "custom". Each c(A, B) tests A minus B
# (A / B when log_transform is TRUE), using condition names exactly as in the CSV.
custom_pairs <- list()
# e.g. custom_pairs <- list(c("Drug 1", "Baseline"), c("Wash", "Drug 1"))

# "auto" uses the correction named above. Alternatives: "none" (uncorrected --
# only defensible for one or two pre-specified comparisons), "holm", "bonferroni".
p_adjust <- "auto"

# Draw p-values on the plot, but only for comparisons below p_threshold.
show_p_on_plot <- TRUE
p_threshold    <- 0.05

# Comparisons panel: each comparison's difference with its adjusted 95% CI --
# the intervals the p-values actually come from.
#   "separate" : saved as its own image (<csv name>_comparisons.png)
#   "none"     : not drawn
comparison_plot <- "separate"

# Decimal places shown in every table, the note and the plot.
num_decimals <- 2
# p-values: below 0.0001 always shows "<0.0001"; a p that rounds to 0 at this
# many places shows as "<0.001". 3 places keeps e.g. p = 0.0511 as "0.051".
p_decimals   <- 3

# What the black marker + bars on the plot show:
#   "model"      : model-estimated mean +/- 95% CI
#                  (geometric mean +/- 95% CI when log_transform = TRUE)
#   "median_iqr" : median with interquartile range, computed from the raw data
#                  (descriptive only -- see note in section 9)
summary_display <- "model"

y_label          <- NULL       # y-axis label; NULL = first word of the CSV name
show_diagnostics <- TRUE       # residual plots to check model assumptions
save_outputs     <- TRUE       # write the note and plots
show_plots       <- FALSE      # TRUE = also draw the plots in RStudio's Plots pane (they are all in the Viewer report)
html_report      <- TRUE       # show the note, with every plot, in RStudio's Viewer to flip through
#   (temporary -- not saved; needs: install.packages("commonmark"))
results_dir      <- "H:/Obsidian Vault/Glykys Labratory/Lab Notebook/Figures"
plot_dir         <- "H:/Obsidian Vault/Glykys Labratory/Admin/Attachments"
experiments_dir  <- "H:/Obsidian Vault/Glykys Labratory/Lab Notebook/Experiments"  # mouse and slice notes (linked in the outlier tables)

# PLOT AXES (all plots) --------------------------------------------------------
# Values that are changes (e.g. Z change, area change) get a dashed line at 0.
#   "auto": yes if the CSV name contains "change" (or compare_to_zero is TRUE)
change_measure <- "auto"   # "auto", TRUE or FALSE
# TRUE: the y axis always includes 0; in the publication figure it also starts
# and ends on a labelled tick, so readers see where it begins. Not used with log_transform.
y_include_zero <- TRUE

# PUBLICATION FIGURE (for Illustrator) -----------------------------------------
# ONE file holding the SuperPlot (cells + slice means, colour = mouse, model
# mean and CI on top of the data) with a small key to its right, and a compact
# comparisons strip below. Drawn at the final printed size: place it in Illustrator
# at 100% so the sizes below are exact. Strokes and text are in points (pt).
make_pub_figures   <- TRUE
pub_formats        <- c("svg")  # "svg": every dot, line and label is its own editable object
#   (needs: install.packages("svglite")); "pdf" and "eps" also possible
pub_dir            <- NULL       # NULL = same folder as plot_dir (inside the vault, so the note can link it)
pub_width_mm       <- 40         # figure width, key included (two fit side by side in an 8 cm column)
pub_key_width_mm   <- NULL       # key to the right of the SuperPlot: NULL = just wide enough for
#   its labels, a number = that many mm, 0 = no key
pub_key_nudge_mm   <- 1          # how far the key reaches over the SuperPlot's empty right edge
#   (bigger = closer to the data; 0 = no overlap)
pub_key_labels     <- c(cell = "Cell", slice = "Slice", mouse = "Colored, Mouse")  # key text, top to bottom
pub_key_pt         <- 4          # key text
pub_main_height_mm <- 30         # SuperPlot (+ key)
pub_comp_height_mm <- 10         # comparisons strip at the bottom (0 = none)
pub_font           <- "Arial"
pub_title_pt       <- 8          # y-axis title (bold)
pub_tick_pt        <- 7          # numbers and labels on the SuperPlot's axes
pub_p_pt           <- 5          # p-values
pub_small_pt       <- 5          # the comparisons strip's labels and numbers
pub_axis_pt        <- 0.25       # x and y axes, ticks, p-value brackets, dashed reference lines
pub_dot_pt         <- 0.25       # outline of the cell dots and slice-mean circles
pub_ci_pt          <- 0.5        # CI lines
pub_mean_pt        <- 1          # mean line
pub_slice_line_pt  <- 0.25       # lines joining the same slice across conditions
pub_cell_size      <- 0.8        # cell dot size (ggplot size units, about mm)
pub_slice_size     <- 1.8        # slice-mean circle size
pub_cell_alpha     <- 0.4        # cell dot opacity (0 = invisible, 1 = solid)
pub_cell_layout    <- "beeswarm" # "beeswarm" (needs the ggbeeswarm package) or "jitter"
pub_summary        <- "mean_ci"  # "mean_ci" (mean line + error bars), "mean" (line only) or "none"
pub_error_bar      <- "ci"       # error bars: "ci" (model 95% CI) or "sem" (model SEM, +/- 1 SE)
pub_summary_offset <- 0          # 0 = mean and error bars on top of the data; e.g. 0.3 = to the right
pub_mean_width     <- 0.25       # length of the mean line (1 = the gap between conditions)
pub_ci_cap_width   <- 0.1        # width of the error-bar end caps (0 = no caps)
# Titles in this figure only (the note and notebook plot keep the CSV's names).
# With pop-ups on, a window offers to change them each run (pub_ask_titles).
pub_ask_titles       <- TRUE
pub_micro            <- TRUE     # TRUE: a "u" before a unit (uM, um, umol, uL, us, uA, uV, ug, uF, uW,
#   uHz) becomes the micro sign, e.g. "10 uM" -> "10 μM"
pub_ion_superscript  <- TRUE     # TRUE: ion charges are raised: "Ca2+", "Cl-", "[K+]o" -> charge as superscript
pub_dash             <- "22"     # dashed lines: "22" = 2 pt dash, 2 pt gap; "44" = 4 pt (R's "dashed")
pub_y_title          <- NULL     # NULL = y_label (with "(log scale)" when logged)
pub_x_title          <- NULL     # NULL = no x-axis title
pub_condition_labels <- NULL     # e.g. c("NaCl" = "Control", "Bumetanide" = "Bumet."); "\n" = line break
pub_mouse_colours  <- NULL       # NULL = colour-blind-safe Okabe-Ito colours (section 9),
# or e.g. c("2025.05.21" = "#0072B2", ...)

# LabStats supplies validated settings and temporary input/output paths here.
if (exists("labstats_overrides", inherits = FALSE)) list2env(labstats_overrides, envir = environment())

# 2b. Pop-up helpers and path check ---------------------------------------------
ask <- interactive_prompts && interactive()
has_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()

# Yes/No question in a pop-up (RStudio), or a small list window otherwise.
ask_yes_no <- function(title, question, default = TRUE) {
  if (!ask) return(default)
  if (has_rstudio) {
    return(isTRUE(rstudioapi::showQuestion(title, question, ok = "Yes", cancel = "No")))
  }
  pick <- utils::select.list(c("Yes", "No"), preselect = if (default) "Yes" else "No",
                             title = paste(title, "-", question), graphics = TRUE)
  if (length(pick) == 0 || identical(pick, "")) default else pick == "Yes"
}

# Pick one (or several) options from a list window. Closing the window keeps the default.
ask_choice <- function(title, choices, default = 1, multiple = FALSE) {
  if (!ask) return(choices[default])
  message("\n", title)   # also in the console, in case the window title is cut short
  pick <- utils::select.list(choices, preselect = choices[default], multiple = multiple,
                             title = title, graphics = TRUE)
  if (length(pick) == 0 || identical(pick, "")) {
    message("Nothing chosen -- using: ", paste(choices[default], collapse = ", "))
    return(choices[default])
  }
  pick
}

# Type a short text in a pop-up (RStudio), or in the console otherwise.
# Cancel keeps the default.
ask_text <- function(title, message, default = "") {
  if (!ask) return(default)
  if (has_rstudio) {
    x <- rstudioapi::showPrompt(title, message, default)
    return(if (is.null(x)) default else x)
  }
  x <- readline(paste0(message, " [", default, "]: "))
  if (identical(x, "")) default else x
}

pick_file <- function(caption, start_dir) {
  if (!dir.exists(start_dir)) start_dir <- getwd()
  if (has_rstudio) return(rstudioapi::selectFile(caption = caption, path = start_dir, filter = "CSV files (*.csv)"))
  if (.Platform$OS.type == "windows") {
    f <- utils::choose.files(default = file.path(start_dir, "*.csv"), caption = caption, multi = FALSE,
                             filters = matrix(c("CSV files (*.csv)", "*.csv"), 1))
    return(if (length(f) == 0) NULL else f)
  }
  file.choose()
}

pick_dir <- function(caption, start_dir) {
  if (!dir.exists(start_dir)) start_dir <- getwd()
  if (has_rstudio) return(rstudioapi::selectDirectory(caption = caption, path = start_dir))
  if (.Platform$OS.type == "windows") {
    d <- utils::choose.dir(default = start_dir, caption = caption)
    return(if (is.na(d)) NULL else d)
  }
  NULL
}

run_choices <- character(0)   # what was picked in the pop-ups, listed in the note

# The CSV: pick it in a file window (pick_csv_each_run), or if the path in
# SETTINGS does not exist.
csv_path <- file.path(r_csv_dir, csv_name)
if (ask && pick_csv_each_run) {
  picked <- pick_file("Choose the CSV to analyse (Cancel = use csv_name in SETTINGS)", r_csv_dir)
  if (is.null(picked) || length(picked) == 0 || identical(picked, "")) {
    message("No CSV chosen -- using csv_name from SETTINGS: ", csv_name)
  } else {
    r_csv_dir <- dirname(picked)
    csv_name  <- basename(picked)
    csv_path  <- picked
    run_choices <- c(run_choices, paste0("CSV chosen in a pop-up: ", csv_name))
  }
}
if (!file.exists(csv_path)) {
  if (!ask) stop("CSV not found: ", csv_path, " -- fix r_csv_dir / csv_name in SETTINGS.")
  message("CSV not found: ", csv_path, " -- choose it in the window that opens.")
  picked <- pick_file("Choose the CSV to analyse", r_csv_dir)
  if (is.null(picked) || length(picked) == 0 || identical(picked, "")) stop("No CSV chosen.")
  r_csv_dir <- dirname(picked)
  csv_name  <- basename(picked)
  run_choices <- c(run_choices, paste0("CSV chosen in a pop-up: ", csv_name))
}

# Output folders: pick another if one does not exist (cancel = create it).
for (dir_var in c("results_dir", "plot_dir")) {
  d <- get(dir_var)
  if (!dir.exists(d)) {
    label  <- c(results_dir = "the Obsidian note", plot_dir = "the images")[[dir_var]]
    picked <- NULL
    if (ask) {
      message("Folder for ", label, " not found: ", d, " -- choose one (or cancel to create it).")
      picked <- pick_dir(paste("Folder for", label), dirname(d))
    }
    if (is.null(picked) || identical(picked, "")) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      message("Created ", d)
    } else {
      assign(dir_var, picked)
      run_choices <- c(run_choices, paste0("Folder for ", label, ": ", picked))
    }
  }
}
if (is.null(pub_dir)) pub_dir <- plot_dir
if (is.null(y_label)) y_label <- strsplit(tools::file_path_sans_ext(csv_name), "_")[[1]][1]

# 3. Import --------------------------------------------------------------------
# check.names = FALSE keeps condition names exactly as written
# ("Drug 1" would otherwise become "Drug.1"). fileEncoding = "UTF-8-BOM" removes
# the invisible byte-order mark Excel puts at the start of some CSVs
# (otherwise the first column reads as "ï..Mouse").
df <- read.csv(file.path(r_csv_dir, csv_name), check.names = FALSE, fileEncoding = "UTF-8-BOM")

id_cols        <- c("Mouse", "Slice", "Cell")
condition_cols <- setdiff(names(df), id_cols)

stopifnot(
  "Mouse/Slice/Cell columns not found -- is this an R CSV?" = all(id_cols %in% names(df)),
  "No condition columns found" = length(condition_cols) >= 1,
  "summary_display must be 'model' or 'median_iqr'" = summary_display %in% c("model", "median_iqr"),
  "comparison_type must be all_pairs, vs_first, consecutive or custom" =
    comparison_type %in% c("all_pairs", "vs_first", "consecutive", "custom"),
  "p_adjust must be auto, none, holm or bonferroni" = p_adjust %in% c("auto", "none", "holm", "bonferroni"),
  "model_engine must be lmer, nlme or auto" = model_engine %in% c("lmer", "nlme", "auto"),
  "df_method must be per_comparison or smallest" = df_method %in% c("per_comparison", "smallest"),
  "pub_summary must be mean_ci, mean or none" = pub_summary %in% c("mean_ci", "mean", "none"),
  "change_measure must be \"auto\", TRUE or FALSE" = identical(change_measure, "auto") || isTRUE(change_measure) || isFALSE(change_measure),
  "outlier_action must be flag, remove or none" = outlier_action %in% c("flag", "remove", "none"),
  "outlier_unit must be cell, slice or both" = outlier_unit %in% c("cell", "slice", "both"),
  "outlier_level must be residuals or raw" = outlier_level %in% c("residuals", "raw"),
  "pub_error_bar must be ci or sem" = pub_error_bar %in% c("ci", "sem"),
  "pub_cell_layout must be beeswarm or jitter" = pub_cell_layout %in% c("beeswarm", "jitter"),
  "pub_formats must be pdf, svg and/or eps" = all(pub_formats %in% c("pdf", "svg", "eps")),
  "comparison_plot must be separate or none" = comparison_plot %in% c("separate", "none"),
  "slice_condition must be shared, per_condition, auto or none" =
    slice_condition %in% c("shared", "per_condition", "auto", "none")
)

if (log_transform && compare_to_zero) {
  warning("log_transform + compare_to_zero tests the GEOMETRIC MEAN against 1 ",
          "(log10 = 0), which is rarely what you want. Set one of them to FALSE.")
}

# 4. Wide -> long --------------------------------------------------------------
# Conditions are detected from the file, not hard-coded. values_drop_na removes
# the blank (missing) entries. Factor levels follow column order, so the first
# condition column is the reference level.
df_long <- df |>
  pivot_longer(
    cols           = all_of(condition_cols),
    names_to       = "Condition",
    values_to      = "Value",
    values_drop_na = TRUE
  ) |>
  mutate(
    Mouse     = as.factor(Mouse),
    Slice     = as.factor(Slice),
    Cell      = as.factor(Cell),
    Condition = factor(Condition, levels = condition_cols)
  ) |>
  droplevels() |>
  as.data.frame()

# 4b. Choices for this dataset (pop-ups) ----------------------------------------
src_log <- "Fixed by you in SETTINGS"
src_cmp <- "Fixed by you in SETTINGS"

# Decision 1: log scale -- only offered when every value is above 0.
if (any(df_long$Value <= 0)) {
  if (log_transform) warning("log_transform is TRUE but the data has zeros or negatives -- analysing raw values.")
  log_transform <- FALSE
  src_log <- "Log not possible: the data has zeros or negatives"
  message("Log scale skipped: the data has zeros or negatives.")
} else if (ask) {
  v    <- df_long$Value
  skew <- mean((v - mean(v))^3) / sd(v)^3
  log_transform <- ask_yes_no("Decision 1: data scale", paste0(
    "Log10-transform the values?\n\nUse log when the values are right-skewed or their spread grows ",
    "with the mean. Skewness of these data: ", round(skew, 2), " (above about 1 = clearly right-skewed)."),
    default = log_transform)
  src_log <- "Chosen in a pop-up when running"
  run_choices <- c(run_choices, paste0("Log10 transform: ", if (log_transform) "yes" else "no"))
}

# Test each condition's mean against 0?
if (ask) {
  compare_to_zero <- ask_yes_no("Compare to 0", paste0(
    "Also test each condition's mean against 0?\n\nUse this for change values, where 0 = no change."),
    default = compare_to_zero || grepl("change", csv_name, ignore.case = TRUE))
  run_choices <- c(run_choices, paste0("Test each mean against 0: ", if (compare_to_zero) "yes" else "no"))
}

# Decision 7: which comparisons (and, for "vs control", which condition is the control).
cond_now <- levels(df_long$Condition)
if (ask && length(cond_now) >= 2) {
  cmp_opts <- c("All pairs (Tukey)"                           = "all_pairs",
                "Each condition vs one control (Dunnett)"      = "vs_first",
                "Each condition vs the one before it (Holm)"   = "consecutive",
                "Only the pairs I pick (Holm)"                 = "custom")
  default_i <- match(comparison_type, cmp_opts)
  cmp_pick  <- ask_choice("Decision 7: which comparisons?", names(cmp_opts),
                          default = if (is.na(default_i)) 1 else default_i)
  comparison_type <- unname(cmp_opts[cmp_pick])
  cmp_text <- cmp_pick
  if (comparison_type == "vs_first") {
    control <- ask_choice("Which condition is the control?", cond_now, default = 1)
    df_long$Condition <- factor(df_long$Condition, levels = c(control, setdiff(cond_now, control)))
    cmp_text <- paste0(cmp_text, "; control: ", control)
  }
  if (comparison_type == "custom") {
    pair_grid  <- t(utils::combn(cond_now, 2))              # columns: earlier, later
    pair_names <- paste(pair_grid[, 2], "vs", pair_grid[, 1])
    picked     <- ask_choice("Pick the comparisons to test", pair_names,
                             default = seq_along(pair_names), multiple = TRUE)
    custom_pairs <- lapply(match(picked, pair_names), function(i) c(pair_grid[i, 2], pair_grid[i, 1]))
    cmp_text <- paste0(cmp_text, ": ", paste(picked, collapse = ", "))
  }
  src_cmp <- "Chosen in a pop-up when running"
  run_choices <- c(run_choices, paste0("Comparisons: ", cmp_text))
}
cond_levels <- levels(df_long$Condition)

# Outliers: flag, remove or skip.
if (ask) {
  out_opts <- c("Flag outliers but keep them (ROUT)" = "flag",
                "Remove outliers (ROUT) -- you confirm the list first" = "remove",
                "No outlier test" = "none")
  default_i  <- match(outlier_action, out_opts)
  out_pick   <- ask_choice(paste0("Outliers: ROUT test, Q = ", outlier_Q * 100, "%"), names(out_opts),
                           default = if (is.na(default_i)) 1 else default_i)
  outlier_action <- unname(out_opts[out_pick])
  run_choices <- c(run_choices, paste0("Outliers: ", out_pick))
  if (outlier_action != "none") {
    unit_opts <- c("Cells" = "cell",
                   "Slices (flagged only, never removed)" = "slice",
                   "Cells and slices (slices flagged only)" = "both")
    default_i <- match(outlier_unit, unit_opts)
    unit_pick <- ask_choice("Outliers: test cells, slices or both?", names(unit_opts),
                            default = if (is.na(default_i)) 1 else default_i)
    outlier_unit <- unname(unit_opts[unit_pick])
    run_choices  <- c(run_choices, paste0("Outlier test on: ", unit_pick))
  }
}

# Illustrator figure: axis titles and condition names. These change only the
# publication figure; the note, tables and notebook plot keep the CSV's names.
if (ask && make_pub_figures && pub_ask_titles) {
  pub_cond_now <- levels(df_long$Condition)
  y_now   <- if (is.null(pub_y_title)) y_label else pub_y_title
  x_now   <- if (is.null(pub_x_title)) "" else pub_x_title
  lab_now <- setNames(pub_cond_now, pub_cond_now)
  if (!is.null(pub_condition_labels)) {
    hit <- pub_cond_now %in% names(pub_condition_labels)
    lab_now[hit] <- pub_condition_labels[pub_cond_now[hit]]
  }
  if (ask_yes_no("Illustrator figure", paste0(
    "Change the axis titles or condition names in the Illustrator figure?\n\n",
    "y axis: ", y_now, "\nx axis: ", if (x_now == "") "(none)" else x_now,
    "\nConditions: ", paste(lab_now, collapse = ", ")), default = FALSE)) {
    pub_y_title <- ask_text("Illustrator figure: y axis", "y-axis title:", y_now)
    pub_x_title <- ask_text("Illustrator figure: x axis", "x-axis title (leave empty for none):", x_now)
    for (cn in pub_cond_now) {
      lab_now[cn] <- ask_text("Illustrator figure: condition name",
                              paste0("Name shown for \"", cn, "\" (type \\n for a line break):"), lab_now[[cn]])
    }
    pub_condition_labels <- lab_now
    run_choices <- c(run_choices, paste0(
      "Illustrator figure titles: y = \"", pub_y_title, "\", x = \"", pub_x_title, "\", conditions = ",
      paste0(names(lab_now), " -> ", lab_now, collapse = "; ")))
  }
}
if (!is.null(pub_y_title))          pub_y_title          <- gsub("\\n", "\n", pub_y_title, fixed = TRUE)
if (!is.null(pub_x_title))          pub_x_title          <- gsub("\\n", "\n", pub_x_title, fixed = TRUE)
if (!is.null(pub_condition_labels)) pub_condition_labels <- gsub("\\n", "\n", pub_condition_labels, fixed = TRUE)

# 5. Transform -----------------------------------------------------------------
if (log_transform) {
  if (any(df_long$Value <= 0)) {
    stop("log_transform = TRUE but the data has zeros or negatives. ",
         "Set log_transform <- FALSE (e.g. for change-from-baseline data).")
  }
  df_long$Response <- log10(df_long$Value)
} else {
  df_long$Response <- df_long$Value
}

multi       <- nlevels(df_long$Condition) >= 2
cond_levels <- levels(df_long$Condition)

# 5b. Outliers (ROUT) ------------------------------------------------------------
# ROUT (Motulsky & Brown 2006): (1) a robust fit that outliers barely pull on --
# here a robust centre (Lorentzian / Cauchy M-estimate); (2) the robust SD of
# the residuals (RSDR = 68.27th percentile of |residual|, scaled by n/(n-1));
# (3) each point's t = |residual| / RSDR and its p-value; (4) a false-discovery
# test, from the largest residual down: the i-th largest is an outlier if
# p < Q * (n - (i - 1)) / n; testing stops at the first point that passes.
# This follows the published method; results can differ slightly from Prism.
rout_test <- function(x, Q = 0.01) {
  n <- length(x)
  none <- list(outlier = rep(FALSE, n), t = rep(NA_real_, n), p = rep(NA_real_, n))
  if (n < 4) return(none)
  rsdr_of <- function(r) stats::quantile(abs(r), 0.6827, names = FALSE) * n / (n - 1)
  centre <- stats::median(x)
  for (iter in 1:200) {                       # iteratively reweighted robust centre
    r    <- x - centre
    rsdr <- rsdr_of(r)
    if (!is.finite(rsdr) || rsdr == 0) break
    w   <- 1 / (1 + (r / rsdr)^2)             # Lorentzian (Cauchy) weights
    new <- sum(w * x) / sum(w)
    done <- abs(new - centre) < 1e-10 * (1 + abs(centre))
    centre <- new
    if (done) break
  }
  r    <- x - centre
  rsdr <- rsdr_of(r)
  if (!is.finite(rsdr) || rsdr == 0) return(none)
  t_val <- abs(r) / rsdr
  p_val <- 2 * stats::pt(t_val, df = n - 1, lower.tail = FALSE)
  out   <- rep(FALSE, n)
  ord   <- order(t_val, decreasing = TRUE)
  for (i in seq_len(n)) {
    if (p_val[ord[i]] < Q * (n - (i - 1)) / n) out[ord[i]] <- TRUE else break
  }
  list(outlier = out, t = t_val, p = p_val)
}

df_long$Outlier   <- FALSE
outlier_tbl       <- slice_outlier_tbl <- flagged_slice_pts <- NULL
outliers_removed  <- 0
outlier_basis     <- slice_basis <- ""
test_cells        <- outlier_action != "none" && outlier_unit %in% c("cell", "both")
test_slices       <- outlier_action != "none" && outlier_unit %in% c("slice", "both")
outlier_note       <- if (test_cells) "" else "Cells were not tested for outliers."
slice_outlier_note <- if (test_slices) "" else "Slices were not tested for outliers."
rout_methods      <- ""

# The model used to judge each cell against its slice, and each slice against
# its mouse (outlier_level = "residuals").
pre_fit <- NULL
if ((test_cells || test_slices) && outlier_level == "residuals") {
  cr <- any(duplicated(df_long[, c("Mouse", "Slice", "Cell")]))
  out_formula <- as.formula(paste0(
    "Response ~ ", if (nlevels(df_long$Condition) > 1) "Condition" else "1",
    " + (1 | ", if (cr) "Mouse/Slice/Cell" else "Mouse/Slice", ")",
    if (nlevels(df_long$Condition) > 1) " + (1 | Mouse:Slice:Condition)" else ""))
  pre_fit <- tryCatch(
    suppressMessages(suppressWarnings(lmer(out_formula, data = df_long, control = lmerControl(optimizer = "bobyqa")))),
    error = function(e) NULL)
  if (is.null(pre_fit)) message("The outlier model could not be fitted -- testing raw values instead.")
}

# --- Slices: flagged for review, never removed automatically -------------------
# Few slices per condition make this test weak, and an unusual slice is often
# real biology (e.g. a non-responding slice). Exclude a slice only for a
# technical reason decided in advance (slice health, drift, failed imaging).
if (test_slices) {
  sl <- df_long |>
    group_by(Mouse, Slice, Condition) |>
    summarise(n_cells = n(),
              Slice_mean = if (log_transform) 10^mean(log10(Value)) else mean(Value),
              resp_mean  = mean(Response), .groups = "drop") |>
    as.data.frame()
  sl_vals     <- sl$resp_mean
  slice_basis <- "each slice's mean in each condition"
  if (!is.null(pre_fit)) {
    re <- lme4::ranef(pre_fit)
    # Each slice's own response to each condition: its slice-by-condition effect.
    grp <- re[["Mouse:Slice:Condition"]]
    key <- paste(sl$Mouse, sl$Slice, sl$Condition, sep = ":")
    slice_basis_model <- "each slice's own response to each condition (its slice-by-condition effect in the model)"
    if (!is.null(grp)) {
      eff <- grp[key, "(Intercept)"]
      if (all(!is.na(eff))) {
        sl_vals     <- eff
        slice_basis <- slice_basis_model
      }
    }
  }
  sl$Outlier <- FALSE
  sl$Robust_t <- sl$p <- NA_real_
  for (idx in split(seq_len(nrow(sl)), sl$Condition)) {
    res <- rout_test(sl_vals[idx], outlier_Q)
    sl$Outlier[idx]  <- res$outlier
    sl$Robust_t[idx] <- res$t
    sl$p[idx]        <- res$p
  }
  n_sl <- sum(sl$Outlier)
  too_few <- names(which(table(sl$Condition) < 4))
  if (n_sl > 0) {
    slice_outlier_tbl <- data.frame(
      Mouse = as.character(sl$Mouse), Slice = as.character(sl$Slice), Condition = as.character(sl$Condition),
      n_cells = sl$n_cells, Slice_mean = sl$Slice_mean, Robust_t = sl$Robust_t, p = sl$p,
      Action = "Flagged, kept -- check imaging notes"
    )[sl$Outlier, , drop = FALSE]
    rownames(slice_outlier_tbl) <- NULL
    cat("\nROUT slice outliers (flagged only, Q = ", outlier_Q * 100, "%):\n", sep = "")
    print(slice_outlier_tbl, row.names = FALSE)
  }
  slice_outlier_note <- paste0(
    "ROUT (Q = ", outlier_Q * 100, "%) on ", slice_basis, ": ",
    if (n_sl == 0) "no slices flagged."
    else paste0(n_sl, " slice(s) flagged for review and kept in the analysis. Check them against your imaging ",
                "notes; exclude a slice only for a technical reason decided in advance."),
    if (length(too_few) > 0) paste0(" Not tested (fewer than 4 slices): ", paste(too_few, collapse = ", "), ".") else ""
  )
}

# --- Cells: flagged, or removed after you confirm --------------------------------
if (test_cells) {
  vals <- df_long$Response
  outlier_basis <- "the values within each condition"
  if (!is.null(pre_fit)) {
    vals <- as.numeric(resid(pre_fit))
    outlier_basis <- "the residuals within each condition, after removing each mouse's and slice's offset"
  }
  rout_t <- rout_p <- rep(NA_real_, nrow(df_long))
  for (idx in split(seq_len(nrow(df_long)), df_long$Condition)) {
    res <- rout_test(vals[idx], outlier_Q)
    df_long$Outlier[idx] <- res$outlier
    rout_t[idx] <- res$t
    rout_p[idx] <- res$p
  }
  n_out <- sum(df_long$Outlier)
  
  if (n_out > 0) {
    outlier_tbl <- data.frame(
      Mouse = as.character(df_long$Mouse), Slice = as.character(df_long$Slice), Cell = as.character(df_long$Cell),
      Condition = as.character(df_long$Condition), Value = df_long$Value, Robust_t = rout_t, p = rout_p
    )[df_long$Outlier, , drop = FALSE]
    rownames(outlier_tbl) <- NULL
    cat("\nROUT cell outliers (Q = ", outlier_Q * 100, "%):\n", sep = "")
    print(outlier_tbl, row.names = FALSE)
    
    if (outlier_action == "remove") {
      confirm <- ask_yes_no("Remove outliers?", paste0(
        "ROUT found ", n_out, " cell outlier(s); they are listed in the console. Remove them from the analysis?"),
        default = TRUE)
      if (confirm) {
        df_long <- droplevels(df_long[!df_long$Outlier, , drop = FALSE])
        outliers_removed   <- n_out
        outlier_tbl$Action <- "Removed"
      } else {
        outlier_tbl$Action <- "Kept (removal declined)"
        outlier_action     <- "flag"
      }
    } else {
      outlier_tbl$Action <- "Flagged, kept"
    }
  }
  
  outlier_note <- paste0(
    "ROUT (Q = ", outlier_Q * 100, "%) on ", outlier_basis, ": ",
    if (n_out == 0) "no outliers found."
    else if (outliers_removed > 0) paste0(n_out, " outlier(s) found and removed from the analysis.")
    else paste0(n_out, " outlier(s) found and kept in the analysis (open diamonds in the notebook plot).")
  )
}

# Positions of flagged slices for the notebook plot (open squares).
if (!is.null(slice_outlier_tbl)) {
  flagged_slice_pts <- data.frame(
    Condition = factor(slice_outlier_tbl$Condition, levels = levels(df_long$Condition)),
    Value     = slice_outlier_tbl$Slice_mean)
}

# Methods sentence for the writing template.
if (test_cells) {
  rout_methods <- paste0(
    "Outliers were identified with the ROUT method (Q = ", outlier_Q * 100, "%; Motulsky & Brown, 2006) applied to ",
    outlier_basis, "; ",
    if (sum(df_long$Outlier) + outliers_removed == 0) "none were found. "
    else if (outliers_removed > 0) paste0(outliers_removed, " measurement(s) were excluded. ")
    else paste0(sum(df_long$Outlier), " were identified and retained. ")
  )
}
if (test_slices) {
  rout_methods <- paste0(
    rout_methods,
    "Slices were screened with ROUT (Q = ", outlier_Q * 100, "%) on ", slice_basis, " and flagged for review; ",
    "none were excluded on this basis. "
  )
}

# 6. Model building blocks -----------------------------------------------------
# Random intercepts: cells nested in slices nested in mice. Cell is only kept
# if a cell is actually measured more than once (repeated measures);
# otherwise it is not estimable. A cell is identified by Mouse + Slice + Cell,
# so "Cell 1" in two different slices is not counted as a repeat.
cells_repeated <- any(duplicated(df_long[, c("Mouse", "Slice", "Cell")]))
base_groups    <- if (cells_repeated) "Mouse/Slice/Cell" else "Mouse/Slice"

# Fixed effects -- what the model's coefficients mean:
#   compare_to_zero = FALSE -> "Condition"
#       R's default coding. Intercept = mean of the FIRST condition column;
#       every other coefficient = (that condition - first condition).
#   compare_to_zero = TRUE, several conditions -> "0 + Condition"
#       "0 +" removes the shared intercept, so there is one coefficient per
#       condition = that condition's mean, each tested against 0.
#   compare_to_zero = TRUE, one condition -> "1"
#       A factor with one level can't be a term, so the model is intercept-only:
#       the intercept is that condition's mean, tested against 0.
# The p-values in the "Fixed Effects" table test exactly these coefficients
# (unadjusted). Between-condition tests come from the Comparisons table.
fixed_part <- if (!compare_to_zero) {
  "Condition"
} else if (!multi) {
  "1"
} else {
  "0 + Condition"
}

if (!compare_to_zero && !multi) {
  stop("Only one condition column -- nothing to compare. ",
       "Set compare_to_zero <- TRUE to test it against 0.")
}

# 0/1 indicator per condition, used for a separate slice SD per condition.
dummy_names  <- paste0("sc_", seq_along(cond_levels))
dummy_labels <- setNames(paste("Slice x", cond_levels), dummy_names)
for (i in seq_along(cond_levels)) {
  df_long[[dummy_names[i]]] <- as.numeric(df_long$Condition == cond_levels[i])
}

# 6a. lme4 models: one cell-level (residual) variance for every condition.
#   shared        : + (1 | Mouse:Slice:Condition)
#   per_condition : + (0 + sc_k | Mouse:Slice) for each condition k
lmer_formula <- function(sc) {
  sc_part <- switch(
    sc,
    none          = NULL,
    shared        = "(1 | Mouse:Slice:Condition)",
    per_condition = paste0("(0 + ", dummy_names, " | Mouse:Slice)", collapse = " + ")
  )
  as.formula(paste(c(paste0("Response ~ ", fixed_part, " + (1 | ", base_groups, ")"), sc_part),
                   collapse = " + "))
}

# Every lmer fit is labelled, so a warning says WHICH model it came from.
# Several models are fitted to run the tests in section 7; only one is used
# for the results (see "Model used"). bobyqa is a more robust optimizer than
# the default for models with many variance terms.
fit_warnings <- character(0)
fit_lmer <- function(sc) {
  label <- paste0("Model ", c(none = "A", shared = "C", per_condition = "E")[[sc]],
                  " (lmer, slice x condition: ", sc, ")")
  tryCatch(
    withCallingHandlers(
      # do.call() makes summary() print the actual formula.
      do.call("lmer", list(formula = lmer_formula(sc), data = quote(df_long),
                           control = lmerControl(optimizer = "bobyqa"))),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, paste0("[", label, "] ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      message = function(m) {   # "boundary (singular) fit" arrives as a message
        fit_warnings <<- c(fit_warnings, paste0("[", label, "] ", trimws(conditionMessage(m))))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      fit_warnings <<- c(fit_warnings, paste0("[", label, "] FAILED: ", conditionMessage(e)))
      NULL
    }
  )
}

# 6b. nlme models: same fixed and random effects, fitted by REML. The slice x
#     condition term is written as condition effects within each slice:
#     pdIdent = one shared SD, pdDiag = one SD per condition.
nlme_random <- function(sc) {
  slice_part <- switch(
    sc,
    none          = quote(~ 1),
    shared        = quote(nlme::pdBlocked(list(nlme::pdIdent(~ 1), nlme::pdIdent(~ 0 + Condition)))),
    per_condition = quote(nlme::pdBlocked(list(nlme::pdIdent(~ 1), nlme::pdDiag(~ 0 + Condition))))
  )
  if (cells_repeated) {
    bquote(list(Mouse = ~ 1, Slice = .(slice_part), Cell = ~ 1))
  } else {
    bquote(list(Mouse = ~ 1, Slice = .(slice_part)))
  }
}

fixed_formula <- as.formula(paste("Response ~", fixed_part))

fit_lme <- function(sc, variance_per_condition) {
  w <- if (variance_per_condition) quote(nlme::varIdent(form = ~ 1 | Condition)) else NULL
  tryCatch(
    eval(bquote(nlme::lme(fixed = .(fixed_formula), random = .(nlme_random(sc)),
                          weights = .(w), data = df_long, method = "REML",
                          control = nlme::lmeControl(maxIter = 200, msMaxIter = 200,
                                                     niterEM = 50)))),
    error = function(e) {
      warning("nlme::lme did not converge (slice x condition = ", sc, ", cell variance ",
              if (variance_per_condition) "per condition" else "shared", "): ",
              conditionMessage(e))
      NULL
    }
  )
}

# Likelihood-ratio test between two REML fits with identical fixed effects.
# Variance terms tested against 0 sit on a boundary, so those p-values are
# conservative (slightly too large).
lr_test <- function(m0, m1) {
  if (is.null(m0) || is.null(m1)) return(c(stat = NA, p = NA))
  ll0  <- logLik(m0)
  ll1  <- logLik(m1)
  stat <- max(0, 2 * (as.numeric(ll1) - as.numeric(ll0)))
  df   <- attr(ll1, "df") - attr(ll0, "df")
  c(stat = stat, p = pchisq(stat, df, lower.tail = FALSE))
}

fit_row <- function(label, m) {
  ll <- logLik(m)
  data.frame(Model = label, n_par = attr(ll, "df"), AIC = AIC(m),
             BIC = BIC(m), logLik = as.numeric(ll))
}

# 7. Choose the variance structure ---------------------------------------------
sc_labels <- c(none = "No slice x condition", shared = "Slice x condition, shared SD",
               per_condition = "Slice x condition, SD per condition")

# 7a. Slice x condition (fitted with lmer).
sc_needed <- if (!multi) {
  "none"
} else if (check_variance_homogeneity || slice_condition == "auto") {
  c("none", "shared", "per_condition")
} else {
  slice_condition
}
lmer_fits <- setNames(lapply(sc_needed, fit_lmer), sc_needed)

slice_test <- NULL
slice_p    <- NA    # shared vs per-condition slice SD
if (length(sc_needed) == 3) {
  ok <- !vapply(lmer_fits, is.null, TRUE)
  t_shared <- lr_test(lmer_fits$none, lmer_fits$shared)
  t_percon <- lr_test(lmer_fits$shared, lmer_fits$per_condition)
  slice_p  <- t_percon[["p"]]
  slice_test <- do.call(rbind, lapply(names(lmer_fits)[ok], function(nm) fit_row(sc_labels[[nm]], lmer_fits[[nm]])))
  lr_vals <- list(none = c(stat = NA, p = NA), shared = t_shared, per_condition = t_percon)[ok]
  slice_test$LR_stat <- vapply(lr_vals, `[[`, 0, "stat")
  slice_test$p       <- vapply(lr_vals, `[[`, 0, "p")
  rownames(slice_test) <- NULL
}

slice_condition_used <- if (!multi) {
  "none"
} else if (slice_condition != "auto") {
  slice_condition
} else if (!is.na(slice_p) && slice_p < variance_alpha) {
  "per_condition"
} else {
  "shared"
}

model_lmer <- lmer_fits[[slice_condition_used]]
if (is.null(model_lmer)) {
  stop("The lmer model with slice_condition = \"", slice_condition_used,
       "\" could not be fitted -- see the warning above.")
}

# 7b. Cell-level (residual) variance: shared vs per condition (fitted with nlme,
#     using the slice x condition structure chosen above).
m_equal <- m_hetero <- variance_test <- NULL
variance_p <- NA

if (multi && (check_variance_homogeneity || model_engine != "lmer")) {
  m_equal  <- fit_lme(slice_condition_used, FALSE)
  m_hetero <- fit_lme(slice_condition_used, TRUE)
}

if (!is.null(m_equal) && !is.null(m_hetero)) {
  t_res <- lr_test(m_equal, m_hetero)
  variance_p <- t_res[["p"]]
  variance_test <- rbind(fit_row("Cell SD shared", m_equal),
                         fit_row("Cell SD per condition", m_hetero))
  variance_test$LR_stat <- c(NA, t_res[["stat"]])
  variance_test$p       <- c(NA, variance_p)
}

# 7d. Final model.
engine_used <- switch(
  model_engine,
  lmer = "lmer",
  nlme = "nlme",
  auto = if (!is.na(variance_p) && variance_p < variance_alpha) "nlme" else "lmer"
)
if (engine_used == "nlme" && !multi) {
  message("Only one condition, so there is only one residual variance -- using lmer.")
  engine_used <- "lmer"
}
if (engine_used == "nlme" && is.null(m_hetero)) {
  stop("The per-condition cell-variance model (nlme) could not be fitted -- see the ",
       "warning above. Set model_engine <- \"lmer\" to continue.")
}

model <- if (engine_used == "nlme") m_hetero else model_lmer

slice_desc <- c(none = "none", shared = "shared SD", per_condition = "SD per condition")[[slice_condition_used]]
cell_desc  <- if (engine_used == "nlme") "SD per condition" else "shared SD"
model_desc <- paste0(
  if (engine_used == "nlme") "nlme::lme" else "lme4::lmer",
  "; random intercepts ", base_groups,
  "; slice x condition: ", slice_desc,
  "; cell-level SD: ", cell_desc
)

# 8. Results tables ------------------------------------------------------------
# All estimates are back-transformed to the ORIGINAL scale, so a log-transformed
# analysis reports geometric means and RATIOS instead of log10 values.
back           <- if (log_transform) function(x) 10^x else function(x) x
estimate_label <- if (log_transform) "Geometric_mean" else "Mean"
effect_label   <- if (log_transform) "Ratio" else "Difference"

# --- Number formatting: used for EVERY table (console and note) and the plot ---
# Tables are kept as full-precision numbers and only formatted when shown, so
# rounding never feeds into a calculation.
#   p-values : p_decimals places; "<0.0001" below 0.0001; a p that would round
#              to 0 (e.g. 0.003 with 2 places) shows as "<0.01".
#   counts   : whole numbers (n, n_mice, n_slice, n_cells, n_par).
#   other    : num_decimals places, trailing zeros kept (1.50, not 1.5).
p_cols_all     <- c("p", "p_adj", "p_value", "p_vs_0", "SW_p", "DAgo_p")
count_cols_all <- c("n", "n_par", "n_mice", "n_slice", "n_cells", "n_observations")

fmt_num <- function(x) ifelse(is.na(x), "", sprintf(paste0("%.", num_decimals, "f"), x))

fmt_p <- function(p) {
  out <- sprintf(paste0("%.", p_decimals, "f"), p)
  floor_txt <- paste0("<", sprintf(paste0("%.", p_decimals, "f"), 10^-p_decimals))
  out[!is.na(p) & round(p, p_decimals) == 0] <- floor_txt
  out[!is.na(p) & p < 0.0001] <- "<0.0001"
  out[is.na(p)] <- ""
  out
}

fmt_tbl <- function(tbl) {
  tbl <- as.data.frame(tbl)   # plain data.frame prints cleaner than a tibble
  for (col in names(tbl)) {
    x <- tbl[[col]]
    if (!is.numeric(x)) next
    tbl[[col]] <- if (col %in% p_cols_all) fmt_p(x)
    else if (col %in% count_cols_all) ifelse(is.na(x), "", as.character(round(x)))
    else fmt_num(x)
  }
  tbl
}

show_table <- function(tbl, title) {
  cat("\n", title, "\n", strrep("=", nchar(title)), "\n", sep = "")
  print(fmt_tbl(tbl), row.names = FALSE, right = TRUE)
}

# lmer and lme store their results differently; these return the same layout.
fixed_effects_of <- function(m) {
  if (inherits(m, "lme")) {
    tt <- summary(m)$tTable
    data.frame(Effect = rownames(tt), Estimate = tt[, "Value"], Std_Error = tt[, "Std.Error"],
               DF = tt[, "DF"], t_value = tt[, "t-value"], p_value = tt[, "p-value"])
  } else {
    tt <- coef(summary(m))
    data.frame(Effect = rownames(tt), Estimate = tt[, "Estimate"], Std_Error = tt[, "Std. Error"],
               DF = tt[, "df"], t_value = tt[, "t value"], p_value = tt[, "Pr(>|t|)"])
  }
}

random_effects_of <- function(m) {
  if (inherits(m, "lme")) {
    vc  <- nlme::VarCorr(m)
    rn  <- rownames(vc)
    grp <- "Group"
    rows <- list()
    for (i in seq_along(rn)) {
      if (grepl(" =$", rn[i])) { grp <- sub(" =$", "", rn[i]); next }   # "Slice =" starts a level
      label <- if (rn[i] == "Residual") {
        if (!is.null(m$modelStruct$varStruct)) "Residual / cell (first condition; see Variance check)" else "Residual / cell"
      } else if (rn[i] == "(Intercept)") {
        grp
      } else {
        paste(grp, "x", sub("^Condition", "", rn[i]))
      }
      rows[[length(rows) + 1]] <- data.frame(Group = label,
                                             Variance = suppressWarnings(as.numeric(vc[i, "Variance"])),
                                             SD = suppressWarnings(as.numeric(vc[i, "StdDev"])))
    }
    out <- do.call(rbind, rows)
    # pdIdent repeats one shared SD on every condition row; show it once.
    if (slice_condition_used == "shared") {
      is_sc <- grepl("^Slice x ", out$Group)
      if (any(is_sc)) {
        first <- which(is_sc)[1]
        out$Group[first] <- "Slice x Condition"
        out <- out[!is_sc | seq_len(nrow(out)) == first, ]
      }
    }
    out
  } else {
    vc <- as.data.frame(VarCorr(m))
    label <- ifelse(vc$var1 %in% dummy_names, dummy_labels[vc$var1],
                    ifelse(vc$grp == "Mouse:Slice:Condition", "Slice x Condition",
                           ifelse(vc$grp == "Residual", "Residual / cell", sub(":.*$", "", vc$grp))))
    data.frame(Group = unname(label), Variance = vc$vcov, SD = vc$sdcor)
  }
}

random_tbl <- random_effects_of(model)
rownames(random_tbl) <- NULL

# Variance tables for decisions 3 and 4 ----------------------------------------
# For each decision the model is fitted twice, and the SDs from BOTH versions
# are shown side by side:
#   assuming equal   : one SD for all conditions
#   allowing unequal : a separate SD per condition
# Both columns are model estimates of the same thing, so they can be compared
# directly. Which version is used is decided by the likelihood-ratio test (does
# the unequal version fit the data clearly better?), not by which SDs are smaller.
slice_sd_tbl <- cell_sd_tbl <- NULL
sigma_of <- function(m) if (inherits(m, "lme")) m$sigma else sigma(m)
cell_sds_of <- function(m) {
  if (is.null(m)) return(rep(NA_real_, length(cond_levels)))
  if (inherits(m, "lme") && !is.null(m$modelStruct$varStruct)) {
    r <- coef(m$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)
    unname(m$sigma * r[cond_levels])
  } else {
    rep(sigma_of(m), length(cond_levels))
  }
}
count_by_cond <- function(x) {
  as.numeric(tapply(x, df_long$Condition, function(v) length(unique(v)))[cond_levels])
}
# Slice SD in a condition = how much slices differ from each other in that
# condition = sqrt(shared slice variance + that condition's slice x condition
# variance). The two parts trade off in the fit, so only the total is shown.
total_slice_sd <- function(m, cnd) {
  if (is.null(m)) return(NA_real_)
  re <- random_effects_of(m)
  v  <- function(g) { x <- re$Variance[re$Group == g]; if (length(x) == 0 || is.na(x[1])) 0 else x[1] }
  sqrt(v("Slice") + v(paste("Slice x", cnd)) + v("Slice x Condition"))
}

if (check_variance_homogeneity && multi) {
  slice_sd_tbl <- data.frame(
    Condition        = cond_levels,
    n_slice          = count_by_cond(paste(df_long$Mouse, df_long$Slice)),
    SD_none          = vapply(cond_levels, function(c) total_slice_sd(lmer_fits$none, c), 0, USE.NAMES = FALSE),
    SD_shared        = vapply(cond_levels, function(c) total_slice_sd(lmer_fits$shared, c), 0, USE.NAMES = FALSE),
    SD_per_condition = vapply(cond_levels, function(c) total_slice_sd(lmer_fits$per_condition, c), 0, USE.NAMES = FALSE)
  )
  cell_sd_tbl <- data.frame(
    Condition           = cond_levels,
    n_cells             = as.numeric(table(df_long$Condition)[cond_levels]),
    SD_assuming_equal   = cell_sds_of(if (!is.null(m_equal)) m_equal else model_lmer),
    SD_allowing_unequal = cell_sds_of(m_hetero)
  )
}

# Comparisons BETWEEN conditions (later - earlier, so signs match the model
# coefficients). Which comparisons and which correction: see SETTINGS.
if (multi) {
  auto_adjust <- c(all_pairs = "tukey", vs_first = "dunnettx",
                   consecutive = "holm", custom = "holm")[[comparison_type]]
  adjust_used <- if (p_adjust == "auto") auto_adjust else p_adjust
}

make_contrasts <- function(e) {
  lv <- cond_levels
  switch(
    comparison_type,
    all_pairs   = pairs(e, reverse = TRUE, adjust = adjust_used),
    vs_first    = contrast(e, "trt.vs.ctrl", ref = if (exists("labstats_control", inherits = TRUE)) match(labstats_control, lv) else 1, adjust = adjust_used),
    consecutive = contrast(e, "consec", adjust = adjust_used),
    custom      = {
      if (length(custom_pairs) == 0) {
        stop("comparison_type = 'custom' but custom_pairs is empty.")
      }
      weights <- lapply(custom_pairs, function(pr) {
        if (length(pr) != 2 || !all(pr %in% lv)) {
          stop("Each custom_pairs entry must be c(A, B) using these names: ",
               paste(lv, collapse = ", "))
        }
        w <- setNames(rep(0, length(lv)), lv)
        w[pr[1]] <- 1
        w[pr[2]] <- -1
        w
      })
      names(weights) <- vapply(custom_pairs, function(pr) paste(pr[1], "-", pr[2]), "")
      contrast(e, method = weights, adjust = adjust_used)
    }
  )
}

# 8a. Estimated mean (or geometric mean) per condition, 95% CI
# Degrees of freedom (decision 6). Satterthwaite df always come from the lmer
# model with the chosen slice x condition structure:
#   * lmer engine: these are simply lmer's own df.
#   * nlme engine: nlme counts df at the CELL level (e.g. 133) even though only
#     ~8-9 slices replicate the condition effect, which makes p-values too
#     small. So nlme's estimates and SEs are combined with lmer's df.
# df_method = "per_comparison": each mean / comparison gets its own df.
# df_method = "smallest"      : every one gets the smallest df (cautious).
df_note <- NULL
if (multi) {
  emm_lmer     <- emmeans(model_lmer, ~ Condition, level = 0.95)
  df_means_vec <- as.data.frame(summary(emm_lmer))$df
  df_diffs_vec <- as.data.frame(make_contrasts(emm_lmer))$df
  df_vs_ref    <- as.data.frame(contrast(emm_lmer, "trt.vs.ctrl", ref = 1, adjust = "none"))$df
  if (df_method == "smallest") {
    df_means_vec[] <- min(df_means_vec)
    df_diffs_vec[] <- min(df_diffs_vec)
    df_vs_ref[]    <- min(df_diffs_vec)
  }
  df_means <- min(df_means_vec)
  df_diffs <- min(df_diffs_vec)
  
  # Estimates and SEs from the model actually used.
  emm_base <- if (engine_used == "nlme") emmeans(model, ~ Condition, level = 0.95) else emm_lmer
  
  # Each row is computed with its own df (update() sets one df for the grid,
  # so it is done row by row).
  emm_df <- do.call(rbind, lapply(seq_along(df_means_vec), function(i) {
    as.data.frame(summary(update(emm_base, df = df_means_vec[i]), infer = c(TRUE, TRUE)))[i, ]
  }))
  
  df_range <- function(x) {
    if (diff(range(x)) < 0.005) fmt_num(x[1]) else paste0(fmt_num(min(x)), " to ", fmt_num(max(x)))
  }
  df_note <- paste0(
    if (engine_used == "nlme") "nlme estimates and SEs with df from the matching lmer model (Satterthwaite). "
    else "Satterthwaite df from lmer. ",
    if (df_method == "per_comparison") "Each mean and comparison has its own df"
    else "Every mean and comparison uses the smallest df",
    ": means df = ", df_range(df_means_vec), ", comparisons df = ", df_range(df_diffs_vec), "."
  )
} else {
  # Single condition: intercept-only model, t-based CI built by hand.
  fe     <- fixed_effects_of(model)
  crit   <- qt(0.975, fe$DF[1])
  emm_df <- data.frame(
    Condition = cond_levels,
    emmean    = fe$Estimate[1],
    lower.CL  = fe$Estimate[1] - crit * fe$Std_Error[1],
    upper.CL  = fe$Estimate[1] + crit * fe$Std_Error[1],
    SE        = fe$Std_Error[1],
    df        = fe$DF[1],
    p.value   = fe$p_value[1]
  )
}

# Each condition's own mean. Its CI includes ALL the variation between mice and
# slices, because it describes the level of one condition. It is NOT used for
# the comparisons: there each slice is its own control, so the slice and mouse
# offsets cancel (see the Comparisons table/plot).
estimates_tbl <- data.frame(
  Condition = as.character(emm_df$Condition),
  Estimate  = back(emm_df$emmean),
  SEM       = emm_df$SE,
  CI_lower  = back(emm_df$lower.CL),
  CI_upper  = back(emm_df$upper.CL),
  DF        = emm_df$df
)
names(estimates_tbl)[2] <- estimate_label
# SEM = the model's standard error of each condition mean (same source as the
# CI). With log_transform it is on the log10 scale, so it is labelled that way.
if (log_transform) names(estimates_tbl)[names(estimates_tbl) == "SEM"] <- "SEM_log10"

# "p vs 0" only means something for un-logged data tested against 0.
if (compare_to_zero && !log_transform) estimates_tbl$p_vs_0 <- emm_df$p.value

# Fixed-effects table, with the same df as above. Rows are either a condition
# mean (intercept, or every row when compare_to_zero) or a condition minus the
# first condition.
fixed_tbl <- fixed_effects_of(model)
rownames(fixed_tbl) <- NULL
if (multi) {
  fixed_tbl$DF      <- if (compare_to_zero) df_means_vec else c(df_means_vec[1], df_vs_ref)
  fixed_tbl$p_value <- 2 * pt(-abs(fixed_tbl$t_value), fixed_tbl$DF)
}

# 8b. Comparisons table (built from make_contrasts above).
comparisons_tbl    <- NULL
comparison_pairs   <- NULL     # (A, B, p) for each comparison, used by the plot
adjust_label       <- NULL

if (multi) {
  lv <- cond_levels
  ct <- make_contrasts(emm_base)   # used for the labels below
  
  # Row i is computed with comparison i's own df. The correction still covers
  # the whole family (Tukey/Dunnett use the family size; Holm/Bonferroni are
  # applied across the raw p-values afterwards).
  ct_rows <- lapply(seq_along(df_diffs_vec), function(i) {
    ct_i <- make_contrasts(update(emm_base, df = df_diffs_vec[i]))
    list(ci    = as.data.frame(confint(ct_i))[i, ],
         p_fam = as.data.frame(ct_i)$p.value[i],
         p_raw = as.data.frame(summary(ct_i, adjust = "none"))$p.value[i])
  })
  ct_ci <- do.call(rbind, lapply(ct_rows, `[[`, "ci"))
  ct_p  <- if (adjust_used %in% c("tukey", "dunnettx")) {
    vapply(ct_rows, `[[`, 0, "p_fam")
  } else {
    p.adjust(vapply(ct_rows, `[[`, 0, "p_raw"), method = adjust_used)
  }
  
  # Which two conditions each comparison involves: the +1 and -1 coefficients.
  cf <- as.data.frame(coef(ct))
  cf <- cf[, grepl("^c\\.", names(cf)), drop = FALSE]
  pair_A <- vapply(cf, function(cc) lv[which(cc > 0)[1]], "")
  pair_B <- vapply(cf, function(cc) lv[which(cc < 0)[1]], "")
  
  comparisons_tbl <- data.frame(
    Comparison = paste(pair_A, if (log_transform) "/" else "-", pair_B),
    Effect     = back(ct_ci$estimate),
    CI_lower   = back(ct_ci$lower.CL),
    CI_upper   = back(ct_ci$upper.CL),
    DF         = ct_ci$df,
    p_adj      = ct_p
  )
  names(comparisons_tbl)[2] <- effect_label
  
  comparison_pairs <- data.frame(A = pair_A, B = pair_B, p = ct_p)
  
  adjust_name <- switch(
    adjust_used,
    tukey = "Tukey", dunnettx = "Dunnett", holm = "Holm", bonferroni = "Bonferroni",
    none = "no correction"
  )
  # emmeans has no Holm-type interval, so intervals use the (more cautious)
  # Bonferroni adjustment whenever p-values use Holm.
  ci_note <- if (adjust_used == "holm") "; CI Bonferroni-adjusted" else ""
  adjust_label <- paste0(adjust_name, ci_note)
}

# 8b-2. Where the degrees of freedom come from ---------------------------------
# Satterthwaite df, computed exactly as lmerTest does it:
#   df = 2 * SE^4 / Var(SE^2)
#   SE^2      = L' V L  -- squared standard error of the comparison
#               (L = the comparison's weights on the model coefficients,
#                V = covariance matrix of the coefficient estimates)
#   Var(SE^2) = g' A g  -- how uncertain SE^2 is, because it is built from
#               ESTIMATED variance parameters (mouse, slice, cell SDs)
#               g = how much SE^2 changes when each variance parameter changes
#               A = covariance matrix of the variance-parameter estimates
# Var(SE^2) is split by parameter (g_k * (A g)_k), so the table shows which
# estimated SD limits the df. Uses the lmer model, which is the df source.
df_source_label <- function(par) {
  vapply(par, function(nm) {
    if (nm == "sigma") return("Cell noise")
    grp  <- sub("\\.[^.]*$", "", nm)
    term <- sub("^.*\\.", "", nm)
    if (term %in% dummy_names) return(dummy_labels[[term]])
    if (grp == "Mouse") "Mouse"
    else if (startsWith(grp, "Cell")) "Cell offset"
    else if (grp == "Mouse:Slice:Condition") "Slice x Condition"
    else if (startsWith(grp, "Slice")) "Slice"
    else grp
  }, "", USE.NAMES = FALSE)
}
df_explain <- paste0(
  "Satterthwaite df = 2 x SE^4^ / Var(SE^2^). The SE of each comparison is built from the estimated SDs ",
  "(mouse, slice, cell). Each of those SDs is itself uncertain -- an SD estimated from a few mice is much less ",
  "certain than one estimated from 100 cells -- and Var(SE^2^) is how much that uncertainty carries into the SE. ",
  "The % columns split Var(SE^2^) by source: the source with the largest share is what limits the df. If one source ",
  "were the only one, the df would equal the df of that SD (roughly its number of units, minus what the model uses up). ",
  "Shares can be slightly negative or above 100 because the SD estimates are correlated. ",
  "Computed from the lmer version of the model (the df source)",
  if (df_method == "smallest") "; the Comparisons table uses the smallest of these DF (df_method = \"smallest\")."
  else "; DF here matches the Comparisons table."
)

df_breakdown_tbl <- NULL
if (multi) {
  df_breakdown_tbl <- tryCatch({
    m   <- model_lmer
    L   <- make_contrasts(emm_lmer)@linfct
    Vb  <- as.matrix(m@vcov_beta)
    A   <- as.matrix(m@vcov_varpar)
    J   <- m@Jac_list
    src <- df_source_label(c(names(getME(m, "theta")), "sigma"))
    lev <- unique(src)
    do.call(rbind, lapply(seq_len(nrow(L)), function(i) {
      l     <- L[i, ]
      se2   <- sum(l * (Vb %*% l))
      g     <- vapply(J, function(Jk) sum(l * (as.matrix(Jk) %*% l)), 0)
      Ag    <- as.vector(A %*% g)
      vse2  <- sum(g * Ag)
      share <- tapply(100 * g * Ag / vse2, factor(src, levels = lev), sum)
      row <- data.frame(Comparison = comparisons_tbl$Comparison[i], DF = 2 * se2^2 / vse2,
                        check.names = FALSE)
      for (k in seq_along(lev)) row[[paste0("% ", lev[k])]] <- unname(share[k])
      row
    }))
  }, error = function(e) {
    fit_warnings <<- c(fit_warnings, paste0("[df breakdown] could not be computed: ", conditionMessage(e)))
    NULL
  })
}

# 8d. Analysis decisions -- one row per choice, in plain words -----------------
# This table goes at the top of the note and the console output, so every
# result can be read knowing exactly which model produced it.
how_chosen <- function(setting_value, p, alpha_word) {
  if (setting_value == "auto") {
    if (is.na(p)) "Auto: test could not be run, used the simpler option"
    else paste0("Auto: likelihood-ratio test p = ", fmt_p(p),
                if (p < variance_alpha) paste0(" < ", variance_alpha, ", so ", alpha_word)
                else paste0(" >= ", variance_alpha, ", so kept the simpler option"))
  } else {
    "Fixed by you in SETTINGS"
  }
}

n_slices <- n_distinct(df_long$Mouse, df_long$Slice)
n_mice   <- n_distinct(df_long$Mouse)
min_df   <- if (multi) min(ct_ci$df) else NA

decisions_tbl <- data.frame(
  Decision = c(
    "1. Data scale",
    "2. Nesting (random intercepts)",
    "3. Slice x condition",
    "4. Cell-level noise",
    "5. Software (lmer vs nlme)",
    "6. Degrees of freedom",
    "7. Comparisons and correction"
  ),
  Choice = c(
    if (log_transform) "log10" else "Raw values",
    base_groups,
    c(none = "Not included", shared = "Included, one SD for all conditions",
      per_condition = "Included, separate SD per condition")[[slice_condition_used]],
    if (engine_used == "nlme") "Separate SD per condition" else "One SD for all conditions",
    if (engine_used == "nlme") "nlme::lme" else "lme4::lmer",
    if (!multi) "Model df"
    else paste0(if (df_method == "per_comparison") "Own df per comparison" else "Smallest df for all",
                " (Satterthwaite", if (engine_used == "nlme") ", from lmer" else "", "): ",
                df_range(df_diffs_vec)),
    if (multi) paste0(gsub("_", " ", comparison_type), ", ", adjust_name) else "None (one condition)"
  ),
  Why = c(
    src_log,
    if (cells_repeated) "Cells were measured in more than one condition"
    else "Each cell measured once, so Cell is not estimable",
    if (!multi) "Only one condition"
    else if (slice_condition == "auto") how_chosen("auto", slice_p, "used a separate SD per condition")
    else "Fixed by you in SETTINGS",
    if (!multi) "Only one condition"
    else how_chosen(if (model_engine == "auto") "auto" else "fixed", variance_p, "used a separate SD per condition"),
    "Follows from decision 4: lmer cannot give each condition its own cell-level SD",
    paste0(if (df_method == "per_comparison") "Field standard (lmerTest/emmeans default). "
           else "Fixed by you: cautious option. ",
           if (engine_used == "nlme") "Taken from lmer because nlme counts df per cell, which overstates replication" else ""),
    src_cmp
  ),
  Meaning = c(
    if (log_transform) "Effects are ratios (fold changes); means are geometric means"
    else "Effects are differences in the original units",
    "Cells are grouped within slices within mice; cells from one slice are not independent",
    if (slice_condition_used == "none") "Assumes every slice responds to the drug by the same amount"
    else paste0("Slices may respond differently, so the drug effect is judged against slice-to-slice ",
                "variation (", n_slices, " slices, ", n_mice, " mice), not against cells"),
    if (engine_used == "nlme") "Noisier conditions get less weight instead of distorting the others"
    else "Assumes cells are equally noisy in every condition",
    "Software only -- the model itself is described by decisions 2-4",
    if (!multi) "--"
    else if (min_df < 30) "Low df: the comparisons are limited by the number of slices, as they should be (see 'Where the degrees of freedom come from')"
    else "High df: the comparisons are limited by the number of cells",
    if (multi) "p-values in the Comparisons table are adjusted for the number of comparisons" else "--"
  )
)

decision_summary <- paste0(
  "Every cell measurement (", nrow(df_long), ") is used, grouped by ",
  if (cells_repeated) "cell, slice and mouse. " else "slice and mouse. ",
  if (slice_condition_used == "none") {
    "Slices are assumed to respond to each condition by the same amount. "
  } else {
    paste0("Each slice is allowed its own response to each condition",
           if (slice_condition_used == "per_condition") ", with that slice-to-slice spread differing by condition" else "",
           ", so the condition comparisons rest on the ", n_slices, " slices from ", n_mice, " mice. ")
  },
  if (engine_used == "nlme") {
    "Cell-to-cell noise is estimated separately for each condition (nlme)."
  } else {
    "Cell-to-cell noise is assumed equal across conditions (lmer)."
  }
)

# 8d-2. The six possible models -------------------------------------------------
# Decisions 3 and 4 are two separate choices (3 x 2), so there are six possible
# models. All six share the fixed effects, the Mouse/Slice/Cell nesting and the
# comparisons; they differ only in how the variability is modelled.
models_tbl <- data.frame(
  Model = c("A", "B", "C", "D", "E", "F"),
  Slice_x_condition = rep(c("None", "Shared SD", "SD per condition"), each = 2),
  Cell_noise = rep(c("One SD", "SD per condition"), times = 3),
  Software = rep(c("lmer", "nlme"), times = 3),
  Assumes = c(
    "Every slice responds by the same amount; cells equally noisy in every condition",
    "Every slice responds by the same amount; cell noise differs by condition",
    "Slices respond differently, by a similar spread in every condition; cells equally noisy",
    "Slices respond differently, by a similar spread in every condition; cell noise differs",
    "Slices respond differently, and that spread differs by condition; cells equally noisy",
    "Slices respond differently, and that spread differs by condition; cell noise differs"
  )
)
used_row <- match(slice_condition_used, c("none", "shared", "per_condition")) * 2 -
  (engine_used == "lmer")
models_tbl$Used <- ifelse(seq_len(nrow(models_tbl)) == used_row, "<- used", "")
model_letter <- models_tbl$Model[used_row]

# 8d-3. Model-fitting checks --------------------------------------------------
# Several models are fitted only to run the tests (section 7). For each, list
# the variance terms estimated at 0 ("boundary (singular) fit") in words, and
# keep any other warning (e.g. convergence) as is.
zero_terms <- function(m) {
  if (is.null(m)) return(character(0))
  re <- random_effects_of(m)
  re <- re[!grepl("^Residual", re$Group), , drop = FALSE]
  re$Group[!is.na(re$SD) & re$SD < 1e-3 * max(c(re$SD, sigma(m)), na.rm = TRUE)]
}
lmer_letter <- c(none = "A", shared = "C", per_condition = "E")
df_source   <- lmer_letter[[slice_condition_used]]

zero_lines <- c(
  vapply(names(lmer_fits), function(nm) {
    z <- zero_terms(lmer_fits[[nm]])
    if (length(z) == 0) return(NA_character_)
    paste0("Model ", lmer_letter[[nm]], " (lmer", if (lmer_letter[[nm]] == df_source && multi) ", source of the df" else "",
           "): ", paste(z, collapse = ", "))
  }, ""),
  if (engine_used == "nlme") {
    z <- zero_terms(model)
    if (length(z) > 0) paste0("Model ", model_letter, " (nlme, USED for the results): ", paste(z, collapse = ", "))
  }
)
zero_lines <- zero_lines[!is.na(zero_lines)]
other_warnings <- fit_warnings[!grepl("boundary \\(singular\\)", fit_warnings)]

fit_check_lines <- c(
  paste0("Results come from model ", model_letter, ". The other models are fitted only for the tests",
         if (engine_used == "nlme") paste0(" (and model ", df_source, " supplies the degrees of freedom)") else "", "."),
  if (length(zero_lines) > 0) c(
    "Variance terms estimated at 0 -- that level adds no detectable extra variation. Common with few mice or slices, and harmless for the comparisons:",
    paste0("  - ", zero_lines)
  ) else "No variance terms were estimated at 0.",
  if (length(other_warnings) > 0) c(
    "Warnings to check (only matter if they are on the model used or the df source):",
    paste0("  - ", other_warnings)
  ) else "No convergence warnings or fitting errors."
)

# 8e. Plain-language reading of the variance checks ---------------------------
lr_glossary <- paste0(
  "How to read these tables: logLik = how well the model fits the data (higher, i.e. less negative, ",
  "is better). Extra parameters always raise it a little, so AIC and BIC add a penalty per ",
  "parameter (n_par; BIC's penalty is larger). Lower AIC/BIC = better fit for its complexity; a drop ",
  "of more than ~10 is strong evidence. LR_stat = 2 x the gain in logLik over the row above; ",
  "p tests whether that gain is bigger than the extra parameters would give by chance."
)

pv <- function(p) {
  if (is.na(p)) return("test not available")
  txt <- fmt_p(p)
  if (startsWith(txt, "<")) paste0("p ", sub("<", "< ", txt, fixed = TRUE)) else paste0("p = ", txt)
}
describe_sds <- function(tbl, eq_col, uneq_col, eq_name = "Assuming equal", uneq_name = "Allowing unequal") {
  eq <- tbl[[eq_col]]; un <- tbl[[uneq_col]]
  txt <- paste0(eq_name, ": ", fmt_num(eq[1]), " in every condition.")
  if (!all(is.na(un))) {
    hi <- which.max(un); lo <- which.min(un)
    txt <- paste0(txt, " ", uneq_name, ": from ", fmt_num(un[lo]), " (", tbl$Condition[lo], ") to ",
                  fmt_num(un[hi]), " (", tbl$Condition[hi], ")",
                  if (un[lo] >= 0.005) paste0(", the largest ", fmt_num(un[hi] / un[lo]), " times the smallest.")
                  else paste0("; ", tbl$Condition[lo], " is essentially 0."))
  }
  txt
}
test_reason <- function(setting_is_auto, p) {
  if (!setting_is_auto) "fixed by you in SETTINGS."
  else if (is.na(p)) "the test could not be run, so the simpler version was kept."
  else if (p < variance_alpha) paste0("the test gave ", pv(p), " (< ", variance_alpha, "), so the unequal version fits the data clearly better.")
  else paste0("the test gave ", pv(p), " (>= ", variance_alpha, "), so there is no clear evidence the SDs differ and the simpler version was kept.")
}
sd_note <- paste0(
  "A larger SD in the 'allowing unequal' column does not mean a worse model: both columns are estimates of the same ",
  "thing under different assumptions. Which version is used is decided by the test (does it fit the data clearly ",
  "better?), not by the size of the SDs. With many cells the test can detect fairly modest differences; allowing ",
  "unequal SDs costs little when they are in fact similar."
)

cell_reading <- cell_choice <- NULL
if (!is.null(cell_sd_tbl)) {
  cell_reading <- describe_sds(cell_sd_tbl, "SD_assuming_equal", "SD_allowing_unequal")
  cell_choice  <- paste0("**Chosen: ", if (engine_used == "nlme") "a separate cell SD per condition (nlme)"
                         else "one cell SD for all conditions (lmer)", "** -- ",
                         test_reason(model_engine == "auto", variance_p))
}

slice_reading <- slice_choice <- NULL
if (!is.null(slice_sd_tbl)) {
  slice_reading <- paste0("If slices all respond alike (no slice x condition term): ", fmt_num(slice_sd_tbl$SD_none[1]),
                          " in every condition. ",
                          describe_sds(slice_sd_tbl, "SD_shared", "SD_per_condition",
                                       "One shared spread", "A spread per condition"))
  p_any <- if (!is.null(slice_test) && nrow(slice_test) == 3) slice_test$p[2] else NA
  slice_choice <- paste0(
    "**Chosen: ", c(none = "no slice x condition term", shared = "slice x condition with one shared spread",
                    per_condition = "slice x condition with a spread per condition")[[slice_condition_used]], "** -- ",
    if (slice_condition != "auto") "fixed by you in SETTINGS."
    else paste0("in \"auto\" the slice x condition term is always included (row 2 vs 1: ", pv(p_any), "); ",
                if (!is.na(slice_p) && slice_p < variance_alpha)
                  paste0("a spread per condition fits the data clearly better (row 3 vs 2: ", pv(slice_p), ").")
                else paste0("there is no clear evidence the spread differs between conditions (row 3 vs 2: ",
                            pv(slice_p), "), so one shared spread was kept."))
  )
}

# 8c. Descriptive statistics from the raw data (original scale)
# Slices are counted as Mouse + Slice, so "Slice 1" in two mice counts as two.
descriptives_tbl <- df_long |>
  group_by(Condition) |>
  summarise(
    n_mice  = n_distinct(Mouse),
    n_slice = n_distinct(Mouse, Slice),
    n_cells = n(),
    Mean    = mean(Value),
    SD      = sd(Value),
    # Cell-level 95% CI and SEM (every cell treated as independent -- descriptive only)
    CI_lower = mean(Value) - qt(0.975, n() - 1) * sd(Value) / sqrt(n()),
    CI_upper = mean(Value) + qt(0.975, n() - 1) * sd(Value) / sqrt(n()),
    SEM     = sd(Value) / sqrt(n()),
    Median  = median(Value),
    Q1      = quantile(Value, 0.25),
    Q3      = quantile(Value, 0.75),
    .groups = "drop"
  ) |>
  mutate(Condition = as.character(Condition)) |>
  as.data.frame()


# 8f. Writing template for a paper ---------------------------------------------
# Draft Methods / Results / figure-legend text filled in from THIS run. Words in
# [square brackets] are for you to fill in or check. Edit before use.
fmt_ci <- function(lo, hi) paste0("95% CI ", fmt_num(lo), " to ", fmt_num(hi))
p_text <- function(p) {
  txt <- fmt_p(p)
  if (startsWith(txt, "<")) paste0("p ", sub("<", "< ", txt, fixed = TRUE)) else paste0("p = ", txt)
}
n_cells_unique <- n_distinct(df_long$Mouse, df_long$Slice, df_long$Cell)
# How many conditions each cell (Mouse + Slice + Cell) was measured in.
conds_per_cell <- df_long |>
  group_by(Mouse, Slice, Cell) |>
  summarise(k = n_distinct(Condition), .groups = "drop") |>
  pull(k)
n_cells_2plus <- sum(conds_per_cell >= 2)
n_cells_all   <- sum(conds_per_cell == length(cond_levels))
pkg_v <- function(pk) as.character(utils::packageVersion(pk))

methods_text <- paste0(
  "Statistical analysis. ", y_label, " [describe the measurement and units] was analysed with linear ",
  "mixed-effects models to account for the nested structure of the data (", nrow(df_long),
  " total measurements [", paste0(descriptives_tbl$n_cells, " ", descriptives_tbl$Condition, collapse = ", "),
  "] from ", n_cells_unique, " different cells",
  if (length(cond_levels) > 2) paste0(" (", n_cells_2plus, " measured in at least 2 conditions, ",
                                      n_cells_all, " in all ", length(cond_levels), ")")
  else if (multi) paste0(" (", n_cells_all, " measured in both conditions)") else "",
  ", in ", n_slices, " slices from ", n_mice, " mice). ",
  if (log_transform) "Values were log10-transformed before analysis; effects are reported as ratios of geometric means. " else "",
  rout_methods,
  "Condition was a fixed effect, with random intercepts for mouse, slice within mouse",
  if (cells_repeated) ", and cell within slice" else "", ". ",
  if (slice_condition_used != "none") paste0(
    "A random slice-by-condition effect allowed each slice to respond differently to each condition",
    if (slice_condition_used == "per_condition") ", with a separate variance for each condition" else "", ". ") else "",
  if (engine_used == "nlme") "The residual (cell-level) variance was estimated separately for each condition. " else "",
  "The variance structure was chosen with a rule set before the comparisons were examined: ",
  "restricted maximum-likelihood (REML) likelihood-ratio tests (alpha = ", variance_alpha, ") of whether slices respond ",
  "differently between conditions and whether the variances differ between conditions. ",
  "Models were fitted in R (", R.version$major, ".", R.version$minor, ") with ",
  if (engine_used == "nlme") paste0("nlme (", pkg_v("nlme"), "); degrees of freedom were estimated with the Satterthwaite ",
                                    "approximation from the equivalent model in lme4 (", pkg_v("lme4"), ") and lmerTest (", pkg_v("lmerTest"), ")")
  else paste0("lme4 (", pkg_v("lme4"), ") and lmerTest (", pkg_v("lmerTest"), "), using Satterthwaite degrees of freedom"),
  ". ",
  if (multi) paste0("Pairwise comparisons between conditions were computed with emmeans (", pkg_v("emmeans"), ") and adjusted ",
                    "for multiple comparisons with ", if (adjust_used == "none") "no correction" else paste0(adjust_name, "'s method"), ". ") else "",
  "Data are presented as model-estimated ", ifelse(log_transform, "geometric means", "means"),
  " with 95% confidence intervals, and differences with ", if (multi && adjust_used != "none") "multiplicity-adjusted " else "",
  "95% confidence intervals. Statistical significance was set at ", if (multi && adjust_used != "none") "adjusted " else "", "p < 0.05."
)

n_text <- function(cond) {
  r <- descriptives_tbl[descriptives_tbl$Condition == cond, ]
  paste0("n = ", r$n_cells, " cells, ", r$n_slice, " slices, ", r$n_mice, " mice")
}

results_text <- character(0)
if (multi) {
  results_text <- vapply(seq_len(nrow(comparisons_tbl)), function(i) {
    r   <- comparisons_tbl[i, ]
    est <- r[[effect_label]]
    sig <- !is.na(r$p_adj) && r$p_adj < 0.05
    up  <- if (log_transform) est > 1 else est > 0
    pair <- c(comparison_pairs$A[i], comparison_pairs$B[i])
    paste0(
      y_label, " in ", pair[1],
      if (sig) paste0(" was significantly ", if (up) "higher" else "lower", " than in ", pair[2])
      else paste0(" did not differ significantly from ", pair[2]),
      " (", tolower(effect_label), " ", fmt_num(est), ", ", fmt_ci(r$CI_lower, r$CI_upper),
      ", df = ", fmt_num(r$DF), ", ", p_text(r$p_adj), ", ", adjust_name, "-adjusted).",
      if (!sig) " [This is not evidence of no effect: read the CI for the range of effects compatible with the data.]" else ""
    )
  }, "")
}
means_text <- paste0(
  "Model-estimated ", ifelse(log_transform, "geometric mean ", "mean "), y_label, " was ",
  paste0(vapply(seq_len(nrow(estimates_tbl)), function(i) {
    r <- estimates_tbl[i, ]
    paste0(fmt_num(r[[estimate_label]]), " (", fmt_ci(r$CI_lower, r$CI_upper), ") in ", r$Condition,
           " (", n_text(r$Condition), ")")
  }, ""), collapse = "; "), "."
)

legend_text <- paste0(
  "[Figure X.] ", y_label, " in each condition. Points show individual cells, coloured by mouse",
  if (multi) "; lines connect the same cell across conditions" else "", ". Black diamonds and bars show ",
  if (summary_display == "model") paste0("model-estimated ", ifelse(log_transform, "geometric means", "means"), " with 95% CI")
  else "medians with interquartile ranges",
  ". ", if (multi && show_p_on_plot) paste0("Brackets show ", adjust_name, "-adjusted p < ", p_threshold, ". ") else "",
  if (multi && comparison_plot != "none") paste0("[Panel/Figure Y.] Differences between conditions with ", adjust_name,
                                                 "-adjusted 95% CI; the dashed line marks no difference. ") else "",
  "Linear mixed-effects model (see Methods). ",
  paste(vapply(descriptives_tbl$Condition, function(cd) paste0(cd, ": ", n_text(cd)), ""), collapse = "; "), "."
)

caution_text <- c(
  "Before using this text:",
  "- Replace or check everything in [square brackets].",
  "- Report the Comparisons p-values, not the Fixed effects table.",
  "- Give n as mice / slices / cells, as above. The comparisons rest mainly on slices and mice, not on the number of cells.",
  if (multi && any(descriptives_tbl$n_mice < max(descriptives_tbl$n_mice)))
    "- Some conditions have fewer mice than others; say so (e.g. a missing washout recording)."
)

# 9. Plot ----------------------------------------------------------------------
# Coloured points/lines = individual cells (lines join the same cell across
# conditions, colours = mouse). Black marker + bars = summary_display.
#
# NOTE on "median_iqr": the median and IQR describe the SPREAD of the raw data;
# the mean +/- CI describes the UNCERTAINTY of an estimated mean. They answer
# different questions, and the p-values above always come from the mixed model,
# whichever summary you display.
if (summary_display == "model") {
  plot_summary <- data.frame(
    Condition = estimates_tbl$Condition,
    center    = estimates_tbl[[estimate_label]],
    lo        = estimates_tbl$CI_lower,
    hi        = estimates_tbl$CI_upper
  )
  plot_caption <- paste0("Black: model-estimated ",
                         ifelse(log_transform, "geometric mean", "mean"),
                         " \u00B1 95% CI of each condition (includes mouse and slice variation;\n",
                         "not the comparison -- see the comparisons panel/table)")
} else {
  plot_summary <- data.frame(
    Condition = descriptives_tbl$Condition,
    center    = descriptives_tbl$Median,
    lo        = descriptives_tbl$Q1,
    hi        = descriptives_tbl$Q3
  )
  plot_caption <- "Black: median with interquartile range (raw data)"
}
lv_all <- cond_levels
plot_summary$Condition <- factor(plot_summary$Condition, levels = lv_all)

# p-value annotations: only p below p_threshold is drawn.
#   * brackets between conditions for the comparisons table (adjusted p)
#   * "vs 0" labels above a condition when compare_to_zero = TRUE
# Same precision as the tables (fmt_p), so plot and tables show identical values.
label_p   <- function(p) {
  txt <- fmt_p(p)
  ifelse(startsWith(txt, "<"), paste0("p ", sub("<", "< ", txt, fixed = TRUE)), paste0("p = ", txt))
}
nudge_x   <- 0.22                                  # black markers sit this far right
to_axis   <- if (log_transform) log10 else function(x) x
from_axis <- if (log_transform) function(x) 10^x else function(x) x

annot_layers <- list()

if (show_p_on_plot) {
  
  rng <- diff(range(to_axis(df_long$Value)))
  if (!is.finite(rng) || rng == 0) rng <- 1
  gap  <- 0.075 * rng       # vertical spacing between stacked labels/brackets
  tick <- 0.02 * rng        # bracket end-tick length
  y_top <- to_axis(max(df_long$Value))
  
  # (a) "vs 0" labels
  vs0_shown <- FALSE
  if ("p_vs_0" %in% names(estimates_tbl)) {
    v0 <- estimates_tbl[!is.na(estimates_tbl$p_vs_0) & estimates_tbl$p_vs_0 < p_threshold, , drop = FALSE]
    if (nrow(v0) > 0) {
      vs0_shown <- TRUE
      annot_layers <- c(annot_layers, list(
        geom_text(data = data.frame(x = match(v0$Condition, lv_all) + nudge_x,
                                    y = from_axis(y_top + gap * 0.6),
                                    label = paste("vs 0:", label_p(v0$p_vs_0))),
                  aes(x = x, y = y, label = label),
                  size = 3.1, inherit.aes = FALSE)
      ))
    }
  }
  
  # (b) brackets
  n_sig <- 0
  if (!is.null(comparison_pairs)) {
    sig <- comparison_pairs[!is.na(comparison_pairs$p) & comparison_pairs$p < p_threshold, , drop = FALSE]
    n_sig <- nrow(sig)
    if (n_sig > 0) {
      sig$xa <- match(sig$A, lv_all) + nudge_x
      sig$xb <- match(sig$B, lv_all) + nudge_x
      sig    <- sig[order(abs(sig$xa - sig$xb), pmin(sig$xa, sig$xb)), , drop = FALSE]  # narrow first
      y_axis <- y_top + gap * (if (vs0_shown) 1.8 else 1) + (seq_len(n_sig) - 1) * gap * 1.6
      y_pl   <- from_axis(y_axis)
      y_tk   <- from_axis(y_axis - tick)
      
      seg <- rbind(
        data.frame(x = sig$xa, xend = sig$xb, y = y_pl, yend = y_pl),
        data.frame(x = sig$xa, xend = sig$xa, y = y_tk, yend = y_pl),
        data.frame(x = sig$xb, xend = sig$xb, y = y_tk, yend = y_pl)
      )
      lab <- data.frame(x = (sig$xa + sig$xb) / 2, y = y_pl, label = label_p(sig$p))
      
      annot_layers <- c(annot_layers, list(
        geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                     linewidth = 0.5, inherit.aes = FALSE),
        geom_text(data = lab, aes(x = x, y = y, label = label),
                  vjust = -0.35, size = 3.3, inherit.aes = FALSE)
      ))
    }
  }
  
  plot_caption <- paste0(
    plot_caption, "\n",
    if (!is.null(comparison_pairs)) {
      if (n_sig > 0) paste0("Brackets: comparisons with p < ", p_threshold, " (", adjust_name, ")")
      else paste0("No comparisons with p < ", p_threshold, " (", adjust_name, ")")
    } else if (vs0_shown) paste0("Labels: p < ", p_threshold, " vs 0")
    else paste0("No p < ", p_threshold, " to show")
  )
}

plot_caption <- paste0(plot_caption, "\nSlice x condition: ", slice_desc,
                       "; cell-level SD: ", cell_desc)

# Captions are wrapped to fit the 6.5-inch-wide image instead of being cut off.
if (any(df_long$Outlier)) {
  plot_caption <- paste0(plot_caption, "\nOpen diamonds: ROUT cell outliers (flagged, kept in the analysis)")
}
if (!is.null(flagged_slice_pts)) {
  plot_caption <- paste0(plot_caption, "\nOpen squares: slice means flagged by ROUT (kept; check imaging notes)")
}

wrap_caption <- function(txt, width = 95) {
  paste(unlist(lapply(strsplit(txt, "\n", fixed = TRUE)[[1]], strwrap, width = width)), collapse = "\n")
}
caption_lines <- function(txt) length(strsplit(txt, "\n", fixed = TRUE)[[1]])

# Colours (all plots): colour-blind-safe Okabe-Ito hues in a fixed order
# (checked for colour-blind separation). Yellow is left out (too faint on white)
# and black is kept for the means and CIs. With more mice than colours, a
# different shape is added so no two mice look the same.
okabe_ito    <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")
mouse_levels <- levels(df_long$Mouse)
mouse_index  <- setNames(seq_along(mouse_levels), mouse_levels)
shape_group  <- pmin((mouse_index - 1) %/% length(okabe_ito) + 1, 5)
if (is.null(pub_mouse_colours)) {
  mouse_colours <- setNames(okabe_ito[(mouse_index - 1) %% length(okabe_ito) + 1], names(mouse_index))
} else {
  mouse_colours <- pub_mouse_colours
}
mouse_shapes_fill  <- setNames(c(21, 22, 24, 23, 25)[shape_group], names(mouse_index))  # outlined (publication)
mouse_shapes_solid <- setNames(c(16, 15, 17, 18, 8)[shape_group], names(mouse_index))   # solid (notebook plot)

# Change measures get a dashed line at 0 (see change_measure in SETTINGS).
if (identical(change_measure, "auto")) {
  is_change <- grepl("change", csv_name, ignore.case = TRUE) || compare_to_zero
} else {
  is_change <- isTRUE(change_measure)
}

p <- ggplot(df_long, aes(x = Condition, y = Value)) +
  # Reference line at 0 when testing against 0.
  { if ((compare_to_zero || is_change) && !log_transform) geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") } +
  # Lines join the same cell across conditions (nothing to join if only one).
  { if (multi) geom_line(aes(group = interaction(Mouse, Slice, Cell), colour = Mouse), alpha = 0.35) } +
  geom_point(aes(colour = Mouse, shape = Mouse), size = 1.8, alpha = 0.75) +
  { if (any(df_long$Outlier)) geom_point(data = df_long[df_long$Outlier, , drop = FALSE], aes(x = Condition, y = Value),
                                         shape = 5, size = 2.6, stroke = 0.8, colour = "black", inherit.aes = FALSE) } +
  { if (!is.null(flagged_slice_pts)) geom_point(data = flagged_slice_pts, aes(x = Condition, y = Value),
                                                shape = 0, size = 3.4, stroke = 0.9, colour = "black", inherit.aes = FALSE) } +
  geom_errorbar(data = plot_summary,
                aes(x = Condition, ymin = lo, ymax = hi),
                width = 0.12, linewidth = 0.7,
                position = position_nudge(x = 0.22), inherit.aes = FALSE) +
  geom_point(data = plot_summary,
             aes(x = Condition, y = center),
             size = 3.2, shape = 18,
             position = position_nudge(x = 0.22), inherit.aes = FALSE) +
  labs(title = tools::file_path_sans_ext(csv_name),
       x = NULL,
       y = if (log_transform) paste0(y_label, " (log scale)") else y_label,
       caption = wrap_caption(plot_caption)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.05),
        plot.caption.position = "plot") +
  annot_layers

# Extra headroom so brackets and their labels are not clipped.
y_expand <- expansion(mult = c(0.05, 0.10))
p <- p + if (log_transform) scale_y_log10(expand = y_expand) else scale_y_continuous(expand = y_expand)
p <- p + scale_colour_manual(values = mouse_colours) + scale_shape_manual(values = mouse_shapes_solid)
if (y_include_zero && !log_transform) p <- p + expand_limits(y = 0)

if (show_plots) print(p)

# Comparisons panel ------------------------------------------------------------
# One row per comparison: the difference (ratio if log_transform) with its
# multiplicity-adjusted 95% CI and adjusted p -- the same numbers as the
# Comparisons table. A CI that crosses the dashed line (0, or 1 for ratios)
# matches p >= 0.05 (exactly for Tukey/Dunnett/Bonferroni; for Holm the CI is
# Bonferroni, so it can be slightly wider than the Holm p implies).
p_comp <- NULL
if (multi && comparison_plot != "none") {
  comp_df <- data.frame(
    Comparison = factor(comparisons_tbl$Comparison, levels = rev(comparisons_tbl$Comparison)),
    est = comparisons_tbl[[effect_label]],
    lo  = comparisons_tbl$CI_lower,
    hi  = comparisons_tbl$CI_upper,
    p   = comparisons_tbl$p_adj
  )
  comp_df$label <- paste0(label_p(comp_df$p), ", df = ", fmt_num(comparisons_tbl$DF))
  null_value <- if (log_transform) 1 else 0
  
  p_comp <- ggplot(comp_df, aes(y = Comparison, x = est)) +
    geom_vline(xintercept = null_value, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, linewidth = 0.7, orientation = "y") +
    geom_point(size = 3.2, shape = 18) +
    geom_text(aes(x = hi, label = label), hjust = -0.1, size = 3.2) +
    labs(title = paste0("Comparisons (", ifelse(log_transform, "ratio", "difference"),
                        ", 95% CI, ", adjust_label, "-adjusted)"),
         x = if (log_transform) paste0("Ratio of ", y_label, " (log scale)") else paste0("Difference in ", y_label),
         y = NULL,
         caption = paste0("CI crossing the dashed line = not significant at 0.05 after correction")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  
  # Room on the right for the p / df labels.
  x_expand <- expansion(mult = c(0.05, 0.45))
  p_comp <- p_comp + if (log_transform) scale_x_log10(expand = x_expand) else scale_x_continuous(expand = x_expand)
  
  if (show_plots) print(p_comp)
}

# Image heights: the main plot grows with its caption so nothing is cut off;
# the comparisons image grows with the number of comparisons.
main_height <- 4.6 + 0.16 * caption_lines(wrap_caption(plot_caption))
comp_height <- if (is.null(p_comp)) 0 else 1.4 + 0.45 * nrow(comparisons_tbl)

# 9b. Publication figure for Illustrator ---------------------------------------
# One file: the SuperPlot (Lord et al. 2020, J Cell Biol) with a small key to
# its right, and a compact comparisons strip below. Every cell is shown, but the slice means --
# the units the comparisons rest on -- are drawn over them, coloured by mouse,
# and the model mean with its CI (or SEM) is drawn on top of everything.
p_pub <- p_comp_pub <- p_key <- pub_figure <- NULL
pub_total_height_mm <- 0
# p-values in the publication figure are written without spaces: "p=0.156", "p<0.001".
label_p_pub <- function(p) gsub(" ", "", label_p(p), fixed = TRUE)
# "u" before a unit -> micro sign (pub_micro). Only a "u" that does not follow a
# letter counts, so words such as "sum" or "Bumetanide" are left alone.
micro <- function(x) {
  if (!pub_micro || is.null(x)) return(x)
  gsub("(?<![A-Za-z])u(?=(?:mol|Hz|M|m|L|l|g|s|A|V|F|W)(?![a-z]))", "\u03bc", x, perl = TRUE)
}
# Ion charges as superscripts (pub_ion_superscript): "Ca2+", "Cl-", "[K+]o".
# A charge counts only right after an ion name that does not follow another
# letter, so "KO-" or "Drug-B" stay as they are.
ion_pattern <- paste0("(?<![A-Za-z])(HCO3|NH4|SO4|PO4|NO3|OH|Ca|Cl|Na|Mg|Zn|Li|Cs|Rb|Ba|Sr|Fe|Cu|Mn|Cd|Br|K|H|I|F)",
                      "([0-9]?[+-])(?![A-Za-z0-9])")
# Text for the figure: micro sign, and ion charges raised (as an R plotmath
# expression, which keeps every piece as editable text in the SVG).
pub_text <- function(x, bold = FALSE) {
  if (is.null(x)) return(x)
  x <- micro(as.character(x))
  if (!pub_ion_superscript || !any(grepl(ion_pattern, x, perl = TRUE))) return(x)
  one_line <- function(s) {
    s <- gsub("([\"\\\\])", "\\\\\\1", s)   # escape quotes and backslashes
    paste0("\"", gsub(ion_pattern, "\\1\"^\"\\2\"*\"", s, perl = TRUE), "\"")
  }
  expr <- vapply(x, function(s) {
    lines <- unname(vapply(strsplit(s, "\n", fixed = TRUE)[[1]], one_line, ""))
    e <- if (length(lines) == 0) "\"\"" else if (length(lines) == 1) lines else
      paste0("atop(displaystyle(", lines[1], "), displaystyle(", paste(lines[-1], collapse = "*\" \"*"), "))")
    if (bold) paste0("bold(", e, ")") else e
  }, "", USE.NAMES = FALSE)
  parse(text = expr)
}
# Condition names as shown in the publication figure (pub_condition_labels).
pub_lab <- function(x) {
  x <- as.character(x)
  if (is.null(pub_condition_labels)) return(micro(x))
  hit <- x %in% names(pub_condition_labels)
  x[hit] <- unname(pub_condition_labels[x[hit]])
  micro(x)
}
# Width of the longest label in mm (for the key), measured with the font itself
# when the systemfonts package is there, otherwise estimated.
text_width_mm <- function(labels, pt) {
  labels <- unlist(strsplit(labels, "\n", fixed = TRUE))
  w <- NULL
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    w <- tryCatch(systemfonts::string_width(labels, family = pub_font, size = pt, res = 72),
                  error = function(e) NULL)
  }
  if (is.null(w) || length(w) == 0 || any(!is.finite(w))) w <- nchar(labels) * pt * 0.55
  max(w) * 25.4 / 72 * 1.05
}
pub_files <- pub_previews <- pub_save_errors <- character(0)

# Unit converters: ggplot measures lines, point outlines and geom text in its
# own units; these turn points (pt) into those units.
pt_lw     <- function(pt) pt * 96 / 72 / ggplot2::.pt               # line width
pt_stroke <- function(pt) pt * 96 / 72 * 2 / ggplot2::.stroke       # point outline
pt_text   <- function(pt) pt / ggplot2::.pt                         # geom_text size

if (make_pub_figures) {
  # Windows' on-screen plot window and PNG writer only know fonts registered
  # here; without this they warn "font family not found in Windows font
  # database". (SVG and PDF writers find the font by themselves.)
  if (.Platform$OS.type == "windows" && !(pub_font %in% names(grDevices::windowsFonts()))) {
    do.call(grDevices::windowsFonts, setNames(list(grDevices::windowsFont(pub_font)), pub_font))
  }
  
  pub_theme <- theme_classic(base_size = pub_tick_pt, base_family = pub_font) +
    theme(
      axis.line         = element_line(linewidth = pt_lw(pub_axis_pt), colour = "black"),
      axis.ticks        = element_line(linewidth = pt_lw(pub_axis_pt), colour = "black"),
      axis.ticks.length = unit(1, "mm"),
      axis.text         = element_text(size = pub_tick_pt, colour = "black", family = pub_font),
      axis.title        = element_text(size = pub_title_pt, colour = "black", family = pub_font, face = "bold"),
      axis.title.y      = element_text(size = pub_title_pt, colour = "black", family = pub_font, face = "bold"),
      legend.position   = "none",
      plot.margin       = margin(1, 1, 1, 1, "mm")
    )
  
  pub_slice_means <- df_long |>
    group_by(Mouse, Slice, Condition) |>
    summarise(Value = if (log_transform) 10^mean(log10(Value)) else mean(Value), .groups = "drop")
  
  # Cells: faint dots, spread sideways so they don't overlap.
  use_beeswarm <- pub_cell_layout == "beeswarm" && requireNamespace("ggbeeswarm", quietly = TRUE)
  if (pub_cell_layout == "beeswarm" && !use_beeswarm) {
    message("ggbeeswarm is not installed -- using jitter for the cell dots (install.packages(\"ggbeeswarm\") for beeswarm).")
  }
  cell_layer <- if (use_beeswarm) {
    ggbeeswarm::geom_quasirandom(aes(fill = Mouse, colour = Mouse, shape = Mouse), size = pub_cell_size,
                                 stroke = pt_stroke(pub_dot_pt), alpha = pub_cell_alpha, width = 0.15)
  } else {
    geom_point(aes(fill = Mouse, colour = Mouse, shape = Mouse), size = pub_cell_size,
               stroke = pt_stroke(pub_dot_pt), alpha = pub_cell_alpha,
               position = position_jitter(width = 0.12, height = 0, seed = 1))
  }
  
  # Mean line and error bars: model mean with its 95% CI (default) or +/- 1 SEM.
  pub_sum <- plot_summary
  if (pub_error_bar == "sem" && summary_display == "model") {
    pub_sum$lo <- back(emm_df$emmean - emm_df$SE)
    pub_sum$hi <- back(emm_df$emmean + emm_df$SE)
  }
  summary_layers <- list()
  if (pub_summary == "mean_ci") {
    summary_layers <- c(summary_layers, list(
      geom_errorbar(data = pub_sum, aes(x = Condition, ymin = lo, ymax = hi),
                    width = pub_ci_cap_width, linewidth = pt_lw(pub_ci_pt), colour = "black",
                    position = position_nudge(x = pub_summary_offset), inherit.aes = FALSE)))
  }
  if (pub_summary %in% c("mean_ci", "mean")) {
    summary_layers <- c(summary_layers, list(
      geom_errorbar(data = pub_sum, aes(x = Condition, ymin = center, ymax = center),
                    width = pub_mean_width, linewidth = pt_lw(pub_mean_pt), colour = "black",
                    position = position_nudge(x = pub_summary_offset), inherit.aes = FALSE)))
  }
  
  # p-value brackets (same comparisons and threshold as the other plots).
  pub_brackets <- list()
  y_span <- diff(range(to_axis(df_long$Value)))
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  bracket_top <- NA
  if (show_p_on_plot && !is.null(comparison_pairs)) {
    sig <- comparison_pairs[!is.na(comparison_pairs$p) & comparison_pairs$p < p_threshold, , drop = FALSE]
    if (nrow(sig) > 0) {
      gap   <- 0.09 * y_span
      tick  <- 0.02 * y_span
      y_top <- to_axis(max(c(df_long$Value, pub_sum$hi), na.rm = TRUE))
      sig$xa <- match(sig$A, cond_levels)
      sig$xb <- match(sig$B, cond_levels)
      sig    <- sig[order(abs(sig$xa - sig$xb), pmin(sig$xa, sig$xb)), , drop = FALSE]
      y_ax   <- y_top + gap + (seq_len(nrow(sig)) - 1) * gap * 1.5
      seg <- rbind(
        data.frame(x = sig$xa, xend = sig$xb, y = from_axis(y_ax), yend = from_axis(y_ax)),
        data.frame(x = sig$xa, xend = sig$xa, y = from_axis(y_ax - tick), yend = from_axis(y_ax)),
        data.frame(x = sig$xb, xend = sig$xb, y = from_axis(y_ax - tick), yend = from_axis(y_ax))
      )
      lab <- data.frame(x = (sig$xa + sig$xb) / 2, y = from_axis(y_ax), label = label_p_pub(sig$p))
      pub_brackets <- list(
        geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                     linewidth = pt_lw(pub_axis_pt), inherit.aes = FALSE),
        geom_text(data = lab, aes(x = x, y = y, label = label), vjust = -0.3,
                  size = pt_text(pub_p_pt), family = pub_font, inherit.aes = FALSE)
      )
      bracket_top <- from_axis(max(y_ax) + 0.08 * y_span)   # room for the p-value text
    }
  }
  
  # y-axis title of the SuperPlot.
  pub_y_base <- if (!is.null(pub_y_title)) pub_y_title else y_label
  pub_y_text <- micro(if (is.null(pub_y_title) && log_transform) paste0(y_label, " (log scale)") else pub_y_base)
  
  p_pub <- ggplot(df_long, aes(x = Condition, y = Value)) +
    { if ((compare_to_zero || is_change) && !log_transform) geom_hline(yintercept = 0, linetype = pub_dash,
                                                                       linewidth = pt_lw(pub_axis_pt), colour = "grey40") } +
    cell_layer +
    # Lines join the same slice across conditions: the paired comparison, made visible.
    { if (multi) geom_line(data = pub_slice_means,
                           aes(x = Condition, y = Value, group = interaction(Mouse, Slice), colour = Mouse),
                           linewidth = pt_lw(pub_slice_line_pt), inherit.aes = FALSE) } +
    geom_point(data = pub_slice_means, aes(x = Condition, y = Value, fill = Mouse, shape = Mouse),
               colour = "black", size = pub_slice_size, stroke = pt_stroke(pub_dot_pt), inherit.aes = FALSE) +
    summary_layers +   # drawn last, so the mean and CI sit on top of the data
    pub_brackets +
    scale_fill_manual(values = mouse_colours) +
    scale_colour_manual(values = mouse_colours) +
    scale_shape_manual(values = mouse_shapes_fill) +
    scale_x_discrete(labels = function(b) pub_text(pub_lab(b))) +
    labs(x = if (is.null(pub_x_title) || pub_x_title == "") NULL else pub_text(pub_x_title, bold = TRUE),
         y = pub_text(pub_y_text, bold = TRUE)) +
    coord_cartesian(clip = "off") +
    pub_theme
  
  # y axis: starts and ends on labelled ticks (and includes 0 if y_include_zero).
  if (log_transform) {
    p_pub <- p_pub + scale_y_log10(expand = expansion(mult = c(0.05, 0.08)))
  } else {
    y_vals   <- c(df_long$Value, pub_sum$lo, pub_sum$hi, bracket_top, if (y_include_zero) 0)
    y_breaks <- pretty(range(y_vals, na.rm = TRUE), n = 4)
    p_pub <- p_pub + scale_y_continuous(breaks = y_breaks, limits = range(y_breaks), expand = c(0, 0))
  }
  
  # Comparisons strip: one row per comparison -- mean line = the difference,
  # bars = its adjusted 95% CI, p on the right. Compact by design.
  if (multi && pub_comp_height_mm > 0) {
    comp_names <- if (is.null(pub_condition_labels)) comparisons_tbl$Comparison else
      gsub("\n", " ", paste(pub_lab(comparison_pairs$A), if (log_transform) "/" else "-", pub_lab(comparison_pairs$B)), fixed = TRUE)
    cd <- data.frame(
      Comparison = factor(comp_names, levels = rev(unique(comp_names))),
      est = comparisons_tbl[[effect_label]], lo = comparisons_tbl$CI_lower, hi = comparisons_tbl$CI_upper,
      label = label_p_pub(comparisons_tbl$p_adj)
    )
    p_comp_pub <- ggplot(cd, aes(y = Comparison, x = est)) +
      geom_vline(xintercept = if (log_transform) 1 else 0, linetype = pub_dash,
                 linewidth = pt_lw(pub_axis_pt), colour = "grey40") +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.35, linewidth = pt_lw(pub_ci_pt), orientation = "y") +
      geom_errorbar(aes(xmin = est, xmax = est), width = 0.6, linewidth = pt_lw(pub_mean_pt), orientation = "y") +
      geom_text(aes(x = Inf, label = label), hjust = 1, size = pt_text(pub_p_pt), family = pub_font) +
      labs(x = NULL, y = NULL) +
      scale_y_discrete(labels = pub_text) +
      pub_theme +
      theme(axis.text         = element_text(size = pub_small_pt, colour = "black", family = pub_font),
            axis.ticks.length = unit(0.6, "mm"),
            plot.margin       = margin(0.5, 1, 0.5, 1, "mm"))
    comp_x_expand <- expansion(mult = c(0.05, 0.6))   # room on the right for the p-values
    p_comp_pub <- p_comp_pub + if (log_transform) scale_x_log10(expand = comp_x_expand) else scale_x_continuous(expand = comp_x_expand)
  }
  
  # Key: what the symbols mean, in a small dashed box to the right of the SuperPlot.
  # Drawn in mm (x: 0 to key_w_mm, y: 0 to key_h_mm) so sizes stay exact. Its
  # width fits the longest label, so the SuperPlot keeps as much width as it can.
  key_h_mm <- 8
  key_w_mm <- 0
  if (is.null(pub_key_width_mm) || pub_key_width_mm > 0) {
    key_labels <- unname(pub_key_labels[c("cell", "slice", "mouse")])
    key_y   <- c(key_h_mm - 1.7, key_h_mm / 2, 1.7)   # rows: cell, slice, mouse colour
    key_dot <- 1.25                                   # x of the symbols
    key_txt <- 2.5                                    # x where the labels start
    key_w_mm <- if (is.null(pub_key_width_mm)) key_txt + text_width_mm(key_labels, pub_key_pt) + 0.6 else pub_key_width_mm
    p_key <- ggplot() +
      # White inside the box, so it reads on a dark background too.
      annotate("rect", xmin = 0.1, xmax = key_w_mm - 0.1, ymin = 0.3, ymax = key_h_mm - 0.3,
               fill = "white", colour = "black", linetype = pub_dash, linewidth = pt_lw(pub_axis_pt)) +
      annotate("point", x = key_dot, y = key_y[1], shape = 21, fill = "grey55", colour = "grey55",
               alpha = pub_cell_alpha, size = pub_cell_size, stroke = pt_stroke(pub_dot_pt)) +
      annotate("point", x = key_dot, y = key_y[2], shape = 21, fill = "grey55", colour = "black",
               size = pub_slice_size, stroke = pt_stroke(pub_dot_pt)) +
      annotate("point", x = key_dot, y = key_y[3], shape = 21, fill = unname(mouse_colours)[1],
               colour = "black", size = pub_cell_size * 1.3, stroke = pt_stroke(pub_dot_pt)) +
      annotate("text", x = key_txt, y = key_y, label = key_labels, hjust = 0,
               size = pt_text(pub_key_pt), family = pub_font) +
      scale_x_continuous(limits = c(0, key_w_mm), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, key_h_mm), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void(base_family = pub_font) +
      theme(plot.margin = margin(0, 0, 0, 0, "mm"))
  }
  
  # Put the parts into ONE figure (one file to place in Illustrator):
  # SuperPlot with the key to its right, comparisons strip across the bottom.
  if (is.null(p_key)) {
    top_grob <- ggplotGrob(p_pub)
  } else {
    # The key sits pub_key_nudge_mm over the SuperPlot's empty right edge, and the
    # SuperPlot has no right margin, so there is little space between them.
    key_ov    <- min(max(pub_key_nudge_mm, 0), key_w_mm)
    main_w_mm <- pub_width_mm - key_w_mm + key_ov
    key_col   <- gridExtra::arrangeGrob(ggplotGrob(p_key), grid::nullGrob(), ncol = 1,
                                        heights = grid::unit(c(key_h_mm, 1), c("mm", "null")))
    top_grob  <- gtable::gtable(widths  = grid::unit(c(main_w_mm - key_ov, key_ov, key_w_mm - key_ov), "mm"),
                                heights = grid::unit(1, "null"))
    top_grob  <- gtable::gtable_add_grob(top_grob, ggplotGrob(p_pub + theme(plot.margin = margin(1, 0, 1, 1, "mm"))),
                                         t = 1, l = 1, r = 2, z = 1, name = "superplot")
    top_grob  <- gtable::gtable_add_grob(top_grob, key_col, t = 1, l = 2, r = 3, z = 2, name = "key")
  }
  pub_parts   <- list(top_grob)
  pub_heights <- pub_main_height_mm
  if (!is.null(p_comp_pub)) { pub_parts <- c(pub_parts, list(ggplotGrob(p_comp_pub))); pub_heights <- c(pub_heights, pub_comp_height_mm) }
  pub_figure <- gridExtra::arrangeGrob(grobs = pub_parts, ncol = 1, heights = grid::unit(pub_heights, "mm"))
  pub_total_height_mm <- sum(pub_heights)
  if (show_plots) {
    grid::grid.newpage()
    grid::grid.draw(pub_figure)
  }
  
  # The figure legend in the writing template describes this figure.
  legend_text <- paste0(
    "[Figure X.] ", y_label, " in each condition. Small dots: individual cells; large circles: slice means; ",
    "colours: mice",
    if (multi) "; lines join the same slice across conditions" else "",
    if (pub_summary == "mean_ci") paste0(". Horizontal line and bars: model-estimated ",
                                         ifelse(log_transform, "geometric mean", "mean"),
                                         if (pub_error_bar == "sem") " ± SEM" else " with 95% CI")
    else if (pub_summary == "mean") paste0(". Horizontal line: model-estimated ",
                                           ifelse(log_transform, "geometric mean", "mean"))
    else "",
    if (is_change && !log_transform) ". Dashed line: no change" else "",
    ". ", if (multi && show_p_on_plot) paste0("Brackets: ", adjust_name, "-adjusted p < ", p_threshold, ". ") else "",
    if (multi && !is.null(p_comp_pub)) paste0("Bottom: difference between conditions (vertical line) with ", adjust_name,
                                              "-adjusted 95% CI (bars) and p; the dashed line marks no difference. ") else "",
    "Linear mixed-effects model (see Methods). ",
    paste(vapply(descriptives_tbl$Condition, function(cd) paste0(cd, ": ", n_text(cd)), ""), collapse = "; "), "."
  )
}

# Saves the publication figure in every format in pub_formats, plus a PNG
# preview for the note. Returns the file names that were saved.
save_pub <- function(plot, stem, height_mm) {
  saved <- character(0)
  for (fmt in pub_formats) {
    f <- paste0(stem, ".", fmt)
    ok <- tryCatch({
      path <- file.path(pub_dir, f)
      if (fmt == "pdf") {
        ggsave(path, plot, width = pub_width_mm, height = height_mm, units = "mm", bg = "white", device = grDevices::cairo_pdf)
      } else if (fmt == "svg") {
        if (!requireNamespace("svglite", quietly = TRUE)) {
          stop("the svglite package is needed for SVG -- run install.packages(\"svglite\") once")
        }
        # fix_text_size = FALSE keeps labels as plain, editable Arial text in Illustrator
        # (otherwise each label is stretched to a fixed length).
        ggsave(path, plot, width = pub_width_mm, height = height_mm, units = "mm", bg = "white",
               device = svglite::svglite, fix_text_size = FALSE)
      } else {
        ggsave(path, plot, width = pub_width_mm, height = height_mm, units = "mm", bg = "white",
               device = grDevices::cairo_ps, fallback_resolution = 600)
      }
      TRUE
    }, error = function(e) {
      pub_save_errors <<- c(pub_save_errors, paste0(f, ": ", conditionMessage(e)))
      warning("Could not save ", f, ": ", conditionMessage(e))
      FALSE
    })
    if (ok) saved <- c(saved, f)
  }
  # PNG preview: ragg (if installed) finds system fonts directly and draws more cleanly.
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(plot_dir, paste0(stem, ".png")), plot, width = pub_width_mm, height = height_mm,
         units = "mm", bg = "white", dpi = 600, device = png_device)
  saved
}

# Residual diagnostics: use these to judge log_transform and model_engine.
# Left: residuals should be a shapeless cloud (no fan/funnel).
# Right: points should follow the line (roughly normal residuals).
# For the nlme model the residuals are standardized (divided by each
# condition's estimated SD), so they should look equally spread if the
# per-condition variances did their job.

scale_note  <- if (log_transform) "(log10 scale)" else "(raw scale)"
resid_type  <- if (inherits(model, "lme")) "pearson" else "response"
resid_label <- if (inherits(model, "lme")) "Standardized residuals" else "Residuals"

normality_of <- function(x) {
  if (exists("labstats_normality", inherits = TRUE)) return(labstats_normality(x))
  n    <- length(x)
  sw   <- if (n >= 3 && n <= 5000) shapiro.test(x) else NULL
  dago <- tryCatch(fBasics::dagoTest(x), error = function(e) NULL)
  data.frame(
    n       = n,
    SW_W    = if (!is.null(sw))   unname(sw$statistic)              else NA,
    SW_p    = if (!is.null(sw))   sw$p.value                        else NA,
    DAgo_K2 = if (!is.null(dago)) unname(dago@test$statistic[1])    else NA,
    DAgo_p  = if (!is.null(dago)) unname(dago@test$p.value[1])      else NA
  )
}

# ============================================================
# VIEW 1: residual diagnostics -- checks the MODEL'S assumption
# ============================================================
# fitted()/resid() return values in the same row order as df_long for both
# lmer and lme, so Condition can be taken straight from df_long.
diag_df <- data.frame(fitted    = as.numeric(fitted(model)),
                      resid     = as.numeric(resid(model, type = resid_type)),
                      Condition = df_long$Condition)

p_resid <- ggplot(diag_df, aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6) +
  labs(title = paste("Residuals vs fitted", scale_note), x = "Fitted", y = resid_label) +
  theme_bw(base_size = 11)

p_qq_resid <- ggplot(diag_df, aes(sample = resid)) +
  stat_qq(alpha = 0.6) + stat_qq_line(fullrange = TRUE) +
  labs(title = paste("Residual QQ plot", scale_note), x = "Theoretical quantiles", y = "Actual quantiles") +
  theme_bw(base_size = 11) +
  coord_flip()

resid_by_condition <- ggplot(diag_df, aes(sample = resid)) +
  stat_qq(alpha = 0.6) + stat_qq_line(fullrange = TRUE) +
  facet_wrap(~ Condition, scales = "free") +
  labs(title = paste("Residual QQ plot by condition", scale_note),
       x = "Theoretical quantiles", y = "Actual quantiles") +
  theme_bw(base_size = 11) +
  coord_flip()

residual_diagnostics <- arrangeGrob(p_resid, p_qq_resid, ncol = 2)

resid_normality_tbl <- normality_of(diag_df$resid)

# ==================================================================
# VIEW 2: raw values by condition -- the check you're used to in Prism
# ==================================================================
qq_df <- df_long
qq_df$QQValue <- if (log_transform) log10(qq_df$Value) else qq_df$Value

condition_qq <- ggplot(qq_df, aes(sample = QQValue)) +
  stat_qq(alpha = 0.6) + stat_qq_line(fullrange = TRUE) +
  facet_wrap(~ Condition, scales = "free") +
  labs(title = paste("Prism: QQ plot by condition", scale_note), x = "Theoretical quantiles", y = "Actual quantiles") +
  theme_bw(base_size = 11) +
  coord_flip()

condition_normality_tbl <- qq_df |>
  group_by(Condition) |>
  group_modify(~ normality_of(.x$QQValue)) |>
  ungroup() |>
  mutate(Condition = as.character(Condition)) |>
  as.data.frame()

resid_row <- data.frame(Condition = "Model (residuals)", normality_of(diag_df$resid))
normality_tbl <- dplyr::bind_rows(resid_row, condition_normality_tbl)


base                   <- tools::file_path_sans_ext(csv_name)
diagnostics_file       <- paste0(base, "_diagnostics.png")
residual_per_condition <- paste0(base, "_residual_condition_diagnostics.png")
qq_by_condition_file   <- paste0(base, "_qq_by_condition.png")

if (show_diagnostics && save_outputs) {
  ggsave(file.path(plot_dir, diagnostics_file), residual_diagnostics, width = 8, height = 4, dpi = 300)
  ggsave(file.path(plot_dir, residual_per_condition), resid_by_condition, width = 8, height = 3.2, dpi = 300)
  ggsave(file.path(plot_dir, qq_by_condition_file), condition_qq, width = 8, height = 3.2, dpi = 300)
}

# 10. Obsidian writing function ------------------------------------------------
df_to_markdown <- function(df, header_names = NULL) {
  
  df <- fmt_tbl(df)   # same number formatting as the console
  
  header <- paste(if (is.null(header_names)) names(df) else header_names, collapse = " | ")
  separator <- paste(rep("---", ncol(df)), collapse = " | ")
  
  rows <- apply(df, 1, function(x)
    paste(trimws(as.character(x)), collapse = " | "))
  
  paste(
    c(
      paste0("| ", header, " |"),
      paste0("| ", separator, " |"),
      paste0("| ", rows, " |")
    ),
    collapse = "\n"
  )
}


# Obsidian links for the note's tables ---------------------------------------------
# Mouse and slice IDs in the outlier tables link to their notes in
# experiments_dir. A mouse note's name starts with the mouse's date (e.g.
# 2025.08.27-HIC-NaCl-Adult) and has "Type: Whole" in its frontmatter; a slice
# note is named exactly like the slice (e.g. 2025.08.27-S2). IDs with no note
# found stay as plain text. ("\|" is how a link's display name is written
# inside an Obsidian table.)
note_index <- NULL
if (!is.null(experiments_dir) && dir.exists(experiments_dir)) {
  md_files   <- list.files(experiments_dir, pattern = "\\.md$", recursive = TRUE, full.names = TRUE)
  note_index <- data.frame(name = tools::file_path_sans_ext(basename(md_files)), path = md_files,
                           stringsAsFactors = FALSE)
}
note_type <- function(path) {
  lines <- tryCatch(readLines(path, n = 60, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  hit <- grep("^Type:", lines, value = TRUE)
  if (length(hit) == 0) "" else trimws(gsub("[\"']", "", sub("^Type:", "", hit[1])))
}
mouse_link <- function(m) {
  if (is.null(note_index)) return(m)
  cand <- note_index[startsWith(note_index$name, m), , drop = FALSE]
  if (nrow(cand) == 0) return(m)
  whole <- vapply(cand$path, note_type, "", USE.NAMES = FALSE) == "Whole"
  cand  <- if (any(whole)) cand[whole, , drop = FALSE] else cand[!grepl("-S[0-9]+$", cand$name), , drop = FALSE]
  if (nrow(cand) == 0) m else paste0("[[", cand$name[1], "\\|", m, "]]")
}
slice_link <- function(s) {
  if (is.null(note_index) || !(s %in% note_index$name)) s else paste0("[[", s, "]]")
}
link_ids <- function(tbl) {
  if (is.null(tbl)) return(NULL)
  tbl$Mouse <- vapply(tbl$Mouse, mouse_link, "", USE.NAMES = FALSE)
  tbl$Slice <- vapply(tbl$Slice, slice_link, "", USE.NAMES = FALSE)
  tbl
}

# Normality table headers link to the notes explaining each test.
normality_header <- c("Condition", "n",
                      "[[Shapiro-Wilk\\|SW W]]", "[[Shapiro-Wilk\\|SW p]]",
                      "[[D'Agostino-Pearson K2 omnibus\\|DAgo K2]]", "[[D'Agostino-Pearson K2 omnibus\\|DAgo p]]")

# 11. Save ---------------------------------------------------------------------
if (save_outputs) {
  
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  plot_file <- paste0(base, ".png")
  
  ggsave(
    file.path(plot_dir, plot_file),
    p,
    width = 6.5,
    height = main_height,
    dpi = 300
  )
  
  comp_file <- paste0(base, "_comparisons.png")
  if (!is.null(p_comp) && comparison_plot == "separate") {
    ggsave(file.path(plot_dir, comp_file), p_comp, width = 6.5, height = comp_height, dpi = 300)
  }
  
  if (!is.null(pub_figure)) {
    dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
    pub_files    <- c(pub_files, save_pub(pub_figure, paste0(base, "_figure"), pub_total_height_mm))
    pub_previews <- c(pub_previews, paste0(base, "_figure.png"))
  }
  
  sample_size_lines <- paste0("- ", descriptives_tbl$n_cells)
  
  note_lines <- c(
    "---",
    "tags:",
    "  - 🖼-figure/my-data",
    paste0("Type: ", strsplit(base, "_")[[1]][1]),
    paste0('R_csv: "[[', csv_name, ']]"'),
    "Figure Link:",
    paste0('   - "[[', plot_file, ']]"'),
    if (length(pub_files) > 0) c("Illustrator_Figures:", paste0('   - "[[', pub_files, ']]"')),
    "Sample_Size_N:",
    sample_size_lines,
    paste0("Normality: ", if (log_transform) "Log transformed" else "Normal"),
    "Analysis_Test: Nested mixed-effects model",
    paste0('Model: "', model_letter, " (", if (engine_used == "nlme") "nlme" else "lmer", "): ",
           models_tbl$Assumes[models_tbl$Model == model_letter], '"'),
    paste0("Slice_x_Condition: ", slice_desc),
    paste0("Cell_variance: ", cell_desc),
    paste0(
      "Post_hoc_Test: ",
      if (is.null(comparisons_tbl)) "None" else adjust_name
    ),
    paste0(
      "Error_Bar: ",
      if (summary_display == "model")
        "95% CI"
      else
        "Median + IQR"
    ),
    paste0("Modified: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
    "---",
    "",
    "# Plot",
    "",
    paste0("![[", plot_file, "]]"),
    if (!is.null(p_comp) && comparison_plot == "separate") paste0("![[", comp_file, "]]"),
    if (length(pub_previews) > 0) c(
      "",
      "## Publication figures (for Illustrator)",
      "",
      paste0(pub_width_mm, " mm wide, ", pub_font, "; place at 100% in Illustrator so font sizes and strokes stay exact. Files: ",
             if (length(pub_files) > 0) paste0("[[", pub_files, "]]", collapse = ", ") else "none saved"),
      if (length(pub_save_errors) > 0) c("", paste0("**Not saved:** ", pub_save_errors)),
      "",
      paste0("![[", pub_previews, "]]")),
    "",
    "# Results",
    ""
  )
  
  if (!is.null(comparisons_tbl)) {
    note_lines <- c(note_lines,
                    paste0("## Comparisons (report these; p adjusted: ", adjust_label, ")"),
                    "",
                    paste0("Each slice is its own control here, so slice and mouse offsets cancel. ",
                           "DF reflects the number of slices behind each comparison (see decision 6)."),
                    "",
                    df_to_markdown(comparisons_tbl),
                    "",
                    if (!is.null(df_breakdown_tbl)) c("### Where the degrees of freedom come from", "", df_explain, "",
                                                      df_to_markdown(df_breakdown_tbl), "")
    )
  }
  
  note_lines <- c(note_lines,
                  paste0("## Model-estimated ", ifelse(log_transform, "geometric mean", "mean"), " per condition (95% CI)"),
                  "",
                  paste0("Describes each condition on its own: the CI includes all variation between mice and slices. ",
                         "**Not used for the comparisons** -- two overlapping CIs here do not mean two conditions are the same."),
                  "",
                  df_to_markdown(estimates_tbl),
                  "",
                  "## Descriptive statistics (raw data)",
                  "",
                  "SD, 95% CI and SEM in this table treat every cell as independent, like Prism's column statistics: descriptive only, and too narrow for inference. Use the model estimates and the Comparisons table for that.",
                  "",
                  df_to_markdown(descriptives_tbl),
                  "",
                  "## Outliers (ROUT)",
                  "",
                  "### Cells",
                  "",
                  outlier_note,
                  if (!is.null(outlier_tbl)) c("", df_to_markdown(link_ids(outlier_tbl))),
                  "",
                  "### Slices",
                  "",
                  slice_outlier_note,
                  if (!is.null(slice_outlier_tbl)) c("", df_to_markdown(link_ids(slice_outlier_tbl)))
  )
  
  # Evidence behind decisions 3 and 4.
  if (!is.null(slice_sd_tbl) || !is.null(cell_sd_tbl)) {
    note_lines <- c(note_lines, "", "# Checks behind the decisions", "",
                    paste0("For decisions 3 and 4 the model is fitted twice -- once assuming every condition has the same SD, ",
                           "once allowing a separate SD per condition -- and a likelihood-ratio test decides which to use. ",
                           "Each table shows the SDs from both versions side by side, then the test, then what was chosen."),
                    "",
                    paste0("*", sd_note, "*"))
    if (!is.null(slice_sd_tbl)) {
      note_lines <- c(note_lines, "",
                      "## Decision 3: slice-to-slice SD",
                      "",
                      "How much slices differ from each other in each condition, under each of the three slice structures.",
                      "",
                      "- SD_none: slices assumed to respond alike (no slice x condition term).",
                      "- SD_shared: slices respond differently, with one shared spread.",
                      "- SD_per_condition: slices respond differently, with a separate spread per condition.",
                      "",
                      df_to_markdown(slice_sd_tbl),
                      "",
                      slice_reading,
                      "",
                      if (!is.null(slice_test)) c(df_to_markdown(slice_test), ""),
                      slice_choice)
    }
    if (!is.null(cell_sd_tbl)) {
      note_lines <- c(note_lines, "",
                      "## Decision 4: cell-to-cell SD",
                      "",
                      "How much cells from the same slice differ from each other (in the same condition).",
                      "",
                      df_to_markdown(cell_sd_tbl),
                      "",
                      cell_reading,
                      "",
                      if (!is.null(variance_test)) c(df_to_markdown(variance_test), ""),
                      cell_choice)
    }
    note_lines <- c(note_lines, "", paste0("*", lr_glossary, "*"))
  }
  
  note_lines <- c(note_lines,
                  "",
                  "# Assumption checks",
                  "",
                  "Residuals should look roughly normal on the QQ plot. With ~200+ points the tests flag even tiny departures, so judge by the plot.",
                  "",
                  df_to_markdown(normality_tbl, normality_header)
  )
  
  if (show_diagnostics) {
    note_lines <- c(note_lines,
                    "",
                    paste0("![[", diagnostics_file, "]]"),
                    paste0("![[", residual_per_condition, "]]"),
                    paste0("![[", qq_by_condition_file, "]]")
    )
  }
  
  note_lines <- c(note_lines,
                  "",
                  "# Model details",
                  "",
                  "## Model-fitting checks",
                  "",
                  sub("^  - ", "- ", fit_check_lines),
                  "",
                  "## Random effects (where the variability sits)",
                  "",
                  df_to_markdown(random_tbl),
                  "",
                  "## Fixed effects (raw coefficients)",
                  "",
                  "Unadjusted p-values for the model's coefficients (see section 6 of the script). Report the Comparisons table instead.",
                  "",
                  df_to_markdown(fixed_tbl),
                  if (!is.null(df_note)) c("", df_note)
  )
  
  note_lines <- c(note_lines,
                  "",
                  "# Analysis decisions",
                  "",
                  "Each row is one choice about the model. Every result above comes from this exact model.",
                  "",
                  df_to_markdown(decisions_tbl),
                  "",
                  paste0("**In words:** ", decision_summary),
                  if (length(run_choices) > 0) c("", "**Chosen when running:**", "", paste0("- ", run_choices)),
                  "",
                  "## The six possible models",
                  "",
                  paste0("Decisions 3 and 4 are separate choices, so there are six possible models. All six use the same ",
                         "Mouse/Slice/Cell nesting and the same comparisons; they differ only in how variability is modelled. ",
                         "The tests in Checks behind the decisions pick one automatically (when set to \"auto\"). This analysis used model ", model_letter, "."),
                  "",
                  df_to_markdown(models_tbl)
  )
  
  note_lines <- c(note_lines,
                  "",
                  "# Writing template",
                  "",
                  "Draft text filled in from this analysis. Edit before use.",
                  "",
                  "## Methods",
                  "",
                  methods_text,
                  "",
                  "## Results",
                  "",
                  means_text,
                  "",
                  if (length(results_text) > 0) c(paste0("- ", results_text), ""),
                  "## Figure legend",
                  "",
                  legend_text,
                  "",
                  caution_text
  )
  
  writeLines(note_lines, file.path(results_dir, paste0(base, ".md")))
  
  cat("Saved markdown report and plot to:", results_dir, "\n")
  
  # 11b. HTML report (for lab members without Obsidian) ---------------------------
  # The same content as the Obsidian note as ONE self-contained web page (images
  # built in), shown in RStudio's Viewer. Flip through it section by section
  # with the arrow keys or the buttons, or show everything at once. It is not
  # saved: it lives in R's temporary folder until R closes.
  if (html_report) {
    if (!requireNamespace("commonmark", quietly = TRUE)) {
      message("HTML report skipped: run install.packages(\"commonmark\") once to turn it on.")
    } else {
      esc_html <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
      
      # Images are built into the page as base64 text (no extra packages needed).
      b64_file <- function(path) {
        bytes <- as.integer(readBin(path, "raw", file.info(path)$size))
        pad   <- (3 - length(bytes) %% 3) %% 3
        m     <- matrix(c(bytes, rep(0L, pad)), nrow = 3)
        x     <- m[1, ] * 65536 + m[2, ] * 256 + m[3, ]
        idx   <- rbind(x %/% 262144, (x %/% 4096) %% 64, (x %/% 64) %% 64, x %% 64) + 1
        out   <- c(LETTERS, letters, 0:9, "+", "/")[idx]
        if (pad > 0) out[(length(out) - pad + 1):length(out)] <- "="
        paste(out, collapse = "")
      }
      img_mime <- function(f) {
        switch(tolower(tools::file_ext(f)), png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg",
               svg = "image/svg+xml", NA_character_)
      }
      # Obsidian-only syntax -> plain Markdown/HTML.
      obsidian_to_md <- function(txt) {
        embeds <- unique(regmatches(txt, gregexpr("!\\[\\[[^\\]]+\\]\\]", txt, perl = TRUE))[[1]])
        for (e in embeds) {
          f    <- sub("\\|.*$", "", sub("^!\\[\\[(.*)\\]\\]$", "\\1", e, perl = TRUE), perl = TRUE)
          path <- file.path(plot_dir, f)
          mime <- img_mime(f)
          new  <- if (file.exists(path) && !is.na(mime)) {
            paste0("![", f, "](data:", mime, ";base64,", b64_file(path), ")")
          } else {
            f
          }
          txt <- gsub(e, new, txt, fixed = TRUE)
        }
        txt <- gsub("\\[\\[[^\\]|\\\\]+\\\\?\\|([^\\]]+)\\]\\]", "\\1", txt, perl = TRUE)   # [[note|text]] -> text
        gsub("\\[\\[([^\\]]+)\\]\\]", "\\1", txt, perl = TRUE)                            # [[note]] -> note
      }
      
      # Frontmatter -> a "Properties" table at the top.
      fm_end <- which(note_lines == "---")[2]
      fm     <- note_lines[2:(fm_end - 1)]
      props  <- list()
      for (ln in fm) {
        if (grepl("^\\s+-\\s*", ln) && length(props) > 0) {
          k <- names(props)[length(props)]
          props[[k]] <- c(props[[k]], sub("^\\s+-\\s*", "", ln))
        } else if (grepl("^[^:]+:", ln)) {
          props[[sub(":.*$", "", ln)]] <- trimws(sub("^[^:]+:", "", ln))
        }
      }
      clean_prop <- function(v) {
        v <- v[nzchar(v)]
        v <- gsub('^"|"$', "", v)
        v <- gsub("\\[\\[[^\\]|]+\\|([^\\]]+)\\]\\]", "\\1", v, perl = TRUE)
        v <- gsub("\\[\\[([^\\]]+)\\]\\]", "\\1", v, perl = TRUE)
        paste(esc_html(v), collapse = "<br>")
      }
      props_html <- paste0(
        '<h2>Properties</h2><table class="props">',
        paste0("<tr><th>", esc_html(gsub("_", " ", names(props))), "</th><td>",
               vapply(props, clean_prop, ""), "</td></tr>", collapse = ""),
        "</table>")
      
      body_md   <- obsidian_to_md(paste(note_lines[(fm_end + 1):length(note_lines)], collapse = "\n"))
      body_html <- commonmark::markdown_html(body_md, extensions = TRUE)
      body_html <- gsub("\\^([^\\^\\s<]+)\\^", "<sup>\\1</sup>", body_html, perl = TRUE)     # ^2^ -> superscript
      
      html <- paste0(
        '<!DOCTYPE html><html><head><meta charset="utf-8"><title>', esc_html(base), '</title><style>',
        r"(
body { font-family: Arial, Helvetica, sans-serif; color: #1f1f1f; background: #ffffff; margin: 0; }
#nav { position: sticky; top: 0; background: #f4f4f2; border-bottom: 1px solid #d8d8d4; padding: 8px 16px;
       display: flex; gap: 8px; align-items: center; flex-wrap: wrap; z-index: 2; }
#nav .title { font-weight: bold; margin-right: auto; font-size: 14px; }
#nav button, #nav select { font-size: 13px; padding: 4px 10px; border: 1px solid #b9b9b4; border-radius: 4px; background: #fff; cursor: pointer; }
#doc { max-width: 1050px; margin: 0 auto; padding: 8px 16px 40px; line-height: 1.45; font-size: 14px; }
h1 { font-size: 22px; border-bottom: 2px solid #d8d8d4; padding-bottom: 4px; margin-top: 18px; }
h2 { font-size: 18px; margin-top: 22px; } h3 { font-size: 15px; }
table { border-collapse: collapse; margin: 8px 0 14px; font-size: 13px; display: block; overflow-x: auto; max-width: 100%; }
th, td { border: 1px solid #d8d8d4; padding: 4px 8px; text-align: left; vertical-align: top; }
th { background: #f4f4f2; }
table.props th { width: 190px; }
img { max-width: 100%; height: auto; display: block; margin: 10px 0; }
code { background: #f4f4f2; padding: 1px 4px; border-radius: 3px; }
)",
'</style></head><body><div id="nav"><span class="title">', esc_html(base), '</span>',
'<button id="prev">&#9664; Previous</button><select id="sec"></select>',
'<button id="next">Next &#9654;</button><button id="all">Show all</button></div>',
'<div id="doc">', props_html, body_html, '</div><script>',
r"(
(function () {
  var doc = document.getElementById('doc'), kids = Array.prototype.slice.call(doc.children), secs = [], cur = null;
  kids.forEach(function (el) {
    if (el.tagName === 'H1' || cur === null) { cur = document.createElement('section'); secs.push(cur); }
    cur.appendChild(el);
  });
  secs.forEach(function (s) { doc.appendChild(s); });
  var sel = document.getElementById('sec');
  secs.forEach(function (s, k) {
    var h = s.querySelector('h1'), o = document.createElement('option');
    o.value = k; o.textContent = (k + 1) + '. ' + (h ? h.textContent : 'Properties'); sel.appendChild(o);
  });
  var i = 0, all = false;
  function show() {
    secs.forEach(function (s, k) { s.style.display = (all || k === i) ? '' : 'none'; });
    sel.value = i;
    document.getElementById('all').textContent = all ? 'One section at a time' : 'Show all';
    window.scrollTo(0, 0);
  }
  document.getElementById('prev').onclick = function () { if (i > 0) { i--; all = false; show(); } };
  document.getElementById('next').onclick = function () { if (i < secs.length - 1) { i++; all = false; show(); } };
  sel.onchange = function () { i = +sel.value; all = false; show(); };
  document.getElementById('all').onclick = function () { all = !all; show(); };
  document.addEventListener('keydown', function (e) {
    if (e.key === 'ArrowRight') document.getElementById('next').click();
    if (e.key === 'ArrowLeft') document.getElementById('prev').click();
  });
  show();
})();
)",
'</script></body></html>')
      
      # Written only to R's temporary folder (cleared when R closes) -- nothing is
      # saved to the vault or opened in a web browser. RStudio's Viewer shows it.
      html_path <- file.path(tempdir(), paste0(base, ".html"))
      writeLines(html, html_path, useBytes = TRUE)
      if (interactive() && has_rstudio) {
        rstudioapi::viewer(html_path)
        cat("Report shown in RStudio's Viewer (Viewer tab; the arrow button there opens it larger).\n")
      } else {
        message("Report not shown: RStudio's Viewer is only available when running in RStudio.")
      }
    }
  }
}

# 12. Console output -----------------------------------------------------------
# Sanity check: did the right data arrive?
cat("\nMice:", n_distinct(df_long$Mouse),
    "| Slices:", n_distinct(df_long$Mouse, df_long$Slice),
    "| Cells:", n_distinct(df_long$Mouse, df_long$Slice, df_long$Cell), "\n")
print(count(df_long, Condition, name = "n_observations"))

show_table(decisions_tbl, "ANALYSIS DECISIONS (the model behind every result below)")
cat("\nIn words: ", decision_summary, "\n", sep = "")
if (length(run_choices) > 0) cat("\nChosen when running:\n", paste0("  - ", run_choices, "\n"), sep = "")
show_table(models_tbl, "The six possible models (decisions 3 x 4)")

if (!is.null(slice_sd_tbl)) {
  show_table(slice_sd_tbl, "Decision 3: slice-to-slice SD under each version of the model")
  cat(strwrap(slice_reading), sep = "\n")
  if (!is.null(slice_test)) show_table(slice_test, "Decision 3 test (REML likelihood-ratio)")
  cat(strwrap(gsub("**", "", slice_choice, fixed = TRUE)), sep = "\n")
}
if (!is.null(cell_sd_tbl)) {
  show_table(cell_sd_tbl, "Decision 4: cell-to-cell SD under each version of the model")
  cat(strwrap(cell_reading), sep = "\n")
  if (!is.null(variance_test)) show_table(variance_test, "Decision 4 test (REML likelihood-ratio)")
  cat(strwrap(gsub("**", "", cell_choice, fixed = TRUE)), sep = "\n")
}
if (!is.null(slice_sd_tbl) || !is.null(cell_sd_tbl)) {
  cat("\n"); cat(strwrap(sd_note), sep = "\n")
  cat("\n"); cat(strwrap(lr_glossary), sep = "\n")
}
cat("\nModel used: ", model_letter, " -- ", model_desc, "\n", sep = "")
cat("\nModel-fitting checks:\n")
cat(paste0("  ", fit_check_lines), sep = "\n")
if (!is.null(df_note)) cat(df_note, "\n")

# Full model output (random-effect variances, fixed effects, fit warnings).
print(summary(model))

show_table(normality_tbl, paste0("Normality ", scale_note))

if (show_diagnostics) {
  show_table(resid_normality_tbl,
             paste0("Normality of residuals ", scale_note, " -- model's assumption"))
  show_table(condition_normality_tbl,
             paste0("Normality by condition ", scale_note, " -- matches Prism's per-column check"))
  
  if (show_plots) {
    grid::grid.newpage()
    grid::grid.draw(residual_diagnostics)   # view 1: one self-contained call
    print(condition_qq)                     # view 2: one self-contained call
  }
}

# Print everything as readable tables
show_table(estimates_tbl,
           paste0("Model-estimated ", ifelse(log_transform, "geometric mean", "mean"),
                  " per condition (95% CI) -- describes each condition, NOT used for comparisons"))

if (!is.null(comparisons_tbl)) {
  show_table(comparisons_tbl,
             paste0("Comparisons (", ifelse(log_transform, "ratios", "differences"),
                    ", 95% CI; p and CI adjusted: ", adjust_label, ")"))
}
if (!is.null(df_breakdown_tbl)) {
  show_table(df_breakdown_tbl, "Where the degrees of freedom come from (% of Var(SE^2) by source)")
  cat(strwrap(df_explain), sep = "\n")
}

show_table(descriptives_tbl, "Descriptive statistics (raw data)")
cat("\nCell outliers: ", outlier_note, "\n", sep = "")
if (!is.null(outlier_tbl)) show_table(outlier_tbl, "ROUT cell outliers")
cat("\nSlice outliers: ", slice_outlier_note, "\n", sep = "")
if (!is.null(slice_outlier_tbl)) show_table(slice_outlier_tbl, "ROUT slice outliers (flagged only)")

if (length(pub_save_errors) > 0) {
  cat("\nPUBLICATION FIGURES NOT SAVED:\n")
  cat(paste0("  ", pub_save_errors), sep = "\n")
}

cat("\nWRITING TEMPLATE (draft -- edit before use)\n===========================================\n")
cat("\nMethods:\n");  cat(strwrap(methods_text), sep = "\n")
cat("\nResults:\n");  cat(strwrap(means_text), sep = "\n")
for (rt in results_text) cat(strwrap(paste("-", rt), exdent = 2), sep = "\n")
cat("\nFigure legend:\n"); cat(strwrap(legend_text), sep = "\n")
cat("\n"); cat(caution_text, sep = "\n")
