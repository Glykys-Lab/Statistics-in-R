# LabStats: local R analysis for Prism users

A local Shiny GUI plus reusable R functions for common laboratory analyses. This is an initial toolkit, not a complete GraphPad Prism replacement. Data stay in the local R process; no cloud service is required to analyze them.

For an explanation of each source file and how the modules connect, see the [file guide](FILE_GUIDE.md).

## Start on this Windows computer

Double-click **Start LabStats.cmd**. It checks/installs packages, then opens the app at `http://127.0.0.1:7319`. Keep its console window open while working; close it to stop the app. Internet is needed for the first package installation, not for subsequent analysis.

Alternatively, open **LabStats.Rproj** in RStudio and run:

```r
source("setup.R")  # first use on each computer
source("launch.R")
```

R 4.1 or newer is required; development and testing used R 4.6.1. The launcher uses port 7319; change `port` in `launch.R` if another service uses it. In RStudio you can also open `app.R` and click **Run App**. Dependencies live in `.R-library` inside the project.

To give this to another lab member, copy the source folder (omit `.R-library`, `.Rproj.user`, and `artifacts`), install R and optionally RStudio on their computer, and run `setup.R`. Packages must be installed for each operating system and R version. `setup.R` uses current compatible CRAN packages; it is not a version lockfile. Save session details with every finalized analysis.

## First analysis

1. Choose **Column** or **XY** in the project navigator. Click **Import CSV / Excel**, then **Browse** in the Data table tab. Select one or more CSV, TSV, XLSX, or XLS files (25 MB total per selection). When selecting several files, use **File to analyze** to choose one at a time; Analyze applies only to the selected file. Each file keeps its own edits, saved table, and analysis results as you switch. A new Browse selection starts a fresh upload selection; saved tables and results remain in the project navigator. For semicolon-separated CSV files, change **Text delimiter**; for Excel, select the worksheet. Headers must be unique and nonempty. Only data values are read; formulas are not recalculated by R. You can also load an example, create a blank table, or paste a spreadsheet with headers.
2. Edit cells directly in the **Data table** sheet. Tab moves across cells; Enter or up/down moves between rows. Blanks represent missing observations. Up to 100 rows are editable; every imported row is analyzed. For Column tables, select measurement columns and an optional subject ID, or switch to long layout.
3. For Column tables, choose two versus several groups in the navigator, then open **Analysis options**. Its five tabs cover experimental design, repeated measures, comparisons, reporting and residuals. Choose independent versus matched measurements, Gaussian/lognormal/rank assumptions, and (for independent Gaussian data) equal versus unequal variance. For repeated Gaussian/lognormal data, choose complete-subject ANOVA, a mixed model retaining incomplete subjects, or **It depends** to select ANOVA for complete data and mixed effects for missing data. XY offers simple linear regression, Pearson correlation, and Spearman correlation. Review the analysis description and parameters, then click **Analyze**. Analyze is disabled when data and calculation settings match the selected saved result; graph formatting and display digits do not require a new analysis.
4. Results use spreadsheet-style **Tabular results**, **Multiple comparisons**, **Descriptive statistics**, **Predicted means** and **Normality & lognormality** sheets where applicable. Dataset columns follow worksheet order; statistics appear as rows. Calculation notes and renaming are in a collapsible section; data audits, coefficients and additional diagnostics are in **Additional results**. Descriptives include all entered values, missing counts, minimum/maximum, quartiles, mean, SD, SEM and mean CIs; complete-case matched models may analyze fewer subjects. Comparison details report the observed or model-predicted means actually used for the contrast, n1/n2, SE, t and df. Overall-model predicted means can differ from pairwise-model means. **Normality & lognormality** and **Identify outliers (ROUT)** remain standalone analysis families; ROUT has **Summary**, **Outliers** and **Cleaned data** sheets preserving original subject rows. Shapiro-Wilk and D'Agostino-Pearson p-values do not automatically choose a test. For paired comparisons, normality concerns differences; for regression, it concerns residuals.
5. Each analysis saves its input snapshot and creates available data, normal QQ and lognormal QQ graphs. One-way options control creation of comparison intervals, ranks, matched differences and residual graphs. Open these through **Data Tables / Results / Graphs / Layouts** in the searchable navigator. **Family** shows sheets belonging to the selected table. Changing current data/settings requires a new analysis; selecting a saved result displays its original snapshot.
6. In **Graph**, select a graph type and format symbols, colors, matched lines, axes, error bars, and the legend. **Save graph** preserves those settings. **Layouts** combines saved graphs as separate panels or overlays compatible Column/XY data graphs. Column overlays show observations, matched lines where enabled, means, and 95% mean CIs; use compatible response units.
7. Save notes and download the complete JSON project in **Project info**. Reopen it through **Open a LabStats project**; all saved analyses are recalculated from their input snapshots. **Export & R code** offers the arranged result sheets as displayed, plus **Raw:** tables retaining numerical precision. Reports and standalone scripts include the arranged sheets; scripts also export raw tables. Export PNG/PDF figures or layouts and session details as before.

Example datasets are synthetic, not experimental evidence. A generated R script embeds the imported data; handle it as a copy of that dataset.

## Data layouts

Long layout (available for imported datasets):

```csv
subject,group,response
S01,Control,9.1
S02,Control,10.4
S01,Treatment,12.1
S02,Treatment,14.2
```

For independent groups, subject IDs are optional and ignored. For paired tests, each ID identifies the same experimental unit across conditions. Groups follow the worksheet's measurement-column order (or first appearance for long-format tables), keeping experimental steps in sequence on graphs and legends. The t-test estimate and its CI represent the first minus the second group, with the contrast printed in the result. Explicit factor levels are preserved for R inputs.

Prism-like wide layout:

```csv
subject,Control,Treatment
S01,9.1,12.1
S02,10.4,14.2
```

Choose **Column: one group per column**, select only the measurement columns, and select the ID column. If you use row numbers as IDs, matching assumes values on the same row belong together. Blank cells are missing measurements, not zeros. Wide conversion produces `subject`, `group`, and `response` columns. **Export analysis CSV** exports this mapped data; the reproducible script also embeds the mapped input.

XY and dose-response: one numeric X/concentration column and one response column, with replicates on separate rows. Dose-response supports **one curve per analysis**, so import the selected curve as its own table. Use actual positive concentrations, not log10 concentrations. Both increasing and decreasing curves are supported; zero-dose controls are rejected explicitly.

Contingency: columns for row category, column category, and integer count; duplicate cells are summed. Omitted combinations count as zero. Survival: time, event (`1` event / `0` censored), and group, with one row per subject. Examples in `examples/` demonstrate these layouts.

## Supported methods

| Workflow | Implementation and key settings |
|---|---|
| Descriptives | n, mean, SD, SEM, median, quartiles, t-based CI of each mean |
| Independent t tests | Welch or explicitly selected equal-variance Student test |
| Paired t test | Subject-ID matching; complete pairs only |
| Nonparametric two-group tests | Mann-Whitney or paired Wilcoxon; normal approximation with continuity correction |
| One-way ANOVA | Ordinary ANOVA; Tukey all-pairs or pooled-error t comparisons using the selected correction |
| Welch one-way ANOVA | Pairwise Welch tests with the selected correction |
| Kruskal-Wallis | Optional tie-corrected Dunn rank comparisons or pairwise Wilcoxon tests |
| Friedman | Complete repeated-measure blocks; optional tie-corrected Dunn comparisons of within-subject ranks or paired Wilcoxon tests |
| Repeated-measures one-way ANOVA | Complete subjects; ordinary and GG-corrected F tests, GG selected by default; matching-effectiveness F test; pairwise paired t tests when sphericity is not assumed, pooled error when assumed |
| Mixed-effects | REML condition effects with compound symmetry; `nlme::lme` for a nonnegative random intercept, `nlme::gls` for signed covariance or removal of a nonpositive subject component; retains incomplete subjects; optional GG correction, comparisons, fixed-effect CIs and matching diagnostic |
| Automatic repeated analysis | Complete data select repeated-measures ANOVA; any missing selected measurement or missing subject-condition record selects mixed effects; Gaussian or lognormal scale |
| Lognormal repeated models | Repeated ANOVA or REML mixed model on natural logs; all observed values must be positive; comparisons report geometric-mean ratios and ratio CIs; raw graphs remain in original units |
| Nested repeated cells | Mouse / Slice / Cell random effects, optional slice × condition effects (none, shared or per-condition SD); shared or per-condition residual noise; Satterthwaite inference |
| Nested between groups | Mouse and within-mouse slice random effects; shared or separate slice SDs and residual SDs per condition; Satterthwaite inference |
| Normality / lognormality | Per-column Shapiro-Wilk (3–5000 values) and D'Agostino-Pearson (at least 8); lognormal tests operate on logs and require positive values; normal and lognormal actual-vs-predicted QQ graphs; counts of zeros/negative values |
| Assumption diagnostics | Shapiro-Wilk and QQ plots of within-subject differences for matched analyses, residuals for regression, and group values for independent comparisons |
| Outlier identification | Column ROUT constant fit with selectable Q; summary, flagged values, all residual details, and cleaned-data sheets; explicit action creates a separate cleaned table |
| Baseline transformations | Row-wise ratio, percent of baseline, change, percent change, or natural log; creates a derived table preserving rows and IDs; undefined divisions/nonpositive logs are blanked and counted |
| Two-way ANOVA | Independent, balanced, fully crossed design with interaction; at least two replicates per cell |
| Correlation | Pearson or Spearman; Spearman uses an asymptotic p-value |
| Linear regression | Intercept + slope, coefficient CIs, R-squared, mean-response confidence band, residual plot |
| 4PL dose-response | Unweighted nonlinear least squares, free plateaus, EC50, Hill slope, approximate EC50 interval, residual plot |
| Contingency | Pearson chi-square without Yates correction, or Fisher exact |
| Survival | Kaplan-Meier, censor marks, confidence intervals and event-time risk counts in the table, log-rank test |

All hypothesis tests are two-sided where applicable. All-pairs, control or selected contrasts form the chosen correction family; **None** omits follow-ups. Sidak uses `1-(1-p)^m` and individual confidence level `confidence^(1/m)` for a family of m comparisons; Bonferroni uses alpha/m. Sidak familywise interpretation assumes independence and is approximate for correlated contrasts; Bonferroni is valid under arbitrary dependence. Holm, Holm-Sidak and FDR methods adjust p-values but leave displayed t intervals explicitly **unadjusted**. FDR choices are two-stage Benjamini-Krieger-Yekutieli (BKY), Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY), with a separate target Q and discovery decision. BKY implements the two-stage linear step-up procedure, not an iterated variant. Stars use conventional thresholds (0.05, 0.01, 0.001, 0.0001); decisions use selected alpha or Q. Group-mean graph intervals are unadjusted and describe individual condition means for paired data. Output digits and p-value styles affect screens/reports; CSV exports retain numerical precision.

For nonnegative mixed-model matching, both the conventional chi-square(1) likelihood-ratio p-value and a boundary-mixture p-value are reported. The signed model reports only the conventional test. Its signed covariance component must leave the observation covariance positive definite; a negative component has no physical variance/SD interpretation. GLS condition tests use repeated-measures containment degrees of freedom, not Satterthwaite. Removing a nonpositive subject component refits independent residuals. GLM-based initial values and a small initial value are both tried when enabled; the converged fit with greatest likelihood is selected. GG uses sample covariance for complete subjects. With missing data it estimates unstructured covariance using observed-data Gaussian EM, including conditional covariance contributions; it does not replace missing outcomes in the input. Every condition pair needs at least two overlapping subjects, covariance must be estimable, and EM must converge. Choose uncorrected sphericity if the GG estimator is unavailable. The estimator is disclosed in the result sheet; exact Prism equivalence is unverified.

## One-way options

- **Multiple Comparisons:** all pairs, a control column, selected named pairs, no comparisons, or linear trend across equally spaced column positions in worksheet order. Repeated-data trend requires sphericity assumed. Swap direction reverses estimates/intervals without changing two-sided p-values.
- **Options:** choose family-wise correction, FDR or no correction; alpha/confidence, graph creation, descriptive statistics, goodness of fit, effect sizes and AICc model comparisons. Gaussian ANOVA effect sizes are eta-squared/partial eta-squared; mixed effects report treatment differences standardized by residual SD. Variance-based mixed R-squared is unavailable for negative covariance components. AICc compares condition means to a common mean (retaining subject terms in repeated ANOVA), refitting with **ML** when fixed effects differ. Its sample size is the number of observed outcomes; clustered-data small-sample AICc is an approximation. AICc is unavailable when observations do not exceed parameter count plus one.
- **Residuals:** residual-versus-fitted, absolute-residual, normal QQ and heatmap graphs; pooled D'Agostino-Pearson, Anderson-Darling, Shapiro-Wilk and Lilliefors-corrected Kolmogorov-Smirnov tests. Independent groups also allow Brown-Forsythe and Bartlett variance tests. Invalid sample sizes or constant residuals are explicitly marked unavailable. These are exploratory diagnostics because residuals share fitted parameters.
- **Defaults:** save/reset analysis choices in browser local storage. Data, column mappings, control columns and selected pairs are excluded. Saved project analyses keep their own calculation settings; browser defaults configure future analyses.

Column ROUT follows the published adaptive Lorentzian constant fit, robust residual scale `quantile(abs(residual), .6827)*N/(N-1)`, and t-tail FDR scan with selected Q. Q applies separately to each column. Too-small, constant, degenerate-scale, or unconverged columns are explicitly marked unavailable. Its numerical behavior has been checked with synthetic data and invariance tests; exact Prism parity has not been validated. Flagging leaves the source table intact. **Create cleaned data table** uses the saved analysis snapshot, keeping rows and IDs while blanking flagged measurements.

The 4PL model is `bottom + span / (1 + exp(hill * (log(dose) - log(EC50)))))`, with positive span. Positive Hill means decreasing; negative means increasing. The Hill magnitude is bounded to [0.01, 10]; log-EC50 is bounded to the observed log-dose range plus/minus 10. The fit tries three initial midpoints and selects the converged fit with least residual sum of squares. EC50 intervals use a t-based local Wald interval on log-EC50 transformed back to concentrations. These are not profile-likelihood intervals. Plateaus must be supported by the dose range. A converged fit alone is insufficient evidence that parameters are reliable.

## Nested mouse / slice / cell workflow

Choose **Nested models** under Table type, then **Analysis options**. The two included example tables demonstrate repeated and between-group designs. For your data, import or paste a table and map the Mouse, Slice and Cell ID columns in the navigator.

- **Repeated:** one row per Mouse / Slice / Cell combination and one numeric column per condition. Select the measurement columns; their worksheet order is the default sequence. Missing measurements are omitted individually.
- **Between groups:** one row per cell with ID, Condition and numeric measurement columns. Each cell occurs once and each slice receives one condition. Reused slice/cell labels are scoped by their parent identifiers.
- **Design & fitting:** raw or log10 responses; candidate slice structures and shared/unequal residual noise; explicit or automatic engine selection; per-comparison Satterthwaite df or the minimum within each set of means/contrasts. The supplied scripts combine nlme estimates with df from the corresponding lmer model. Candidate-model LRTs follow the reference scripts' conventional approximations; there is no preliminary normality-test gate.
- **Comparisons:** all pairs, a selected control, consecutive conditions or selected pairs. Defaults are Tukey, the emmeans Dunnett approximation (`dunnettx`), or Holm as appropriate. Holm p-values use Bonferroni intervals in this path. Intervals are 95%. Optional mean-against-zero tests are unadjusted; on the log10 scale their original-scale null is a geometric mean of **1**.
- **Outliers:** ROUT at the chosen Q on raw values or nested residuals; cells, slices or both. Flagging keeps observations. Selecting removal excludes flagged cell measurements from this analysis while preserving the source worksheet; slice outliers are always kept. IDs and actions are reported.
- **Figures:** model means and CIs or descriptive medians/IQR; mouse-coloured cells, slice means, matching lines, significance brackets and comparison intervals. Publication controls specify dimensions in millimetres, fonts/strokes in points, colours, labels, micro units, ion superscripts, key placement, beeswarm/jitter, CI/SEM and summary placement. Use **Publication SuperPlot** in the Graph tab and download the ZIP for final-sized SVG/PDF/EPS files. A separate-panel layout can include nested figures; overlays are rejected to preserve their hierarchy and model summaries. EPS font availability depends on R's PostScript font database; any export warnings are reported.
- **Results & export:** decisions, selected/candidate models, variance tests and components, fixed effects, model means/CIs, df sources, mouse/slice/cell counts, missing values, diagnostics, outliers and draft Methods/Results/legend text. Tables keep numeric precision in raw CSV exports. The full report includes section navigation and embedded images; download its HTML, an Obsidian Markdown note, or a ZIP containing notes, tables and all figures. Exported R code embeds both reference scripts and runs outside this project; in interactive RStudio it also opens the HTML report in the Viewer when `rstudioapi` is available.

The model calculations come from the supplied `linear-mixed-model-repeated.R` and `linear-mixed-model-between-groups.R`, bundled as `R/nested-reference-*.inc`. LabStats overrides paths before file access, disables desktop prompts and uses a separate temporary directory for each run. Upload and browser downloads replace the scripts' desktop file pickers and personal output folders. The shared normality helper provides Shapiro-Wilk and D'Agostino-Pearson diagnostics, including constant/small-sample handling. Figure labels and key names survive project JSON round-trips. Choosing a control does not reorder the conditions. Publication figures default to 80 mm wide for readable labels; the guide's 40 / 30 / 10 mm dimensions remain available in the controls.

These models are compatible with the existing analyses but use a different experimental design. The single-subject path retains natural-log / Geisser-Greenhouse options and compound symmetry; nested models use log10 / Satterthwaite and separate mouse/slice/cell effects. Do not expect identical uncertainties or p-values from those different models. Conventional Holm intervals also differ by path: the original one-way path labels them unadjusted, while the supplied nested scripts use Bonferroni intervals. Package references: [lme4 model specification](https://lme4.github.io/lme4/reference/lmer.html), [emmeans model df](https://github.com/rvlenth/emmeans/blob/main/vignettes/models.Rmd), [emmeans Dunnett approximation](https://github.com/rvlenth/emmeans/blob/main/R/summary.R).

Run `source('setup.R')` after updating on another computer to install the additional local dependencies. Nested settings are part of saved tables/analyses and exported projects.

## Limitations and migration checks

- Native Prism `.pzfx`/`.prism` projects are not imported. Export raw tables to CSV or Excel.
- General multi-factor repeated/mixed designs beyond the supplied nested templates, multiple t-test families, curve comparisons, weighting, constrained plateaus, other nonlinear models, and power analysis are not yet implemented. The original single-subject mixed path uses compound symmetry. Nested control comparisons use the Dunnett approximation, not exact multivariate-t Dunnett inference. Inference may differ from Prism.
- Two-way ANOVA rejects unbalanced or missing cells. It does not silently choose a sums-of-squares convention.
- Outlier cleaning and normalization are explicit actions creating new tables. There is no technical-replicate aggregation or substitution of missing outcomes. Missing selected values are audited. Matching rejects duplicate subject-condition records. Paired, Friedman, and repeated-ANOVA tests exclude incomplete subjects; the mixed model retains their available observations and assumes ignorable missingness given the model. Repeated conditions must overlap across subjects to identify the mixed effect.
- Independence and experimental-unit selection require study knowledge. Residuals vs fitted is a diagnostic aid, not a full assumption check. The app does not select a test using preliminary normality tests.
- Reproduce a known lab analysis in both programs before adopting a workflow. Match pairing, group order, tails, equal-variance assumption, comparison family/correction, missing-value rules, nonlinear equation/constraints/weighting, and confidence-interval method. This implementation has been checked against R references, not against a licensed Prism installation.
- Sessions are in memory until you download a project JSON file. Save the project before closing R. Analyses store immutable input snapshots; there is no general worksheet undo system. Native Prism project files remain unsupported.

## Use the scripts without the GUI

```r
source("R/analysis.R")
source("R/plots.R")
dat <- read.csv("examples/independent_groups.csv")
result <- run_analysis(dat, list(
  method = "anova", y = "response", group = "group", confidence = 0.95
))
result$tables
writeLines(analysis_report(result), "report.txt")
plot <- make_plot(result)
ggplot2::ggsave("figure.pdf", plot, width = 7, height = 5)
```

After using `setup.R`, add `.libPaths(c(normalizePath('.R-library'), .libPaths()))` in a fresh R session to use this project's installed packages. Scripts exported by the GUI are self-contained and do not require the source project; they still require their R packages.

`R/analysis.R` holds the analysis engine and `R/options.R` adds comparison/model/diagnostic/reporting helpers; `R/plots.R` holds figures; `R/data.R` handles examples, editing, layout conversion, and script export. `R/project.R` handles JSON persistence and transformations; `R/workspace.inc` connects project navigation to the server. `R/ui.R` and `R/options_ui.R` supply the interface; `app.R` and `R/options_server.inc` handle reactive settings and workflows.

## Verification

From the project directory:

```r
source("tests/run_tests.R")
```

Checks compare results against R references, verify subject matching and invalid-input rejection, check complete/missing-data GG and REML behavior, verify published normality transformations, adjusted comparison CIs and ratio models, ROUT synthetic outliers, baseline transformations, project round-trips, all graph types, overlays and layouts. Exported scripts are executed independently. Server checks exercise editing, stale/saved results, missing-data models, XY/paste, normalization, lognormal analyses, saved graphs, layouts and reopening projects. Test outputs go to ignored `artifacts/`.

Official references: [Shiny app execution](https://shiny.posit.co/r/reference/shiny/latest/runapp.html), [R t tests](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html), [ANOVA](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/aov.html), [nonlinear least squares](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/nls.html).

Repeated-measures references: [GraphPad sphericity and Geisser-Greenhouse correction](https://www.graphpad.com/guides/prism/latest/statistics/stat_sphericity_and_compound_symmet.htm), [nlme package reference](https://stat.ethz.ch/CRAN/web/packages/nlme/nlme.pdf). These document the methods; this implementation has not been validated for exact equivalence to Prism.

Additional primary references: [GraphPad mixed-model interpretation](https://www.graphpad.com/guides/prism/latest/statistics/stat_intepreting-results-mixed-mode.htm), [GraphPad repeated-measures comparisons](https://www.graphpad.com/guides/prism/latest/statistics/stat_multiple_comparisons_after_rep.htm), [Motulsky & Brown ROUT paper](https://link.springer.com/article/10.1186/1471-2105-7-123), [SciPy normality transformations](https://github.com/scipy/scipy/blob/main/scipy/stats/_stats_py.py).

One-way options references: [GraphPad BKY algorithms](https://www.graphpad.com/guides/prism/latest/statistics/stat_algorithms_for_controlling_the.htm), [mixed-model starting values](https://www.graphpad.com/support/faq/updated-reml-calculation-method-920/), [residual options](https://www.graphpad.com/guides/prism/latest/statistics/stat_1wayanova_residuals-tab.htm), [trend and sphericity](https://www.graphpad.com/support/faq/test-for-trend-1-way-RM-ANOVA-no-sphericity/), [nortest normality tests](https://cran.r-project.org/web/packages/nortest/nortest.pdf).
