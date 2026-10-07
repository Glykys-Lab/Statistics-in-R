# LabStats: local R analysis for Prism users

A local Shiny GUI plus reusable R functions for common laboratory analyses. This is an initial toolkit, not a complete GraphPad Prism replacement. Data stay in the local R process; no cloud service is required to analyze them.

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

1. Choose **Column** or **XY** in the project navigator. To import a file, click **Import CSV / Excel** in the top toolbar, then **Browse** in the Data table tab. CSV, TSV, XLSX, and XLS are supported (25 MB maximum). The table loads automatically. For semicolon-separated CSV files, change **Text delimiter**; for Excel, select the worksheet. Headers must be unique and nonempty. Only data values are read; formulas are not recalculated by R. You can also load an example, create a blank table, or paste a spreadsheet with headers.
2. Edit cells directly in the **Data table** sheet. Tab moves across cells; Enter or up/down moves between rows. Blanks represent missing observations. Up to 100 rows are editable; every imported row is analyzed. For Column tables, select measurement columns and an optional subject ID, or switch to long layout.
3. For Column tables, choose two versus several groups, independent versus matched measurements, Gaussian versus rank-based assumptions, and (for independent Gaussian data) equal versus unequal variance. For repeated Gaussian data, choose complete-subject ANOVA or a mixed model retaining incomplete subjects. XY offers simple linear regression, Pearson correlation, and Spearman correlation. Review the analysis description and parameters, then click **Analyze**. Analyze is disabled when the current data and analysis settings match the selected saved result; changing graph formatting alone does not require a new analysis.
4. Read the tabbed data audit, warnings, distribution tests, and result sheets. **Normality & lognormality** and **Identify outliers (ROUT)** are also standalone analysis families. For matched repeated data, choose Gaussian, rank-based, or lognormal assumptions. Shapiro-Wilk and D'Agostino-Pearson p-values do not automatically choose a test. For paired comparisons, normality concerns differences; for regression, it concerns residuals.
5. Each analysis saves its input snapshot and creates available data, normal QQ, lognormal QQ, and comparison-interval graphs. Open these through **Data Tables / Results / Graphs / Layouts** in the searchable navigator. **Family** shows sheets belonging to the selected table. Changing current data/settings requires a new analysis; selecting a saved result displays its original snapshot.
6. In **Graph**, select a graph type and format symbols, colors, matched lines, axes, error bars, and the legend. **Save graph** preserves those settings. **Layouts** combines saved graphs as separate panels or overlays compatible Column/XY data graphs. Column overlays show observations, matched lines where enabled, means, and 95% mean CIs; use compatible response units.
7. Save notes and download the complete JSON project in **Project info**. Reopen it through **Open a LabStats project**; all saved analyses are recalculated from their input snapshots. Export tables, reports, PNG/PDF figures or layouts, a reproducible R script of the selected graph, and session details.

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

For independent groups, subject IDs are optional and ignored. For paired tests, each ID identifies the same experimental unit across conditions. Group names are sorted alphabetically; the t-test estimate and its CI represent the first minus the second group, with the contrast printed in the result.

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
| Welch one-way ANOVA | Pairwise Welch tests; Sidak, Holm, Bonferroni, BH/FDR, or no correction |
| Kruskal-Wallis | Pairwise Wilcoxon tests with chosen correction; not Dunn tests |
| Friedman | Complete repeated-measure blocks; omnibus test only |
| Repeated-measures one-way ANOVA | Complete subjects; ordinary and GG-corrected F tests, GG selected by default; matching-effectiveness F test; pairwise paired t tests when sphericity is not assumed, pooled error when assumed |
| Mixed-effects | `nlme::lme`, REML, fixed condition effect + subject random intercept, common residual variance; retains incomplete subjects; GG-corrected marginal F test, comparisons, fixed-effect CIs, variance components and matching-effectiveness likelihood-ratio test |
| Lognormal repeated models | Repeated ANOVA or REML mixed model on natural logs; all observed values must be positive; comparisons report geometric-mean ratios and ratio CIs; raw graphs remain in original units |
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

All hypothesis tests are two-sided where applicable. All-pairs or control contrasts form the chosen correction family. Sidak uses `1-(1-p)^m` and individual confidence level `confidence^(1/m)` for a family of m comparisons; Bonferroni uses alpha/m. Sidak familywise interpretation assumes independence and is approximate for correlated contrasts; Bonferroni is valid under arbitrary dependence. Holm/BH adjust p-values but their displayed t intervals remain explicitly **unadjusted**. `BH` controls false discovery rate rather than familywise error. Table stars use conventional thresholds (0.05, 0.01, 0.001, 0.0001); the separate decision column uses the chosen alpha. Group-mean graph intervals are unadjusted. For paired data they describe individual condition means rather than the within-subject difference. GUI p-values use scientific notation or `<1e-15`; CSV exports retain numerical precision.

For mixed-model matching, both the conventional chi-square(1) likelihood-ratio p-value and a boundary-mixture p-value are reported because subject variance is bounded at zero. GG uses sample covariance for complete subjects. With missing data it estimates unstructured covariance using observed-data Gaussian EM, including conditional covariance contributions; it does not replace missing outcomes in the input. Every condition pair needs at least two overlapping subjects, covariance must be estimable, and EM must converge. Choose the explicitly uncorrected sphericity setting if the GG estimator is unavailable. The estimator is disclosed in the result sheet; exact Prism equivalence is unverified.

Column ROUT follows the published adaptive Lorentzian constant fit, robust residual scale `quantile(abs(residual), .6827)*N/(N-1)`, and t-tail FDR scan with selected Q. Q applies separately to each column. Too-small, constant, degenerate-scale, or unconverged columns are explicitly marked unavailable. Its numerical behavior has been checked with synthetic data and invariance tests; exact Prism parity has not been validated. Flagging leaves the source table intact. **Create cleaned data table** uses the saved analysis snapshot, keeping rows and IDs while blanking flagged measurements.

The 4PL model is `bottom + span / (1 + exp(hill * (log(dose) - log(EC50)))))`, with positive span. Positive Hill means decreasing; negative means increasing. The Hill magnitude is bounded to [0.01, 10]; log-EC50 is bounded to the observed log-dose range plus/minus 10. The fit tries three initial midpoints and selects the converged fit with least residual sum of squares. EC50 intervals use a t-based local Wald interval on log-EC50 transformed back to concentrations. These are not profile-likelihood intervals. Plateaus must be supported by the dose range. A converged fit alone is insufficient evidence that parameters are reliable.

## Limitations and migration checks

- Native Prism `.pzfx`/`.prism` projects are not imported. Export raw tables to CSV or Excel.
- Nested designs, multi-factor repeated/mixed designs, Dunnett/Dunn tests, multiple t-test families, curve comparisons, weighting, constrained plateaus, other nonlinear models, and power analysis are not yet implemented. The mixed model uses a single subject random intercept, not an arbitrary covariance structure. Its inference and the repeated-ANOVA paired post-hoc tests may differ from Prism.
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

`R/analysis.R` holds the analysis functions; `R/plots.R` holds figures; `R/data.R` handles examples, editing, layout conversion, and script export. `R/project.R` handles JSON persistence and transformations; `R/workspace.inc` connects project navigation to the server. `R/ui.R` supplies the interface and `app.R` the Shiny server.

## Verification

From the project directory:

```r
source("tests/run_tests.R")
```

Checks compare results against R references, verify subject matching and invalid-input rejection, check complete/missing-data GG and REML behavior, verify published normality transformations, adjusted comparison CIs and ratio models, ROUT synthetic outliers, baseline transformations, project round-trips, all graph types, overlays and layouts. Exported scripts are executed independently. Server checks exercise editing, stale/saved results, missing-data models, XY/paste, normalization, lognormal analyses, saved graphs, layouts and reopening projects. Test outputs go to ignored `artifacts/`.

Official references: [Shiny app execution](https://shiny.posit.co/r/reference/shiny/latest/runapp.html), [R t tests](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html), [ANOVA](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/aov.html), [nonlinear least squares](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/nls.html).

Repeated-measures references: [GraphPad sphericity and Geisser-Greenhouse correction](https://www.graphpad.com/guides/prism/latest/statistics/stat_sphericity_and_compound_symmet.htm), [nlme package reference](https://stat.ethz.ch/CRAN/web/packages/nlme/nlme.pdf). These document the methods; this implementation has not been validated for exact equivalence to Prism.

Additional primary references: [GraphPad mixed-model interpretation](https://www.graphpad.com/guides/prism/latest/statistics/stat_intepreting-results-mixed-mode.htm), [GraphPad repeated-measures comparisons](https://www.graphpad.com/guides/prism/latest/statistics/stat_multiple_comparisons_after_rep.htm), [Motulsky & Brown ROUT paper](https://link.springer.com/article/10.1186/1471-2105-7-123), [SciPy normality transformations](https://github.com/scipy/scipy/blob/main/scipy/stats/_stats_py.py).
