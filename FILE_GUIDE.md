# LabStats file guide

This guide explains the maintained files in this project and how they work together. Paths are relative to the project root. Generated package files, caches and exported artifacts are described at the end rather than treated as application source.

## How the application fits together

1. `Start LabStats.cmd` runs package setup and starts `launch.R`. You can also run those R scripts yourself from the project directory.
2. `launch.R` starts Shiny, which loads `app.R` and its UI and analysis modules.
3. `R/ui.R` builds the page. `app.R` and the server `.inc` files respond to input, manage the worksheet and save analysis snapshots.
4. `R/data.R` imports or reshapes a table. `R/analysis.R` runs the selected calculation; nested methods dispatch to `R/nested.R` and the two bundled references.
5. `R/results.R` defines result sheets, `R/results_ui.R` renders them, and `R/plots.R` creates graphs. The same results feed exports.
6. `R/project.R` stores data and settings as JSON. Opening a saved project recalculates its analyses from their saved inputs.

The `.inc` files contain R code. The server includes are sourced inside `server()` so they use that session's inputs, outputs and state. The nested reference includes contain top-level calculations and file operations; keeping them as `.inc` prevents Shiny from automatically sourcing them as ordinary files in `R/`. Only the nested adapter should execute them in the app.

## Startup files

| File | What it does | When to use or change it |
| --- | --- | --- |
| [`Start LabStats.cmd`](Start%20LabStats.cmd) | Windows launcher. Finds `Rscript.exe` on PATH or under Program Files, runs `setup.R`, stops on setup failure, and launches the app. | Double-click to start on Windows; edit when changing R discovery or launcher behavior. |
| [`setup.R`](setup.R) | Checks and installs the application and nested-model dependencies into the project's `.R-library`. Verifies they can be loaded. | Run on a new computer or after dependency changes. It installs compatible CRAN versions; it is not a version lockfile. |
| [`launch.R`](launch.R) | Starts Shiny on `127.0.0.1:7319` and opens the browser, using the project library first. | Run after setup; change the port here if necessary. |
| [`app.R`](app.R) | Main Shiny entry point. Loads modules, creates session state, imports and edits data, chooses methods, validates Analyze clicks, saves results, and binds tables, graphs and downloads. | Change application behavior that connects several modules, such as input routing or result freshness. |

## Data, statistics and saved projects

| File | What it does | Important relationships |
| --- | --- | --- |
| [`R/data.R`](R/data.R) | Provides ordinary example data; imports CSV/TSV/Excel and pasted tables; validates edits; converts wide worksheets to analysis rows; generates standalone R exports with embedded data and code. | Used by `app.R`, tests and exported scripts. Nested examples are generated in `nested.R`. Worksheet column order determines the experimental sequence. |
| [`R/analysis.R`](R/analysis.R) | Main statistical dispatcher and shared validation. Implements descriptive statistics, tests/ANOVA, repeated and mixed analyses, normality, ROUT, regression, dose response, contingency and survival analyses. Returns numerical tables, fits and audit notes. | Loads `options.R`, `results.R` and `nested.R`. Can be sourced without starting Shiny. Selected tests differ in whether they omit missing measurements or incomplete matched subjects. |
| [`R/t_tests.R`](R/t_tests.R) | Runs normal/lognormal one-sample, two-sample and multiple t-tests, validates tails and hypotheses, matches pairs, and corrects the planned family. | Called by `analysis.R`; embedded in standalone exports and displayed by `results.R`. |
| [`R/options.R`](R/options.R) | Shared one-way statistical helpers: comparison families, corrections, signed or nonnegative covariance policies, rank/trend follow-ups, residual diagnostics, optional model tables and display formatting. | Called by `analysis.R` and included in standalone exports. Its method-specific calculations do not depend on Shiny inputs. |
| [`R/project.R`](R/project.R) | Creates project records and stable IDs, encodes typed tables/settings as JSON, validates and recalculates loaded analyses, and supplies baseline transformations and cleaned-data reconstruction. | Used by `workspace.inc` and tests. Project JSON saves inputs and settings, not executable serialized R fit objects. |
| [`R/workspace.inc`](R/workspace.inc) | Session-level project operations: save/open/duplicate tables, selected-row highlighting, archived results, graph/layout records, derived tables, notes, project downloads and restoration. | Sourced inside the `app.R` server. Calls `project.R`; archived snapshots remain separate from the current editable worksheet. |

## User interface, result tables and graphs

| File | What it does | Important relationships |
| --- | --- | --- |
| [`R/ui.R`](R/ui.R) | Defines the main page: project navigator, data sources, import/edit table, analysis options, result tabs, graphs, layouts, project information and exports. | Uses the option/result UI helpers. Input IDs connect its widgets to `app.R` and the server includes. |
| [`R/options_ui.R`](R/options_ui.R) | Builds ordinary one-way/repeated-measures settings tabs and valid correction choices. Includes experimental design, comparison, reporting and residual controls. | Used by `ui.R`; nested options have their own UI. |
| [`R/options_server.inc`](R/options_server.inc) | Converts applicable ordinary settings to analysis configuration, builds dynamic comparison selectors, synchronizes alpha/confidence, applies reporting choices and validates browser defaults. | Sourced inside the server. Cooperates with `worksheet.js`; browser defaults contain choices, not measurements or table mappings. |
| [`R/results.R`](R/results.R) | Defines read-only sheet rows, headers, formats and highlights. Arranges descriptive, comparison, model, normality, outlier and nested results, and supplies arranged/raw CSV export choices. | Numeric calculations remain in the result's raw tables. Plain data-frame sheets omit the otherwise unused row-label column. |
| [`R/results_ui.R`](R/results_ui.R) | Turns result-sheet definitions into Shiny HTML tables, with row numbers, spreadsheet letters and highlights. | Uses `sheet_frame()` from `results.R`; performs no statistical calculations. |
| [`R/plots.R`](R/plots.R) | Constructs data and diagnostic plots, QQ plots, comparison intervals, rank/difference graphs and saved-graph layouts. `make_graph()` selects the requested kind. | Uses result objects from the analysis engine. Nested graphs come from the nested adapter; separate-panel layouts support both ggplots and grid grobs. |
| [`www/worksheet.js`](www/worksheet.js) | Browser event bridge for edits, navigator clicks, keyboard movement and storage/restoration of ordinary analysis defaults. Sends events to Shiny with worksheet revisions and event nonces. | Server code validates changes. Browser storage does not contain worksheet measurements. |
| [`www/style.css`](www/style.css) | Controls page and navigator layout, readable selected/hover states, editable worksheet cells, result sheets and responsive sizes. | Applied by `ui.R`. Later workspace rules refine earlier base rules, so cascade order matters. |

## Nested model files

| File | What it does | Important relationships |
| --- | --- | --- |
| [`R/nested.R`](R/nested.R) | Defines nested defaults/examples; validates Mouse/Slice/Cell mappings and comparisons; runs a bundled reference with temporary paths; collects tables, plots, files and warnings. Builds HTML/Markdown/ZIP exports and preserves named settings through JSON. | Called by `analysis.R`. Reference runs use a separate variable environment and restore emmeans options afterward. Reports embed images as single-line data URIs. |
| [`R/nested_ui.R`](R/nested_ui.R) | Builds nested design, fitting, comparison, outlier, figure and export controls. Parses named labels/colours and turns widget values into nested configuration. | Uses `nested_defaults()`; selecting a control does not change the experimental column sequence. |
| [`R/nested_server.inc`](R/nested_server.inc) | Manages nested mappings and pair/control selectors, pending restored selections, the publication PNG preview, the embedded report and nested downloads. | Sourced inside the server. Pending values allow restoration even when a tab's dynamic outputs are hidden. |
| [`R/nested-reference-repeated.inc`](R/nested-reference-repeated.inc) | Bundled repeated-cell model implementation for wide data: Mouse/Slice/Cell plus condition columns. Handles slice-by-condition variability, residual noise, inference, ROUT, diagnostics, plots and report text. | Executed by `nested.R` with injected settings and temporary paths. Its existing section comments explain the detailed calculations. |
| [`R/nested-reference-between-groups.inc`](R/nested-reference-between-groups.inc) | Bundled between-group implementation for long data: Mouse/Slice/Cell/Condition/Value. Models mouse and within-mouse slice effects with shared or separate condition variances, then produces comparisons, diagnostics, figures and reports. | Executed by `nested.R`. A slice must not appear under several conditions in this design. |

## Test files

Run the complete suite from the project root after package setup:

```r
source("tests/run_tests.R")
```

Or from a terminal where `Rscript` is available:

```text
Rscript tests/run_tests.R
```

| File | What it verifies | How it is run |
| --- | --- | --- |
| [`tests/run_tests.R`](tests/run_tests.R) | Main regression entry point. Defines shared assertions, compares statistical results to references, checks invalid inputs, builds figures, runs exported scripts and sources focused suites. | Run directly from the project root. Outputs go into `artifacts/`. |
| [`tests/t_tests.R`](tests/t_tests.R) | Checks normal/lognormal t-tests and tails against R references, multiple-test families, directional intervals, saved projects and independent exports. | Included by `run_tests.R` or run directly. |
| [`tests/t_test_server_tests.R`](tests/t_test_server_tests.R) | Checks t-test controls, hypothesis freshness, mappings, defaults and archived results. | Included by `server_tests.R` or run directly. |
| [`tests/order_tests.R`](tests/order_tests.R) | Worksheet order survives reversed selector order and missing first measurements in summaries, comparisons, axes and legends. | Included by `run_tests.R`; uses its helpers. |
| [`tests/project_tests.R`](tests/project_tests.R) | Numerical normality/GG/matching/ROUT references, transformations, typed JSON, saved snapshots and figure-layout persistence. | Included by `run_tests.R`; uses shared fixtures/helpers. |
| [`tests/options_tests.R`](tests/options_tests.R) | Automatic method routing, comparison choices, FDR/rank procedures, covariance policies, diagnostics and project/standalone reproducibility. | Included by `run_tests.R`; uses its helpers. |
| [`tests/results_tests.R`](tests/results_tests.R) | Entered versus fitted counts, arranged result sheets, original subject IDs in cleaned exports, reporting switches and numeric export precision. | Included by `run_tests.R`; uses its helpers. |
| [`tests/server_tests.R`](tests/server_tests.R) | Shiny editing, stale results, repeated Analyze clicks, archived navigation/highlighting, transformations, graphs and project workflows. Uses `shiny::testServer()` without a browser. | Included by `run_tests.R`; also loads the next three server suites. |
| [`tests/upload_tests.R`](tests/upload_tests.R) | Per-file upload/edit states, bad-file isolation, CSV decimal choices, Excel, project restoration and readable missing-input errors. | Included by `server_tests.R`; can also initialize its own test app. |
| [`tests/option_server_tests.R`](tests/option_server_tests.R) | Advanced controls become the intended configuration, graphs and reporting choices; restored defaults are validated. | Included by `server_tests.R`; requires its `app_env`. |
| [`tests/refresh_tests.R`](tests/refresh_tests.R) | Correction menus update correctly as methods/comparison approaches change; valid selections and pending defaults survive refreshes. | Included by `server_tests.R`; can initialize `app_env` when run independently. |
| [`tests/nested_tests.R`](tests/nested_tests.R) | Direct-fit nested estimates/contrasts, IDs and invalid pairs, missing measurements, ROUT, report images, publication exports, JSON settings and standalone reproduction. | Included by `run_tests.R`, or run directly from the project root. Outputs go into `artifacts/nested-tests/`. |

## Example data and supporting files

CSV examples are data, so they contain headers and measurements rather than code comments. The nested templates are generated by `nested_demo()` and therefore have no separate CSV here.

| File | Purpose |
| --- | --- |
| [`examples/independent_groups.csv`](examples/independent_groups.csv) | Independent-group measurements for group summaries and comparisons. |
| [`examples/paired.csv`](examples/paired.csv) | Matched before/after measurements with subject identifiers. |
| [`examples/xy.csv`](examples/xy.csv) | X/Y measurements for regression and correlation. |
| [`examples/dose.csv`](examples/dose.csv) | Concentration/response measurements for the four-parameter dose-response model. |
| [`examples/contingency.csv`](examples/contingency.csv) | Category counts for chi-square or Fisher tests. |
| [`examples/factorial.csv`](examples/factorial.csv) | Two-factor measurements for the balanced independent two-way ANOVA example. |
| [`examples/survival.csv`](examples/survival.csv) | Event times, censoring indicators and groups for survival analysis. |
| [`LabStats.Rproj`](LabStats.Rproj) | RStudio project settings: UTF-8, indentation, code indexing and no automatic workspace restoration/saving. Open it to work from the correct project directory. |
| [`README.md`](README.md) | User-facing installation, workflow, statistical conventions, nested-model instructions and current limitations. |
| [`FILE_GUIDE.md`](FILE_GUIDE.md) | This file: source responsibilities, connections and test entry points. |
| [`.gitignore`](.gitignore) | Excludes local dependencies, R session state, generated test artifacts, caches and logs from Git. |
| [`LICENSE`](LICENSE) | Project licensing terms. |

## Generated directories

| Directory | What it contains |
| --- | --- |
| `.R-library/` | Installed package code and dependencies. Managed by setup; not maintained application source. |
| `.cache/` | Runtime/Sass caches and local development checks. It can include generated scripts and screenshots. |
| `artifacts/` | Test-generated tables, plots, projects, reports and reproduction scripts. Tests may create or replace these files. |
| `.Rproj.user/` | Local RStudio session/editor state. |
| `app_cache/` | Local cache directory excluded from Git; not required maintained source. |

For a statistical change, start with `analysis.R`, `options.R` or the nested adapter/references and the matching tests. For a table-layout change, start with `results.R`, `results_ui.R` and `style.css`. For interaction/state bugs, start with `app.R`, the server includes and `worksheet.js`.
