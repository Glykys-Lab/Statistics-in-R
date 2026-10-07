if (dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"), .libPaths()))
source("R/analysis.R"); source("R/data.R"); source("R/plots.R")
invisible(parse("app.R"))
invisible(parse("R/ui.R"))
checks <- 0L
check <- function(ok, description) {
  if (!isTRUE(ok)) stop("FAILED: ", description)
  checks <<- checks + 1L; cat("PASS:", description, "\n")
}
near <- function(x, y, tolerance = 1e-9) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
fails <- function(expr, pattern) { e <- tryCatch({ force(expr); NULL }, error = identity); !is.null(e) && grepl(pattern, conditionMessage(e), ignore.case = TRUE) }
dir.create("artifacts", showWarnings = FALSE)
excel <- readxl::read_excel(readxl::readxl_example("datasets.xlsx"), sheet = 1)
check(nrow(excel) > 0 && ncol(excel) > 0, "Excel reader works with an actual workbook")
d <- demo_data("groups"); pair <- demo_data("paired")
cfg <- list(y = "response", group = "group")
r <- run_analysis(d, c(cfg, method = "anova"))
reference <- aov(response ~ group, d)
check(near(r$tables$anova[["Pr(>F)"]][1], summary(reference)[[1]][["Pr(>F)"]][1]), "ANOVA p-value matches R reference")
check(near(r$tables$posthoc[["p adj"]], TukeyHSD(reference)$group[, "p adj"]), "Tukey adjusted p-values match R reference")
two <- d[d$group != "Drug B", ]; v <- split(two$response, two$group)
for (method in c("welch", "student", "mann_whitney")) {
  r <- run_analysis(two, c(cfg, method = method))
  ref <- if (method == "mann_whitney") wilcox.test(v[[1]], v[[2]], exact = FALSE) else t.test(v[[1]], v[[2]], var.equal = method == "student")
  check(near(r$tables$test$p_value, ref$p.value), paste(method, "matches R reference"))
}
set.seed(912)
shuffled <- pair[sample(nrow(pair)), ]
r <- run_analysis(shuffled, c(cfg, method = "paired_t", id = "subject"))
vp <- split(pair$response, pair$group)
check(near(r$tables$test$p_value, t.test(vp[[1]], vp[[2]], paired = TRUE)$p.value), "Paired test matches by ID after row shuffling")
check(near(r$tables$test$estimate, mean(vp[[1]] - vp[[2]])), "Paired difference sign follows named contrast")
partial <- pair; partial$response[1] <- NA
r <- run_analysis(partial, c(cfg, method = "paired_t", id = "subject"))
check(r$tables$data_audit$incomplete_subjects == 1 && r$tables$data_audit$used_rows == 14, "Incomplete subject is excluded from both conditions")
missing_condition <- d; missing_condition$response[missing_condition$group == "Drug B"] <- NA
check(fails(run_analysis(missing_condition, c(cfg, method = "friedman", id = "subject")), "complete matched"), "A fully missing repeated condition cannot silently disappear")
check(fails(run_analysis(rbind(pair, pair[1, ]), c(cfg, method = "paired_t", id = "subject")), "Duplicate"), "Duplicate paired records rejected")
check(fails(run_analysis(d, c(cfg, method = "welch")), "exactly two"), "Two-group test rejects three groups")
bad <- two; bad$response[1] <- Inf
check(fails(run_analysis(bad, c(cfg, method = "welch")), "infinite"), "Infinite observations rejected")
bad <- two; bad$response <- as.character(bad$response)
check(fails(run_analysis(bad, c(cfg, method = "welch")), "numeric"), "Text responses rejected without silent conversion")
for (method in c("summary", "welch_anova", "kruskal", "friedman")) {
  r <- run_analysis(d, c(cfg, method = method, id = "subject"))
  check(r$tables$data_audit$used_rows == 24, paste(method, "runs on example data"))
}
r <- run_analysis(pair, c(cfg, method = "paired_wilcox", id = "subject"))
check(near(r$tables$test$p_value, wilcox.test(vp[[1]], vp[[2]], paired = TRUE, exact = FALSE)$p.value), "Paired Wilcoxon matches reference")
xy <- demo_data("xy")
for (method in c("pearson", "spearman", "linear")) {
  r <- run_analysis(xy, list(method = method, x = "concentration", y = "response"))
  if (method != "linear") check(near(r$tables$test$p_value, cor.test(xy$concentration, xy$response, method = method, exact = FALSE)$p.value), paste(method, "matches R reference"))
  else check(near(coef(r$fit), coef(lm(response ~ concentration, xy))), "Linear regression matches R reference")
}
dose <- demo_data("dose")
r <- run_analysis(dose, list(method = "dose", x = "concentration", y = "response"))
check(abs(r$tables$ec50_interval$EC50 - 1.8) < .01, "4PL recovers known EC50")
increasing <- dose; increasing$response <- 105 - dose$response
ri <- run_analysis(increasing, list(method = "dose", x = "concentration", y = "response", direction = "increasing"))
check(abs(ri$tables$ec50_interval$EC50 - 1.8) < .01, "Increasing 4PL recovers known EC50")
dose$concentration[1] <- 0
check(fails(run_analysis(dose, list(method = "dose", x = "concentration", y = "response")), "positive"), "4PL rejects zero concentrations")
cont <- demo_data("contingency")
for (method in c("fisher", "chisq")) {
  r <- run_analysis(cont, list(method = method, group = "group", factor2 = "outcome", count = "count"))
  tab <- xtabs(count ~ group + outcome, cont)
  ref <- if (method == "fisher") fisher.test(tab) else chisq.test(tab, correct = FALSE)
  check(near(r$tables$test$p_value, ref$p.value), paste(method, "matches R reference"))
}
cont$count[1] <- -1
check(fails(run_analysis(cont, list(method = "fisher", group = "group", factor2 = "outcome", count = "count")), "nonnegative"), "Negative counts rejected")
f <- demo_data("factorial")
r <- run_analysis(f, c(cfg, method = "two_way", factor2 = "sex"))
check(nrow(r$tables$anova) == 4, "Two-way ANOVA includes interaction")
check(fails(run_analysis(f[-1, ], c(cfg, method = "two_way", factor2 = "sex")), "equal cell"), "Unbalanced two-way design rejected explicitly")
s <- demo_data("survival")
r <- run_analysis(s, list(method = "survival", y = "time", event = "event", group = "group"))
ref <- survival::survdiff(survival::Surv(time, event) ~ group, s)
check(near(r$tables$test$p_value, pchisq(ref$chisq, 1, lower.tail = FALSE)), "Log-rank matches R reference")
s$event[1] <- 2
check(fails(run_analysis(s, list(method = "survival", y = "time", event = "event", group = "group")), "0 .*1"), "Invalid survival event coding rejected")
wide <- data.frame(id = c("a", "b"), Control = c(1, 2), Treatment = c(3, 4))
long <- wide_to_long(wide, c("Control", "Treatment"), "id")
check(identical(long$subject, c("a", "b", "a", "b")) && identical(long$response, c(1, 2, 3, 4)), "Wide data preserves subject alignment")

# Design routing must cover every requested independent/matched variant.
check(unname(design_methods(family = "two")) == "welch", "Independent Gaussian two-group design selects Welch")
check(unname(design_methods(family = "two", variance = "equal")) == "student", "Equal-variance two-group design selects Student")
check(unname(design_methods(family = "two", matched = TRUE)) == "paired_t", "Matched Gaussian two-group design selects paired t")
check(unname(design_methods(family = "two", distribution = "rank")) == "mann_whitney", "Independent rank design selects Mann-Whitney")
check(unname(design_methods(family = "two", matched = TRUE, distribution = "rank")) == "paired_wilcox", "Matched rank design selects Wilcoxon")
check(unname(design_methods(variance = "equal")) == "anova" && unname(design_methods()) == "welch_anova", "Multi-group Gaussian designs distinguish ordinary and Welch ANOVA")
check(unname(design_methods(distribution = "rank")) == "kruskal" && unname(design_methods(matched = TRUE, distribution = "rank")) == "friedman", "Rank multi-group designs distinguish Kruskal-Wallis and Friedman")
check(unname(design_methods(matched = TRUE)) == "rm_anova" && unname(design_methods(matched = TRUE, repeated = "mixed")) == "mixed", "Repeated Gaussian design offers ANOVA and mixed-effects")
check(setequal(unname(design_methods("xy")), c("linear", "pearson", "spearman")), "XY table offers regression and correlation")

rm_cfg <- c(cfg, method = "rm_anova", id = "subject")
rm_result <- run_analysis(d, rm_cfg)
rm_reference <- summary(aov(response ~ group + Error(subject/group), data = transform(d, subject = factor(subject), group = factor(group))))[[2]][[1]]
check(near(rm_result$tables$test$statistic[1], rm_reference[["F value"]][1]), "Repeated-measures F statistic matches R Error(subject/group) reference")
check(near(rm_result$tables$test$p_value[1], rm_reference[["Pr(>F)"]][1]), "Uncorrected repeated-measures p-value matches reference")
check(rm_result$tables$test$epsilon[2] >= .5 && rm_result$tables$test$epsilon[2] <= 1 && rm_result$tables$test$p_value[2] >= rm_result$tables$test$p_value[1], "GG epsilon is bounded and correction does not lower the p-value")
# Independent covariance expression for GG epsilon, on a nontrivial four-condition design.
set.seed(487)
n <- 20; k <- 4
scores <- outer(rnorm(n, sd = 2), rep(1, k)) + matrix(rnorm(n * k), n) %*% diag(c(1, 2, 3, 5)) + rep(c(0, .5, 1, 2), each = n)
four <- data.frame(subject = rep(paste0("P", seq_len(n)), k), group = rep(LETTERS[1:k], each = n), response = c(scores))
gg <- run_analysis(four, rm_cfg)
centering <- diag(k) - matrix(1/k, k, k); centered_cov <- centering %*% cov(scores) %*% centering
expected_epsilon <- sum(diag(centered_cov))^2 / ((k - 1) * sum(centered_cov^2))
check(near(gg$tables$test$epsilon[2], expected_epsilon), "GG epsilon matches centered-covariance formula on four conditions")
check(near(gg$tables$test$p_value[2], pf(gg$tables$test$statistic[2], (k - 1) * expected_epsilon, (n - 1) * (k - 1) * expected_epsilon, lower.tail = FALSE)), "GG p-value uses corrected numerator and denominator degrees of freedom")
set.seed(123)
shuffled_rm <- run_analysis(four[sample(nrow(four)), ], rm_cfg)
check(near(shuffled_rm$tables$test$p_value, gg$tables$test$p_value), "Repeated-measures results preserve ID matching after row shuffle")
rm_pair <- run_analysis(pair, rm_cfg)
check(near(rm_pair$tables$test$statistic[1], t.test(vp[[1]], vp[[2]], paired = TRUE)$statistic^2) && near(rm_pair$tables$test$epsilon[2], 1), "Two-condition repeated ANOVA equals paired t squared with epsilon one")
check(near(rm_result$tables$posthoc$adjusted_p_value, pairwise_table(pairwise.t.test(d$response, d$group, paired = TRUE, pool.sd = FALSE, p.adjust.method = "holm"))$adjusted_p_value), "Repeated ANOVA paired post-hoc p-values match reference")
rm_missing <- d; rm_missing$response[c(2, 13)] <- NA
complete_rm <- run_analysis(rm_missing, rm_cfg)
check(complete_rm$tables$data_audit$incomplete_subjects == 2 && complete_rm$tables$data_audit$used_rows == 18, "Repeated ANOVA excludes entire incomplete subjects")
mixed_cfg <- c(cfg, method = "mixed", id = "subject")
mixed_result <- run_analysis(rm_missing, mixed_cfg)
mixed_reference <- nlme::lme(response ~ group, random = ~1 | subject, data = transform(rm_missing, group = factor(group), subject = factor(subject)), method = "REML", na.action = na.omit)
check(near(mixed_result$tables$test$uncorrected_p_value, anova(mixed_reference, type = "marginal")["group", "p-value"]), "Uncorrected mixed-effects marginal F test matches nlme REML reference")
check(near(mixed_result$tables$coefficients$estimate, nlme::fixef(mixed_reference)), "Mixed-effects coefficients match reference")
check(mixed_result$tables$data_audit$used_rows == 22 && mixed_result$tables$data_audit$retained_incomplete_subjects == 2, "Mixed-effects model retains observations from incomplete subjects")
check(fails(run_analysis(rbind(rm_missing, rm_missing[2, ]), mixed_cfg), "Duplicate"), "Mixed model rejects duplicate records even when missing")
check(fails(run_analysis(missing_condition, mixed_cfg), "Every condition"), "Mixed model rejects a fully unobserved condition")
# One observation per subject-condition; no subject bridges the A/B pair and C.
disconnected <- data.frame(subject = c("s1", "s2", "s1", "s2", "s3", "s4"), group = c("A", "A", "B", "B", "C", "C"), response = c(1, 2, 2, 3, 4, 5))
check(fails(run_analysis(disconnected, mixed_cfg), "disconnected"), "Mixed model rejects disconnected repeated conditions")
normal <- run_analysis(pair, c(cfg, method = "paired_t", id = "subject"))$tables$normality
check(normal$sample == "Within-pair differences" && near(normal$p_value, shapiro.test(vp[[1]] - vp[[2]])$p.value), "Paired normality test uses within-pair differences")
lr <- run_analysis(xy, list(method = "linear", x = "concentration", y = "response"))
check(near(lr$tables$normality$p_value, shapiro.test(residuals(lr$fit))$p.value), "Regression normality test uses residuals")
check(is.na(normality_table(data.frame(y = c(1, 1, 1), group = factor(rep("A", 3))), "summary")$p_value), "Constant samples do not trigger invalid Shapiro tests")
check(is.na(normality_table(data.frame(y = 1:5001, group = factor(rep("A", 5001))), "summary")$p_value), "Normality diagnostic reports large samples without subsampling")
check(near(wide_to_long(column_demo("groups"), c("Control", "Drug A", "Drug B"), "Subject")$response, d$response), "Column example preserves all original measurements")
pasted <- paste_table("Control\tTreatment\n1\t2\n3\t\n5\t6")
check(is.numeric(pasted$Control) && is.na(pasted$Treatment[2]), "Spreadsheet paste keeps blank numeric cells missing")
check(near(paste_table("A;B\n1,5;2,5\n3,5;4,5", ",")$A, c(1.5, 3.5)), "Paste supports semicolon delimiter and decimal commas")
check(fails(paste_table("A,A\n1,2"), "unique"), "Paste rejects duplicate headers")
check(is.na(edit_cell(wide, 1L, 2L, "")$Control[1]), "Clearing a numeric worksheet cell creates NA")
check(fails(edit_cell(wide, 1L, 2L, "abc"), "number"), "Worksheet edits reject text in numeric measurements")
for (res in list(rm_result, mixed_result, lr)) {
  ggplot2::ggplot_build(diagnostic_plot(res)); ggplot2::ggplot_build(normality_plot(res))
}
check(TRUE, "Repeated/mixed residual and normal QQ figures render")

# Render every supported method, including group-style and paired variants.
for (method in unname(analysis_methods)) {
  name <- if (method %in% c("welch", "student", "paired_t", "mann_whitney", "paired_wilcox")) "paired" else if (method %in% c("pearson", "spearman", "linear")) "xy" else if (method %in% c("chisq", "fisher")) "contingency" else if (method == "two_way") "factorial" else if (method %in% c("dose", "survival")) method else "groups"
  config <- list(method = method, y = if (method == "survival") "time" else "response", group = "group", x = "concentration", id = "subject", event = "event", count = "count", factor2 = if (method == "two_way") "sex" else "outcome")
  result <- run_analysis(demo_data(name), config)
  ggplot2::ggsave(file.path("artifacts", paste0(method, ".png")), make_plot(result), width = 7, height = 5, dpi = 100)
  check(file.exists(file.path("artifacts", paste0(method, ".png"))), paste("Figure renders:", method))
}
result <- run_analysis(d, c(cfg, method = "anova"))
for (style in c("box", "violin")) ggplot2::ggplot_build(make_plot(result, style = style))
check(TRUE, "Box and violin figures render")

# Exported script must run independently of app files and reconstruct the result.
writeLines(reproduction_script(d, c(cfg, method = "anova")), "artifacts/reproduce_analysis.R")
script_path <- normalizePath("artifacts/reproduce_analysis.R")
old <- getwd(); setwd("artifacts")
env <- new.env(parent = globalenv()); sys.source(script_path, envir = env)
setwd(old)
check(near(env$result$tables$anova[["Pr(>F)"]][1], result$tables$anova[["Pr(>F)"]][1]), "Standalone exported script reproduces analysis and figures")
for (config in list(rm_cfg, mixed_cfg)) {
  script <- reproduction_script(rm_missing, config)
  path <- file.path("artifacts", paste0("reproduce_", config$method, ".R")); writeLines(script, path)
  script_path <- normalizePath(path); old <- getwd(); setwd("artifacts")
  env <- new.env(parent = globalenv()); sys.source(script_path, envir = env); setwd(old)
  check(near(env$result$tables$test$p_value, run_analysis(rm_missing, config)$tables$test$p_value), paste("Standalone export reproduces", config$method))
}
source("tests/project_tests.R", encoding = "UTF-8")
cat("\n", checks, "checks passed.\n")
source("tests/server_tests.R")
