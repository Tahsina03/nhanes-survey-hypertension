# High blood pressure in US adults, NHANES 2009-2012:
# why survey weights and the complex sampling design matter.
#
# Data: NHANESraw from the R package `NHANES` (public NHANES 2009-10 and 2011-12
#       data, with the original MEC exam weights, strata and PSUs)
# Methods: design-based estimation with the `survey` package
#   - combined 4-year MEC weights (WTMEC2YR / 2)
#   - Taylor-series variance with strata (SDMVSTRA) and PSUs (SDMVPSU)
#   - correct subpopulation analysis (subset the design, not the data)
#   - weighted vs unweighted prevalence, design effects
#   - survey-weighted logistic regression (svyglm) vs ordinary glm
#
# Run from the repository root:  Rscript R/analysis.R

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
options(survey.lonely.psu = "adjust")

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
load_nhanes <- function() {
  if (requireNamespace("NHANES", quietly = TRUE)) return(NHANES::NHANESraw)
  if (file.exists("data/NHANESraw.rda")) {              # local copy, not committed
    e <- new.env(); load("data/NHANESraw.rda", envir = e); return(e$NHANESraw)
  }
  stop("Install the NHANES package first: install.packages('NHANES')")
}
d <- as.data.frame(load_nhanes())

# Combined 2-cycle weight (NCHS guidance: divide 2-year weight by number of cycles)
d$WTMEC4YR <- d$WTMEC2YR / 2

# Outcome: measured high blood pressure, 2017 ACC/AHA threshold (>= 130/80)
d$high_bp <- as.integer(d$BPSysAve >= 130 | d$BPDiaAve >= 80)
# Sensitivity outcome: older JNC 7 threshold (>= 140/90)
d$high_bp_140 <- as.integer(d$BPSysAve >= 140 | d$BPDiaAve >= 90)

d$age_grp <- cut(d$Age, c(20, 40, 60, Inf), right = FALSE,
                 labels = c("20-39", "40-59", "60+"))
d$sex     <- factor(d$Gender, levels = c("female", "male"), labels = c("Female", "Male"))
d$race    <- factor(d$Race1, levels = c("White", "Black", "Mexican", "Hispanic", "Other"),
                    labels = c("Non-Hispanic White", "Non-Hispanic Black", "Mexican American",
                               "Other Hispanic", "Other / multiracial"))
d$bmi_cat <- factor(d$BMI_WHO, levels = c("18.5_to_24.9", "12.0_18.5", "25.0_to_29.9", "30.0_plus"),
                    labels = c("Normal", "Underweight", "Overweight", "Obese"))
d$pir_cat <- cut(d$Poverty, c(0, 1, 2, 4, Inf), right = FALSE,
                 labels = c("< 1", "1 to < 2", "2 to < 4", ">= 4"))
d$pir_desc <- d$pir_cat                       # natural order, for descriptive tables
d$pir_cat <- relevel(d$pir_cat, ref = ">= 4")  # highest income as regression reference
d$ever_smoke  <- factor(d$Smoke100, levels = c("No", "Yes"))
d$phys_active <- factor(d$PhysActive, levels = c("Yes", "No"))
d$diabetes    <- factor(d$Diabetes, levels = c("No", "Yes"))

# Analysis population: adults 20+ examined in the MEC with a BP reading
d$in_pop <- with(d, !is.na(Age) & Age >= 20 & WTMEC2YR > 0 & !is.na(high_bp))

model_vars <- c("high_bp", "age_grp", "sex", "race", "bmi_cat", "pir_cat",
                "ever_smoke", "phys_active", "diabetes")
d$in_model <- d$in_pop & complete.cases(d[, model_vars])

cat("Total records:", nrow(d), "\n")
cat("Adults 20+ with BP measured:", sum(d$in_pop), "\n")
cat("Complete cases for regression:", sum(d$in_model), "\n")

# ---------------------------------------------------------------------------
# 2. Survey design on the FULL data, then subset (correct subpopulation analysis)
# ---------------------------------------------------------------------------
d_mec <- d[d$WTMEC2YR > 0, ]   # zero-weight records (interview only) carry no info
des_full <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR,
                      nest = TRUE, data = d_mec)
des <- subset(des_full, in_pop)

# Wrong design for comparison: weights only, ignoring strata and clustering
des_wt_only <- subset(svydesign(ids = ~1, weights = ~WTMEC4YR, data = d_mec), in_pop)

adults <- d[d$in_pop, ]

# ---------------------------------------------------------------------------
# 2b. Why weighting matters: sample composition vs weighted population composition
# ---------------------------------------------------------------------------
comp <- do.call(rbind, lapply(list(c("age_grp", "Age group"), c("race", "Race/ethnicity"),
                                   c("pir_desc", "Income-to-poverty ratio")), function(v) {
  w <- svymean(as.formula(paste0("~", v[1])), des, na.rm = TRUE)
  s <- prop.table(table(adults[[v[1]]]))
  data.frame(variable = v[2], level = names(s),
             sample_pct = round(100 * as.numeric(s), 1),
             population_pct = round(100 * as.numeric(coef(w)), 1), row.names = NULL)
}))
write.csv(comp, "output/tables/sample_vs_population.csv", row.names = FALSE)
print(comp)

# ---------------------------------------------------------------------------
# 3. Overall prevalence: unweighted vs weighted, and design effects
# ---------------------------------------------------------------------------
p_unw <- mean(adults$high_bp)
se_unw <- sqrt(p_unw * (1 - p_unw) / nrow(adults))
w_design  <- svymean(~high_bp, des, deff = TRUE)
w_wtonly  <- svymean(~high_bp, des_wt_only)

overall <- data.frame(
  approach = c("Unweighted (simple random sample assumption)",
               "Weighted, ignoring strata and PSUs",
               "Weighted, full complex design"),
  prevalence_pct = round(100 * c(p_unw, coef(w_wtonly), coef(w_design)), 1),
  se_pct = round(100 * c(se_unw, SE(w_wtonly), SE(w_design)), 2)
)
overall$ci95 <- sprintf("%.1f to %.1f", overall$prevalence_pct - 1.96 * overall$se_pct,
                        overall$prevalence_pct + 1.96 * overall$se_pct)
write.csv(overall, "output/tables/overall_prevalence.csv", row.names = FALSE)
print(overall)
deff_overall <- round(deff(w_design), 2)
cat("Design effect (overall):", deff_overall, "\n")

w140 <- svymean(~high_bp_140, des)
cat("Weighted prevalence at >=140/90:", round(100 * coef(w140), 1), "%\n")

# ---------------------------------------------------------------------------
# 4. Prevalence by subgroup: unweighted vs design-based
# ---------------------------------------------------------------------------
by_group <- function(var, label) {
  f <- as.formula(paste0("~", var))
  w <- svyby(~high_bp, f, des, svymean, deff = TRUE, na.rm = TRUE)
  u <- aggregate(adults$high_bp, list(level = adults[[var]]), mean)
  data.frame(
    variable = label,
    level = as.character(w[[var]]),
    unweighted_pct = round(100 * u$x[match(w[[var]], u$level)], 1),
    weighted_pct = round(100 * w$high_bp, 1),
    lower95 = round(100 * (w$high_bp - 1.96 * w$se), 1),
    upper95 = round(100 * (w$high_bp + 1.96 * w$se), 1),
    deff = round(w[, grep("^DEff", names(w))], 2),
    row.names = NULL
  )
}
sub_tab <- rbind(
  by_group("age_grp", "Age group"),
  by_group("sex", "Sex"),
  by_group("race", "Race/ethnicity"),
  by_group("pir_desc", "Income-to-poverty ratio")
)
write.csv(sub_tab, "output/tables/prevalence_by_group.csv", row.names = FALSE)
print(sub_tab)

plot_df <- rbind(
  data.frame(sub_tab[, c("variable", "level")], pct = sub_tab$unweighted_pct,
             lo = NA, hi = NA, estimate = "Unweighted"),
  data.frame(sub_tab[, c("variable", "level")], pct = sub_tab$weighted_pct,
             lo = sub_tab$lower95, hi = sub_tab$upper95, estimate = "Survey-weighted (95% CI)")
)
plot_df$level <- factor(plot_df$level, levels = rev(unique(sub_tab$level)))
plot_df$variable <- factor(plot_df$variable, levels = unique(sub_tab$variable))
p1 <- ggplot(plot_df, aes(pct, level, colour = estimate, shape = estimate)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, na.rm = TRUE) +
  geom_point(size = 2.4) +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("Unweighted" = "grey55", "Survey-weighted (95% CI)" = "#1b6ca8")) +
  scale_shape_manual(values = c("Unweighted" = 1, "Survey-weighted (95% CI)" = 16)) +
  labs(x = "Prevalence of high blood pressure (%)", y = NULL, colour = NULL, shape = NULL,
       title = "High blood pressure (>= 130/80 mmHg) in US adults, NHANES 2009-2012",
       subtitle = "Unweighted sample estimates vs design-based population estimates") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0))
ggsave("output/figures/prevalence_weighted_vs_unweighted.png", p1,
       width = 8.5, height = 6.5, dpi = 150, bg = "white")

# ---------------------------------------------------------------------------
# 5. Logistic regression: svyglm vs ordinary glm
# ---------------------------------------------------------------------------
fml <- high_bp ~ age_grp + sex + race + bmi_cat + pir_cat + ever_smoke + phys_active + diabetes
des_model <- subset(des_full, in_model)
fit_svy <- svyglm(fml, design = des_model, family = quasibinomial())
fit_glm <- glm(fml, data = d[d$in_model, ], family = binomial())

or_table <- function(fit, label) {
  est <- coef(fit); se <- sqrt(diag(vcov(fit)))
  data.frame(term = names(est), model = label,
             OR = exp(est), lower = exp(est - 1.96 * se), upper = exp(est + 1.96 * se),
             se_logOR = se, row.names = NULL)[-1, ]
}
or_tab <- rbind(or_table(fit_svy, "Survey-weighted (svyglm)"),
                or_table(fit_glm, "Unweighted (glm)"))
pretty <- c(
  "age_grp40-59" = "Age 40-59 (vs 20-39)", "age_grp60+" = "Age 60+ (vs 20-39)",
  "sexMale" = "Male (vs female)",
  "raceNon-Hispanic Black" = "Non-Hispanic Black (vs White)",
  "raceMexican American" = "Mexican American (vs White)",
  "raceOther Hispanic" = "Other Hispanic (vs White)",
  "raceOther / multiracial" = "Other / multiracial (vs White)",
  "bmi_catUnderweight" = "Underweight (vs normal BMI)",
  "bmi_catOverweight" = "Overweight (vs normal BMI)",
  "bmi_catObese" = "Obese (vs normal BMI)",
  "pir_cat< 1" = "Income/poverty < 1 (vs >= 4)",
  "pir_cat1 to < 2" = "Income/poverty 1 to < 2 (vs >= 4)",
  "pir_cat2 to < 4" = "Income/poverty 2 to < 4 (vs >= 4)",
  "ever_smokeYes" = "Ever smoked 100 cigarettes",
  "phys_activeNo" = "Not physically active",
  "diabetesYes" = "Diabetes")
or_tab$term_label <- factor(pretty[or_tab$term], levels = rev(pretty))
or_out <- or_tab
or_out[, c("OR", "lower", "upper")] <- round(or_out[, c("OR", "lower", "upper")], 2)
or_out$se_logOR <- round(or_out$se_logOR, 3)
write.csv(or_out[, c("term_label", "model", "OR", "lower", "upper", "se_logOR")],
          "output/tables/odds_ratios.csv", row.names = FALSE)

se_ratio <- data.frame(
  term = pretty[names(coef(fit_svy))[-1]],
  se_svy = sqrt(diag(vcov(fit_svy)))[-1],
  se_glm = sqrt(diag(vcov(fit_glm)))[-1],
  row.names = NULL
)
se_ratio$ratio <- round(se_ratio$se_svy / se_ratio$se_glm, 2)
se_ratio[, 2:3] <- round(se_ratio[, 2:3], 3)
write.csv(se_ratio, "output/tables/se_ratio_svy_vs_glm.csv", row.names = FALSE)
print(or_out[, c("term_label", "model", "OR", "lower", "upper")])
cat("Median SE ratio (svyglm / glm):", median(se_ratio$ratio), "\n")

p2 <- ggplot(or_tab, aes(OR, term_label, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.6),
                  size = 0.3) +
  scale_x_log10() +
  scale_colour_manual(values = c("Survey-weighted (svyglm)" = "#1b6ca8",
                                 "Unweighted (glm)" = "grey55")) +
  labs(x = "Odds ratio for high blood pressure (log scale, 95% CI)", y = NULL, colour = NULL,
       title = "Survey-weighted vs unweighted logistic regression") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("output/figures/odds_ratios_svyglm_vs_glm.png", p2,
       width = 8.5, height = 6.5, dpi = 150, bg = "white")

# Design-based Wald tests for each block of terms
wald <- do.call(rbind, lapply(c("age_grp", "race", "bmi_cat", "pir_cat"), function(v) {
  t <- regTermTest(fit_svy, as.formula(paste0("~", v)))
  data.frame(term = v, F = round(as.numeric(t$Ftest), 2), df = t$df,
             ddf = t$ddf, p = signif(as.numeric(t$p), 3))
}))
write.csv(wald, "output/tables/wald_tests.csv", row.names = FALSE)
print(wald)

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat("\nDone. Results written to output/\n")
