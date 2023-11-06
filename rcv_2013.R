# Agreement between different scores for the estimation of 
# cardiovascular risk in Brazil, 2013
#

# Load packages and functions ----

# https://github.com/boyercb/globorisk
load("contrib/sysdata.rda") # coefs, cvdr, rf
source("contrib/globorisk.R")

library(PooledCohort)
library(CVrisk)

source("functions.R")

escores <- c(framingham = "Framingham",
             pooled_cohort = "Pooled Cohort Equations",
             globorisk = "Globorisk-LAC")


# Read data ----

# https://www.pns.icict.fiocruz.br/bases-de-dados/
pns2013_lab <- read.csv2("data/pns2013_exames_2023-05-05.csv", 
                         na.strings = "", fileEncoding = "UTF-8")

d <- with(pns2013_lab, data.frame(
  sex = Z001 - 1,
  age_years = Z002,
  race = ifelse(Z003 <= 5, Z003, NA), # 9 is "ignored"
  smoke_current = as.integer(P050 %in% 1:2),
  chol_total_mgdl = Z031,
  chol_hdl_mgdl = Z032,
  bp_sys_mmhg = W00407,
  bp_meds = impute_bp_meds(Q001, Q002, Q006),
  diabetes = diagnose_diabetes(Q029, Q030, Z034),
  cvd = Q063 == 1 | Q068 == 1, # cardiac disease or stroke
  survey_weight = peso_lab
))

rm(pns2013_lab)


# Estimate CVD risk ----

# Framingham: 30 to 74 years of age
# Pooled Cohort Equations: 40 to 79 years of age
# Globorisk-LAC: 40 to 74 years of age
d$ok_age <- d$age_years |> between(40, 74)
d$ok_complete <- complete.cases(d)
d$ok_nocvd <- !d$cvd
d$ok_notatypical <- # only in Globorisk-LAC
  between(d$bp_sys_mmhg, 70, 270) & 
  between(d$chol_total_mgdl, 67.77, 773.3)
d$ok <- d$ok_age & d$ok_complete & d$ok_nocvd & d$ok_notatypical

# ascvd_10y_frs() returns risk as a percentage, while the other two
# return risk as a proportion.
d$framingham[d$ok] <- with(subset(d, ok), ascvd_10y_frs(
  gender = ifelse(sex == 0, "male", "female"),
  age = age_years,
  hdl = chol_hdl_mgdl,
  totchol = chol_total_mgdl,
  sbp = bp_sys_mmhg,
  bp_med = bp_meds,
  smoker = smoke_current,
  diabetes = diabetes
) / 100)

d$pooled_cohort[d$ok] <- with(subset(d, ok), predict_10yr_ascvd_risk(
  age_years = age_years,
  race = race,
  sex = sex,
  smoke_current = smoke_current,
  chol_total_mgdl = chol_total_mgdl,
  chol_hdl_mgdl = chol_hdl_mgdl,
  bp_sys_mmhg = bp_sys_mmhg,
  bp_meds = bp_meds,
  diabetes = diabetes,
  equation_version = "Yadlowsky_2018",
  override_boundary_errors = TRUE,
  race_levels = list(white = c(1, 3, 5), black = c(2, 4)),
  sex_levels = list(male = 0, female = 1),
  smoke_current_levels = list(no = 0, yes = 1),
  bp_meds_levels =  list(no = 0, yes = 1),
  diabetes_levels = list(no = 0, yes = 1)
))

d$globorisk[d$ok] <- with(subset(d, ok), globorisk(
  sex = sex,
  # globorisk always throws an error when the age is lower than 40
  age = ifelse(age_years >= 40, age_years, NA_integer_),
  sbp = bp_sys_mmhg,
  tc = chol_total_mgdl * 10 / 386.65354,
  dm = diabetes,
  smk = smoke_current,
  iso = "BRA",
  year = 2013,
  version = "lab",
  type = "risk",
  updated_lac = TRUE
))

stopifnot(with(subset(d, ok), all(
  !is.na(framingham) & !is.na(pooled_cohort) & !is.na(globorisk)
)))

for (score in c("framingham", "pooled_cohort", "globorisk")) {
  new_name <- sprintf("%s_cat", score)
  d[[new_name]] <- d[[score]] |> 
    cut(breaks = c(0, 0.1, 0.2, Inf), 
        labels = c("Low", "Intermediate", "High"),
        right = FALSE, ordered_results = TRUE) |> 
    ordered() # somehow cut(..., ordered_results = TRUE) is not enough
}


# Describe sample ----

sample_size <- list(
  n_total = nrow(d),
  # prop_total is 1.0 or 100%
  n_age = sum(d$ok_age),
  prop_age = weighted.mean(d$ok_age, d$survey_weight),
  n_complete = d |> 
    with(sum(ok_age & ok_complete)),
  prop_complete = d |> 
    with(weighted.mean(ok_age & ok_complete, survey_weight)),
  n_nocvd = d |> 
    with(sum(ok_age & ok_complete & ok_nocvd)),
  prop_nocvd = d |> 
    with(weighted.mean(ok_age & ok_complete & ok_nocvd, survey_weight)),
  # d$ok <- d$ok_age & d$ok_complete & d$ok_nocvd & d$ok_notatypical
  n_ok = sum(d$ok),
  prop_ok = weighted.mean(d$ok, d$survey_weight)
)

sample_descr <- with(subset(d, ok), list(
  n_female = sum(sex),
  prop_female = weighted.mean(sex, survey_weight),
  mean_age = weighted.mean(age_years, survey_weight),
  sd_age = weighted.sd(age_years, survey_weight),
  n_black = sum(race %in% c(2, 4)),
  prop_black = weighted.mean(race %in% c(2, 4), survey_weight),
  n_smoke = sum(smoke_current),
  prop_smoke = weighted.mean(smoke_current, survey_weight),
  mean_chol_total = weighted.mean(chol_total_mgdl, survey_weight),
  sd_chol_total = weighted.sd(chol_total_mgdl, survey_weight),
  mean_chol_hdl = weighted.mean(chol_hdl_mgdl, survey_weight),
  sd_chol_hdl = weighted.sd(chol_hdl_mgdl, survey_weight),
  mean_bp = weighted.mean(bp_sys_mmhg, survey_weight),
  sd_bp = weighted.sd(bp_sys_mmhg, survey_weight),
  n_bp_meds = sum(bp_meds),
  prop_bp_meds = weighted.mean(bp_meds, survey_weight),
  n_diabetes = sum(diabetes),
  prop_diabetes = weighted.mean(diabetes, survey_weight)
))

risk_descr <- (\(d, w) {
  cbind(
   t(sapply(d, weighted.quantile, w = w, probs = c(0, 0.25, 0.5, 0.75, 1.0))),
   mean = sapply(d, weighted.mean, w = w),
   sd = sapply(d, weighted.sd, w = w)
)}) (d[d$ok, names(escores)], d$survey_weight[d$ok])

risk_cat_descr <- with(subset(d, ok), {
  res <- list("Framingham" = framingham_cat, 
              "Pooled Cohort Equations" = pooled_cohort_cat, 
              "Globorisk-LAC" = globorisk_cat) |> 
    lapply(\(x) data.frame(
      level = levels(x),
      n = tapply(survey_weight, x, length),
      prop = tapply(survey_weight, x, sum) / sum(survey_weight)
    ))
  for (nm in names(res)) res[[nm]] <- cbind(score = nm, res[[nm]])
  res <- do.call(rbind.data.frame, c(res, make.row.names = FALSE))
  res
})

# Compare estimates ----

d$ratio_fg <- d$framingham / d$globorisk
d$ratio_fp <- d$framingham / d$pooled_cohort
d$ratio_pg <- d$pooled_cohort / d$globorisk

agree_cont <- data.frame(
  A = escores[c("framingham", "framingham", "pooled_cohort")],
  B = escores[c("globorisk", "pooled_cohort", "globorisk")],
  A_lower = with(subset(d, ok), 
                 sapply(list(ratio_fg, ratio_fp, ratio_pg), 
                        \(x) weighted.mean(x < (1/1.25), survey_weight))),
  Agreement = with(subset(d, ok), 
                   sapply(list(ratio_fg, ratio_fp, ratio_pg), 
                          \(x) between(x, 1/1.25, 1.25) |> 
                            weighted.mean(survey_weight))),
  A_higher = with(subset(d, ok), 
                  sapply(list(ratio_fg, ratio_fp, ratio_pg), 
                         \(x) weighted.mean(x > 1.25, survey_weight)))
)

agree_cat <-  data.frame(
  A = agree_cont$A,
  B = agree_cont$B,
  pa = with(subset(d, ok), 
            list(framingham_cat == globorisk_cat,
                 framingham_cat == pooled_cohort_cat,
                 pooled_cohort_cat == globorisk_cat) |> 
              sapply(weighted.mean, survey_weight)),
  ac1 = with(subset(d, ok), c(
    weighted.ac1(framingham_cat, globorisk_cat, survey_weight),
    weighted.ac1(framingham_cat, pooled_cohort_cat, survey_weight),
    weighted.ac1(pooled_cohort_cat, globorisk_cat, survey_weight)
  ))
)

# Tables and figures ----

tab1 <- sample_descr |> 
  with(data.frame(
    Variable = c("Female", "Age (mean, SD)", "Black", "Smoker", 
                 "Total cholesterol (mg/dL) (mean, SD)", 
                 "HDL cholesterol (mg/dL) (mean, SD)",
                 "Blood pressure (mean, SD)", "Using blood pressure medication",
                 "Diabetes"),
    n = c(n_female, mean_age, n_black, n_smoke, 
          mean_chol_total, mean_chol_hdl, 
          mean_bp, n_bp_meds, n_diabetes),
    proportion = c(prop_female, sd_age, prop_black, prop_smoke, 
                   sd_chol_total, sd_chol_hdl, 
                   sd_bp, prop_bp_meds, prop_diabetes)))

with(subset(d, ok), {
  plot(NA, xlim = c(0, 1), ylim = c(0, 15), 
       xlab = "Estimated cardiovascular risk", ylab = NA, yaxt = "n")
  abline(v = 0.1, lty = 2, col = "gray")
  abline(v = 0.2, lty = 2, col = "gray")
  lines(density(framingham, bw = 0.05, from = 0, to = 1,
                weights = survey_weight/sum(survey_weight)), col = 1)
  lines(density(pooled_cohort, bw = 0.05, from = 0, to = 1,
                weights = survey_weight/sum(survey_weight)), col = 2)
  lines(density(globorisk, bw = 0.05, from = 0, to = 1,
                weights = survey_weight/sum(survey_weight)), col = 3)
  legend("topright", col = 1:3, lty = 1, legend = escores)
})



risk_cat_descr$prop |> 
  matrix(nrow = 3, dimnames = list(unique(risk_cat_descr$level), 
                                   unique(risk_cat_descr$score))) |> 
  barplot()

with(subset(d, ok), {
  plot(NA, xlim = c(0.2, 5), ylim = c(0, 1), log = "x",
       xlab = "Ratio between estimates", ylab = NA, yaxt = "n")
  abline(v = 0.8, lty = 2, col = "gray")
  abline(v = 1.25, lty = 2, col = "gray")
  rect(0.8, -3, 1.25, 2, col = "lightgray", border = NA)
  lines(density(ratio_fg, bw = 0.25, 
                weights = survey_weight/sum(survey_weight)), col = 1)
  lines(density(ratio_fp, bw = 0.25, 
                weights = survey_weight/sum(survey_weight)), col = 2)
  lines(density(ratio_pg, bw = 0.25, 
                weights = survey_weight/sum(survey_weight)), col = 3)
  legend("topright", col = 1:3, lty = 1,
         legend = c("Framingham \uf7 Globorisk-LAC", 
                    "Framingham \uf7 PCE", 
                    "PCE \uf7 Globorisk-LAC"))
})

plot_categorical_agreement <- function(x, y, weight = 1, xlab = NULL, ylab = NULL) {
  stopifnot(is.ordered(x), is.ordered(y))
  stopifnot(length(x) == length(y))
  if (length(weight) != length(x)) {
    stop("weight must be missing, length 1, or as long as x and y")
  } 
  if (is.null(xlab)) warning("Label for X axis not informed")
  if (is.null(ylab)) warning("Label for Y axis not informed")
  
  xy <- xtabs(weight ~ y + x)
  totals <- colSums(xy)
  barplot(height = xy / rep(totals, each = nrow(xy)), 
          width = totals, 
          xlab = xlab, ylab = ylab)
}
plot_categorical_agreement(d$globorisk_cat, d$framingham_cat, d$survey_weight, 
                           escores["globorisk"], escores["framingham"])
plot_categorical_agreement(d$pooled_cohort_cat, d$framingham_cat, d$survey_weight, 
                           escores["pooled_cohort"], escores["framingham"])
plot_categorical_agreement(d$globorisk_cat, d$pooled_cohort_cat, d$survey_weight, 
                           escores["globorisk"], escores["pooled_cohort"])


# Write output ----

write.csv2(tab1, "tab1.csv", row.names = FALSE, fileEncoding = "UTF-8")
