# Agreement between different scores for the estimation of 
# cardiovascular risk in Brazil, 2013
#

# Load packages and functions ----

# https://github.com/boyercb/globorisk
load("sysdata.rda") # coefs, cvdr, rf
source("globorisk.R")

library(PooledCohort)
library(CVrisk)

source("functions.R")


# Read data ----

# https://www.pns.icict.fiocruz.br/bases-de-dados/
pns2013_lab <- read.csv2("pns2013_exames_2023-05-05.csv", 
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
  cvd = Q063 == 1 | Q068 == 1, # cardiac heart disease or stroke
  survey_weight = peso_lab
))

rm(pns2013_lab)


# Estimate CDV risk ----

# Framingham: 30 to 74
# Pooled Cohort Equations: 40 to 79
# Globorisk LAC: 40 to 74
d$ok_age <- d$age_years |> between(40, 74)
d$ok_complete <- complete.cases(d)
d$ok_nocvd <- !d$cvd
d$ok_notatypical <- # only in Globorisk
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
        right = FALSE, ordered_results = TRUE)
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
)}) (d[d$ok, c("framingham", "pooled_cohort", "globorisk")], 
     d$survey_weight[d$ok])
