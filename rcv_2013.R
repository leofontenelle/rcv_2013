# Agreement between Globorisk LAC and Pooled Cohort Equations for the 
# estimation of cardiovascular risk in Brazil, 2013
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
  survey_weight = peso_lab
))

rm(pns2013_lab)


# Estimate CDV risk ----

# This function returns risk as a percentage, while the other two
# return risk as a proportion.
d$framingham <- with(d, ascvd_10y_frs(
  gender = ifelse(sex == 0, "male", "female"),
  age = age_years,
  hdl = chol_hdl_mgdl,
  totchol = chol_total_mgdl,
  sbp = bp_sys_mmhg,
  bp_med = bp_meds,
  smoker = smoke_current,
  diabetes = diabetes
) / 100)

d$pooled_cohort <- with(d, predict_10yr_ascvd_risk(
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

d$globorisk <- with(d, globorisk(
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
