# Agreement between Globorisk LAC and Pooled Cohort Equations for the 
# estimation of cardiovascular risk in Brazil, 2013
#

# Load packages and functions ----

# https://github.com/boyercb/globorisk
load("sysdata.rda") # coefs, cvdr, rf
source("globorisk.R")

library(PooledCohort)

source("functions.R")

# Read data ----

# https://www.pns.icict.fiocruz.br/bases-de-dados/
pns2013_lab <- read.csv2("pns2013_exames_2023-05-05.csv", 
                         na.strings = "", fileEncoding = "UTF-8")

d <- with(pns2013_lab, data.frame(
  sex = Z001 - 1,
  age_years = Z002,
  race = ifelse(Z003 <= 5, Z003, NA),
  smoke_current = as.integer(P050 %in% 1:2),
  chol_total_mgdl = Z031,
  chol_hdl_mgdl = Z032,
  bp_sys_mmhg = W00407,
  bp_meds = impute_bp_meds(Q001, Q002, Q006),
  diabetes = diagnose_diabetes(Q029, Q030, Z034),
  survey_weight = peso_lab
))

rm(pns2013_lab)
