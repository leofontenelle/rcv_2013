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
  gender = factor(Z001, labels = c("Male", "Female")),
  age_years = Z002,
  race = Z003 |> 
    factor(labels = c("white", "black", "white", "black", "white", NA)),
  smoke_current = factor(P050, labels = c("yes", "yes", "no")),
  chol_total_mgdl = Z031,
  chol_total_mmoll = Z031 * 10 / 386.65354,
  chol_hdl_mgdl = Z032,
  bp_sys_mmHg = W00407,
  bp_meds = factor(Q060, labels = c("yes", "no")),
  diabetes = diagnose_diabetes(Q029, Q030, Z034) |> 
    factor(labels = c("no", "yes")),
  survey_weight = peso_lab
))

rm(pns2013_lab)