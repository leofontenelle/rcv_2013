between <- function(x, left, right) {
  x >= left & x <= right
}

diagnose_diabetes <- function(Q029, Q030, Z034) {
  stopifnot(length(Q029) == length(Q030), length(Q030) == length(Z034))
  
  # Q030: Did any doctor diagnose you as having diabetes?
  #       1) yes; 2) only during pregnancy; 3) No
  reported <- Q030 == 1
  # Q029: when was the last time you had your blood glucose measured?
  #       1 through 5) various time frames; 6) never
  reported[is.na(reported) & Q029 == 6] <- FALSE
  
  # Z034: glycated hemoglobin
  measured <- Z034 >= 6.5
  
  # When either the reported or the measured values are not missing, 
  # use what we have.
  res <- rep(NA, times = length(Q029))
  res <- measured %in% TRUE | reported %in% TRUE
  res[is.na(measured) & is.na(reported)] <- NA
  
  return(as.integer(res))
}

impute_bp_meds <- function(Q001, Q002, Q006) {
  stopifnot(length(Q001) == length(Q002), length(Q002) == length(Q006))
  
  # Q001: When was your blood pressure last measured?
  #       6 is never; the questionnaire skips the other questions.
  # Q002: Did any doctor diagnose you as having hypertension?
  #       1) yes, 2) only during pregnancy, 3) no
  #       if 2 or 3 skip the other question.
  # Q006: During the last two weeks, did you take any medicines 
  #       because of hypertension?
  #       1) yes, 2) no
  res <- 2L - as.integer(Q006)
  res[is.na(res) & Q002 > 1] <- 0
  res[is.na(res) & Q001 == 6] <- 0
  
  return(res)
}
