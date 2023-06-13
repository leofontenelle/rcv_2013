diagnose_diabetes <- function(Q029, Q030, Z034) {
  stopifnot(length(Q029) == length(Q030), length(Q030) == length(Z034))
  
  # Q030: Did any doctor diagnosed you as having diabetes?
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
  
  return(res)
}
