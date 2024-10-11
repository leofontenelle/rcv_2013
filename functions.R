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


#' Format a vector of quantiles as an interquartile range
#'
#' @param quantiles output of quantile() or weighted.quantile()
#' @param OutDec character to use as the decimal mark
#'
#' @return a character string
#'
iqr <- function(quantiles, OutDec = getOption("OutDec")) {
  if (!("scales" %in% installed.packages())) {
    stop("Function iqr() requires installing package \"scales\"")
  }
  stopifnot(all(c("25%", "75%") %in% names(quantiles)))
  range_char <- scales::number(x = quantiles[c("25%", "75%")], 
                               accuracy = 0.1,
                               decimal.mark = getOption("OutDec"))
  paste0("(", range_char[1], "; ", range_char[2], ")")
}


tabulate_fig3 <- function(from, to, weight = 1) {
  if (!("ggalluvial" %in% installed.packages())) {
    stop("Function tabulate_fig3() requires installing package \"ggalluvial\"")
  }
  
  stopifnot(is.ordered(from), is.ordered(to))
  stopifnot(length(from) == length(to))
  if (length(weight) != length(from)) {
    stop("weight must be missing, length 1, or as long as from and to")
  } 
  
  xtabs(weight ~ to + from) |> 
    as.data.frame() |> 
    transform(Freq = Freq / sum(Freq), 
              fill = to) |> 
    ggalluvial::to_lodes_form(res, axes = c("from", "to")) |> 
    transform(stratum = ordered(stratum, risk_levels),
              fill = ordered(fill, risk_levels))
}

# Function improved from: https://doi.org/10.5281/zenodo.7121199
weighted.quantile <- function(x, w = 1, probs = seq(0, 1, 0.25),
                              na.rm = FALSE, names = TRUE) {
  
  if (any(probs > 1) | any(probs < 0)) stop("'probs' outside [0,1]")
  
  if (length(w) == 1) w <- rep(w, length(x))
  if (length(w) != length(x)) stop("w must have length 1 or be as long as x")
  
  not.missing <- !is.na(w) & !is.na(x)
  if (isTRUE(na.rm)) {
    if (!identical(is.na(x), is.na(w))) warning(paste(
      "One or more observations with missing x but not missing w",
      "or vice-versa! Using only observations with non-missing values",
      "for both."
    ))
    w <- w[not.missing]
    x <- x[not.missing]
  } else if (any(is.na(x) | is.na(w))) {
    warning("One or more missing values in x or w.")
    return(rep(NA_real_, length(x)))
  }
  
  w <- w[order(x)] / sum(w)
  x <- x[order(x)]
  
  if (length(x) == 1) {
    res <- rep(x, length(probs))
  } else {
    cum_w <- cumsum(w) - w * (1 - (seq_along(w) - 1) / (length(w) - 1))
    res <- approx(x = cum_w, y = x, xout = probs, ties = "ordered")$y
  }
  # Accessing non-exported functions is not great practice, generally,
  # but we want to ensure the labels are identical to those of
  # non-weighted quantiles.
  if (isTRUE(names)) names(res) <- stats:::format_perc(probs, "%")
  
  res
}


#' Gwet's generalized first-order agreement coefficient (AC1) 
#' but with survey weights
#'
#' @param a assessment of the first rater
#' @param b assessment of the second rater
#' @param w survey weights
#'
#' @return Gwet's AC1
#'
weighted.ac1 <- function(a, b, w = 1) {
  
  stopifnot(length(a) == length(b))
  if (length(w) == 1) w <- rep(w, length(a))
  stopifnot(length(a) == length(w))
  
  ok <- !is.na(a) & !is.na(b) & !is.na(w)
  if (any(!ok)) {
    warning("Removing observations with missing values")
    a <- a[ok]
    b <- b[ok]
    w <- w[ok]
  }
  
  if (is.factor(a) & is.factor(b)) {
    stopifnot(levels(a) == levels(b))
    ll <- levels(a)
  } else {
    stopifnot(is.factor(a) == is.factor(b))
    ll <- unique(c(a, b))
  }
  
  # Observed agreement
  pa <- weighted.mean(a == b, w)
  # Expected chance agreement by category
  pihat_k <- sapply(ll, \(x) weighted.mean(((a==x) + (b==x)) / 2, w))
  # Overall expected chance agreement
  pe <- sum(pihat_k * (1 - pihat_k)) / (length(ll) - 1)
  
  (pa - pe) / (1 - pe)
}


weighted.var <- function(x, w) {
  # Implementation rom MSWD_w in:
  # https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
  # Also see:
  # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Frequency_weights
  sum(w * (x - weighted.mean(x, w))^2) * sum(w) / (sum(w)^2 - sum(w^2))
}

weighted.sd <- function(x, w) sqrt(weighted.var(x, w))