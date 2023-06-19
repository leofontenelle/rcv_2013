source("functions.R")

# between() ----

df_b <- list(x = rnorm(100), left = rnorm(1), right = rnorm(1))
stopifnot(all(identical(
  between(df_b$x, df_b$left, df_b$right),
  df_b$x >= df_b$left & df_b$x <= df_b$right
)))

# diagnose_diabetes() ----

df_dd <- expand.grid(Q029 = c(1:6, NA), 
                     Q030 = c(1:3, NA), 
                     Z034 = c(6.4, 6.5, 6.6, NA),
                     diabetes = NA_integer_)
df_dd$diabetes <- c(
  1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
  0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
  1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
  0L, 0L, 0L, NA, NA, NA, NA, NA, 0L, NA)

stopifnot(with(df_dd, identical(
  diagnose_diabetes(Q029, Q030, Z034), 
  diabetes
)))


# impute_bp_meds() ----

df_ibm <- expand.grid(Q001 = c(1:6, NA_integer_), 
                      Q002 = c(1:3, NA_integer_), 
                      Q006 = c(1:2, NA_real_),
                      bp_meds = NA_integer_)
df_ibm$bp_meds <- c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, NA, NA, NA, 0, NA, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, NA, NA, NA, 0, NA
)

stopifnot(with(df_ibm, identical(
  impute_bp_meds(Q001, Q002, Q006), 
  bp_meds
)))

# weighted.quantile() ----

# Tests reused from: https://doi.org/10.5281/zenodo.7121199
local({
  x <- rnorm(100)
  stopifnot(isTRUE(all.equal(weighted.quantile(x, w = 1), quantile(x))))
})
local({
  # Example from man('weighted.mean')
  x <- c(3.7, 3.3, 3.5, 2.8)
  w <- c(5,   5,   4,   1)/15
  stopifnot(isTRUE(all.equal(
    weighted.quantile(x, w, 0:4/4, names = FALSE),
    c(2.8, 3.33611111111111, 3.46111111111111, 3.58157894736842,
      3.7)
  )))
})