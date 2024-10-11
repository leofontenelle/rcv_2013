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


# tabulate_fig3() ----

risk_levels <- c("Low", "Intermediate", "High")

d_t3 <- data.frame(
  from = ordered(sample(risk_levels, 1e3, replace = TRUE)),
  to = ordered(sample(risk_levels, 1e3, replace = TRUE)),
  weight = rlnorm(1e3)
)
res_t3 <- with(d_t3, tabulate_fig3(from, to, weight))

stopifnot(
  isTRUE(all.equal(2, sum(res_t3$Freq))),
  isTRUE(all.equal(1, sum(subset(res_t3, x == "from", select = Freq)))),
  isTRUE(all.equal(1, sum(subset(res_t3, x == "to", select = Freq))))
)
for (column in c("from", "to")) {
  if (!isTRUE(all.equal(
    tapply(d_t3$weight, d_t3[[column]], sum)[risk_levels] |> 
    prop.table(),
    with(subset(res_t3, x == column),
         tapply(Freq, stratum, sum))
  ))) {
    stop(sprintf("Subtotals don't match in column '%s'", column))
  }
}

# weighted.ac1() ----

if (!("irrCAC" %in% installed.packages())) {
  stop("Testing function weighted.ac1() requires installing package \"irrCAC\"")
}

df_wac1_bg <- irrCAC::cac.ben.gerry
# One missing value is "" and somehow made it into the levels list
df_wac1_bg$Gerry <- factor(df_wac1_bg$Gerry, levels = letters[1:5])
df_wac1_bg <- df_wac1_bg[complete.cases(df_wac1_bg), ]
stopifnot(identical(
  # Gwet's AC1 is rounded to the 5th digit
  round(weighted.ac1(df_wac1_bg$Ben, df_wac1_bg$Gerry), 5),
  irrCAC::gwet.ac1.raw(df_wac1_bg[, c("Ben", "Gerry")])$est$coeff.val
))

df_wac1_g1g2 <- irrCAC::cac.raw.g1g2
for (column in paste0("Rater", 1:4)) {
  df_wac1_g1g2[[column]] <- as.character(df_wac1_g1g2[[column]])
  missing <- !(df_wac1_g1g2[[column]] %in% letters)
  df_wac1_g1g2[[column]][missing] <- NA_character_
}
stopifnot(with(
  subset(df_wac1_g1g2, complete.cases(Rater3, Rater4)),
  identical(round(weighted.ac1(Rater3, Rater4), 5),
            irrCAC::gwet.ac1.raw(cbind(Rater3, Rater4))$est$coeff.val)))


# weighted.quantile() ----

# Tests reused from: https://doi.org/10.5281/zenodo.7121199
local({
  x <- rnorm(100)
  stopifnot(isTRUE(all.equal(weighted.quantile(x, w = 1), quantile(x))))
})
local({
  # From example('weighted.mean')
  x <- c(3.7, 3.3, 3.5, 2.8)
  w <- c(5,   5,   4,   1)/15
  stopifnot(isTRUE(all.equal(
    weighted.quantile(x, w, 0:4/4, names = FALSE),
    c(2.8, 3.33611111111111, 3.46111111111111, 3.58157894736842,
      3.7)
  )))
})