# Agreement between different scores for the estimation of 
# cardiovascular risk in Brazil, 2013
#

# Load packages and functions ----

# https://github.com/boyercb/globorisk
load("contrib/sysdata.rda") # coefs, cvdr, rf
source("contrib/globorisk.R")

library(CVrisk)
library(ggalluvial)
library(ggplot2)
library(gridExtra)
library(PooledCohort)
library(scales)

source("functions.R")

scores <- c(framingham = "Framingham",
            pooled_cohort = "Pooled Cohort Equations",
            globorisk = "Globorisk-LAC")
risk_levels <- c("Baixo", "Intermediário", "Alto")


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
stopifnot(all(names(scores %in% names(d))))

for (score_name in names(scores)) {
  new_name <- sprintf("%s_cat", score_name)
  d[[new_name]] <- d[[score_name]] |> 
    cut(breaks = c(0, 0.1, 0.2, Inf), labels = risk_levels,
        right = FALSE, ordered_result = TRUE)
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
)}) (d[d$ok, names(scores)], d$survey_weight[d$ok])

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
  A = scores[c("framingham", "framingham", "pooled_cohort")],
  B = scores[c("globorisk", "pooled_cohort", "globorisk")],
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

## Table 1 ----

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

## Figure 1 ----

d_fig1 <- d |> 
  subset(subset = ok, select = c(names(scores), "survey_weight")) |> 
  reshape(direction = "long", 
          varying = list(names(scores)),
          timevar = "score",
          times = scores)
# Make sure the scores show up in the intended order
d_fig1$score <- ordered(d_fig1$score, scores)

fig1 <- d_fig1 |> 
  ggplot(aes(x = framingham, # the variable name is picked by reshape()
             color = score, 
             weight = survey_weight)) + 
  geom_vline(xintercept = c(0.1, 0.2), col = "gray25", lty = 2) + 
  geom_density(bw = 0.05, lwd = 1) +
  scale_x_continuous(NULL, 
                     breaks = 0:10/10, minor_breaks = NULL, 
                     labels = label_percent()) + 
  scale_y_continuous(NULL, labels = NULL) + 
  scale_color_manual(name = "Escore", values = palette("Okabe-Ito")) + 
  theme_light() + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1))

## Figure 2 ----

ratios_fig2 <- c(ratio_fg = "Framingham \uf7 Globorisk-LAC", 
            ratio_fp = "Framingham \uf7 Pooled Cohort Equations", 
            ratio_pg = "Pooled Cohort Equations \uf7 Globorisk-LAC")
d_fig2 <- d |> 
  subset(subset = ok, select = c(names(ratios_fig2), "survey_weight")) |> 
  reshape(direction = "long", 
          varying = list(names(ratios_fig2)),
          timevar = "ratio",
          times = ratios_fig2)
# Make sure the scores show up in the intended order
d_fig2$ratio <- ordered(d_fig2$ratio, ratios_fig2)

fig2 <- d_fig2 |> 
  ggplot(aes(x = ratio_fg, # the variable name is picked by reshape()
             color = ratio, 
             weight = survey_weight)) + 
  scale_x_continuous(name = NULL, 
                     breaks = 1.25/(0.8 ^ (2 * -4:4)),
                     minor_breaks = NULL,
                     labels = label_number(0.01, decimal.mark = ","), 
                     trans = "log") + 
  scale_y_continuous(NULL, labels = NULL) + 
  scale_color_brewer(name = "Raz\ue3o", type ="qual", palette = "Dark2") + 
  geom_rect(inherit.aes = FALSE, data = data.frame(xmin = 0.8, xmax = 1.25, ymin = -0.2, ymax = 1.4), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightgray") + 
  geom_vline(xintercept = c(0.8, 1.25), col = "gray25", lty = 2) + 
  geom_density(bw = 0.1, lwd = 1) +
  coord_cartesian(ylim = c(0, 1.3)) + 
  theme_light() + 
  theme(legend.position = "inside", 
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0, 1))

## Figure 3 ----

plot_fig3 <- function(d_fig3, xlabfrom, xlabto, tag) {
  xlabels <-  setNames(c(xlabfrom, xlabto), c("from", "to"))
    
  res <- d_fig3 |> 
    # Reorder levels to make sure "High" is high in "Low" is low
    transform(fill = factor(fill, rev(levels(fill))), 
              stratum = factor(stratum, rev(levels(stratum)))) |> 
    ggplot(aes(x = x, y = Freq, alluvium = alluvium, stratum = stratum)) + 
    geom_flow(aes(color = stratum, fill = fill), alpha = 0.5) + 
    geom_stratum(aes(color = stratum, fill = stratum), alpha = 0.5) + 
    geom_text(aes(label = stratum), stat = "stratum") +
    scale_x_discrete(NULL, labels = xlabels, expand = c(0.15, 0.15)) + 
    scale_y_continuous(NULL, labels = label_percent(accuracy = 1)) + 
    labs(tag = tag) +
    scale_fill_viridis_d(guide = NULL, option = "E") + 
    scale_color_viridis_d(guide = NULL, option = "E") + 
    theme_light()
  res
}

fig3a <- tabulate_fig3(d$framingham_cat, d$globorisk_cat, d$survey_weight) |> 
  plot_fig3(scores["framingham"], scores["globorisk"], tag = "A")
fig3b <- tabulate_fig3(d$framingham_cat, d$pooled_cohort_cat, d$survey_weight) |> 
  plot_fig3(scores["framingham"], scores["pooled_cohort"], tag = "B")
fig3c <- tabulate_fig3(d$pooled_cohort_cat, d$globorisk_cat, d$survey_weight) |> 
  plot_fig3(scores["pooled_cohort"], scores["globorisk"], tag = "C")

# Write output ----

write.csv2(tab1, "tab1.csv", row.names = FALSE, fileEncoding = "UTF-8")
ggsave(filename = "fig1.png", 
       plot = fig1 + theme(text = element_text(size = 20)), 
       width = 1500, height = 1500 / 2, units = "px", 
       dpi = 96)
ggsave(filename = "fig2.png", 
       plot = fig2 + theme(text = element_text(size = 20)), 
       width = 1500, height = 1500 / 2, units = "px", 
       dpi = 96)
ggsave(filename = "fig3a.png", 
       plot = fig3a + theme(text = element_text(size = 20)), 
       width = 1500 / 2, height = 1500 / 2, units = "px", 
       dpi = 96)
ggsave(filename = "fig3b.png", 
       plot = fig3b + theme(text = element_text(size = 20)), 
       width = 1500 / 2, height = 1500 / 2, units = "px", 
       dpi = 96)
ggsave(filename = "fig3c.png", 
       plot = fig3c + theme(text = element_text(size = 20)), 
       width = 1500 / 2, height = 1500 / 2, units = "px", 
       dpi = 96)
ggsave(filename = "fig3.png",
       plot = arrangeGrob(fig3a + theme(text = element_text(size = 20)), 
                          fig3b + theme(text = element_text(size = 20)), 
                          fig3c + theme(text = element_text(size = 20)),
                          nrow = 2, ncol = 2),
       width = 1500, height = 1500, units = "px",
       dpi = 96)
