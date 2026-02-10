# setwd("C:/Users/P70096364/OneDrive - Maastricht University/MetaNB/R")

lapply(c("MA_NB_tri.R", "model_text_weak.R", "model_text_wishart.R", "helpers.R", "plot_forest.R", "MA_NB_tri_voi.R"), source)
load("data_examples.RData")

# load libraries
library(rjags)
library(coda)
library(tidyverse)
library(forestploter)
library(rlang)

MA_NBCA125_weak <- MA_NB_tri(data = data_ADNEXCA125_full, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                             prior_type = "weak", t = 0.1,
                             return_vars = c("NB", "RU", "probuseful", "NBnew", "NBnew_TA", "sens", "spec"),
                             return_ref = TRUE, return_fit = TRUE)


summary(MA_NBCA125_weak$samples)

sum_user <- summarize_for_users(
  MA_NBCA125_weak,
  data = data_ADNEXCA125_full,
  label_cols = c("Publication", "Country", "N", "Prev"),
  metrics = c("NB", "RU", "probuseful", "sens", "spec"),
  per_study_metrics = c("NB", "RU", "sens", "spec"),
  return_ref = FALSE,
  include_per_study = FALSE)
sum_user

sum <- summarize_for_forestplot(
  MA_NBCA125_weak,
  data = data_ADNEXCA125_full,
  label_cols = c("Publication", "Country", "N", "Prev"),
  targets = c("NB", "RU", "probuseful", "sens", "spec"),
  targets_per_study = c("NB", "RU", "sens", "spec"),
  return_ref = TRUE
)
sum


plot_forest_metric_forestploter(
  sum,
  label_cols = c("Publication", "Country", "N", "Prevalence"),
  prev_col = "Prevalence",
  metric = "NB",
#  center = "Median",
  xlim = c(-0.1, 0.7),
  file_png = "forest_nb.png",
  file_pdf = "forest_nb.pdf"
)

plot_forest_metric_forestploter(
  sum,
  label_cols = c("Publication", "Country", "N", "Prev"),
  prev_col = "Prev",
  metric = "RU",
  xlim = c(-0.1, 1),
  file_png = "forest_ru.png",
  file_pdf = "forest_ru.pdf"
)

plot_forest_metric_forestploter(
  sum,
  label_cols = c("Publication", "Country", "N", "Prev"),
  prev_col = "Prev",
  metric = "sens",
  # center = "Median",
  xlim = c(0.7, 1),
  xticks = seq(0.7, 1, by = 0.05),
  file_png = "forest_sens.png",
  file_pdf = "forest_sens.pdf"
)


plot_forest_metric_forestploter(
  sum,
  label_cols = c("Publication", "Country", "N", "Prev"),
  prev_col = "Prev",
  metric = "spec",
  # center = "Median",
  xlim = c(0.3, 1),
  xticks = seq(0.3, 1, by = 0.1),
  file_png = "forest_spec.png",
  file_pdf = "forest_spec.pdf"
)

# Wishart prior
MA_NBCA125_wishart <- MA_NB_tri(data_ADNEXCA125_full, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                                t = 0.1, prior_type = "wishart")

attributes(MA_NBCA125_wishart)

summary(MA_NBCA125_wishart)


# ----------------- Stats Med ------------------------------
# Fever example, weak realistic prior
MA_fever_weak <- MA_NB_tri(data_fever, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                           prior_type = "weak", t = 0.2,
                               return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                               "NBnew_ref","probuseful","probuseful_ref"))
summary(MA_fever_weak)

# Fever example, wishart prior
MA_fever_wishart <- MA_NB_tri(data_fever, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                              prior_type = "wishart", t = 0.2,
                                  return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                  "NBnew_ref","probuseful","probuseful_ref"))
summary(MA_fever_wishart)
# Ovarian cancer example, weak realistic prior
MA_ovarian_t005_weak <- MA_NB_tri(data_ovarian_t005, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                                  prior_type = "weak", t = 0.05, prev_ref = 0.15,
                                      weak_priors = list(
                                        mu_etap=-1.04, tau_etap=1/(0.23^2),
                                        mu_lambdasens0=2.99, tau_lambdasens0=1/(0.36^2),
                                        mu_lambdaspec0=1.37, tau_lambdaspec0=1/(1.73^2),
                                        mu_zss=0.07, tau_zss=1/(0.49^2),
                                        mu_zsp=-0.63, tau_zsp=1/(0.33^2),
                                        a_csp=-0.41, b_csp=0.41,
                                        tau_varprev=0.25, tau_varsens=0.25, tau_varspec=0.25 # not changed, still half-normal
                                      ), # table A1 in stats med paper
                                      return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                      "NBnew_ref","probuseful","probuseful_ref"))
summary(MA_ovarian_t005_weak)

MA_ovarian_t01_weak <- MA_NB_tri(data_ovarian_t01, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                                 prior_type = "weak", t = 0.1, prev_ref = 0.15,
                                     weak_priors = list(
                                       mu_etap=-1.04, tau_etap=1/(0.22^2),
                                       mu_lambdasens0=2.31, tau_lambdasens0=1/(0.31^2),
                                       mu_lambdaspec0=2.12, tau_lambdaspec0=1/(2.26^2),
                                       mu_zss=-0.16, tau_zss=1/(0.42^2),
                                       mu_zsp=-0.74, tau_zsp=1/(0.33^2),
                                       a_csp=-0.23, b_csp=0.38,
                                       tau_varprev=0.25, tau_varsens=0.25, tau_varspec=0.25 # not changed, still half-normal
                                     ), # table A1 in stats med paper
                                     return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                     "NBnew_ref","probuseful","probuseful_ref"))
summary(MA_ovarian_t01_weak)

# Ovarian cancer example, wishart prior
MA_ovarian_t005_wishart <- MA_NB_tri(data_ovarian_t005, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                                     prior_type = "wishart", t = 0.05, prev_ref = 0.15,
                                         return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                         "NBnew_ref","probuseful","probuseful_ref"))
summary(MA_ovarian_t005_wishart)

MA_ovarian_t01_wishart <- MA_NB_tri(data_ovarian_t01, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                                    prior_type = "wishart", t = 0.1, prev_ref = 0.15,
                                        return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                        "NBnew_ref","probuseful","probuseful_ref"))
summary(MA_ovarian_t01_wishart)


##################################  VOI  #####################################
MA_NBCA125_weak_forVOI <- MA_NB_tri(data = data_ADNEXCA125_full, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                             prior_type = "weak", t = 0.1,
                             return_vars = c("NBnew", "NBnew_TA", "ENBnew", "ENBnew_TA", "prevnew", "pooledNB", "pooledNB_TA", "probuseful"),
                             compute_EVPI = TRUE)

res <- compute_voi_metrics(MA_NBCA125_weak_forVOI)
res$metrics
res$diagnostics

fit_voi <- sample_voi_draws(
  data = data_ADNEXCA125_full,
  tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
  prior_type = "weak",
  t = 0.1,
  J = 1000,
  burnin = 3000,
  sample_by = 2000,
  auto_resample = TRUE,
  max_rounds = 20
)

fit_voi$metrics
fit_voi$diagnostics

voi_res <- MA_NB_tri_voi(
  data_ADNEXCA125_full,
  tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
  # center_row = 5,
  # center_label_cols = c("Publication", "Country"),
  prior_type = "weak",
  t = 0.1,
  auto_resample = TRUE
)

options(scipen = 999)
glimpse(voi_res$metrics)
c(voi_res$metrics$EVPI_population_perfectinfo,
  voi_res$metrics$EVPPI_cluster_perfectprevalenceinfo,voi_res$metrics$EVPI_cluster_perfectinfo)
glimpse(voi_res$diagnostics)



# Make a copy so we don't overwrite our original
dat_test <- data_ADNEXCA125_full

# Initialize fake CI columns with NA
dat_test$sens_ci_low  <- NA_real_
dat_test$sens_ci_high <- NA_real_
dat_test$spec_ci_low  <- NA_real_
dat_test$spec_ci_high <- NA_real_

# Choose some rows to "pretend" we have reported CIs for
# Pick spread-out indices so you can find them easily
idx_rep <- c(2, 7, 13, 21, 29, 37)

# Compute observed sens/spec (for sanity checks / potential point override)
obs_sens <- with(dat_test, tp / n_event)
obs_spec <- with(dat_test, tn / n_nonevent)

# Create fake reported CIs that are clearly different from your Bayesian CrIs:
# Identical, very narrow "reported" intervals
dat_test$sens_ci_low[idx_rep]  <- 0.910
dat_test$sens_ci_high[idx_rep] <- 0.915

dat_test$spec_ci_low[idx_rep]  <- 0.980
dat_test$spec_ci_high[idx_rep] <- 0.985

# Optional: quick sanity display (should show numbers only in idx_rep rows)
dat_test[idx_rep, c("sens_ci_low","sens_ci_high","spec_ci_low","spec_ci_high")]

plot_forest_metric_forestploter(
  sum,
  data = dat_test,
  label_cols = c("Publication", "Country", "N", "Prev"),
  prev_col = "Prev",
  metric = "sens",
  use_reported_ci = TRUE,
  reported_low_col = "sens_ci_low",
  reported_high_col = "sens_ci_high",
  # optional for easier viewing:
  xlim = c(0.7, 1),
  xticks = seq(0.7, 1, by = 0.05),
  file_png = "forest_sens.png",
  file_pdf = "forest_sens.pdf"
)


