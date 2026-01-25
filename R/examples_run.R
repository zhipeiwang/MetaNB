# setwd("C:/Users/P70096364/OneDrive - Maastricht University/MetaNB/R")

lapply(c("MA_NB_tri.R", "model_text_weak.R", "model_text_wishart.R", "helpers.R", "plot_forest.R"), source)
load("data_examples.RData")

# load libraries
library(rjags)
library(coda)
library(tidyverse)
library(forestploter)

MA_NBCA125_weak <- MA_NB_tri(data_ADNEXCA125, prior = "weak", t = 0.1, return_vars = c("NB", "RU", "probuseful", "NBnew", "NBnew_TA"))

attributes(MA_NBCA125_weak)

summary(MA_NBCA125_weak)

sum <- summarize_mcmc_outputs(
  MA_NBCA125_weak,
  study_info = study_info_ADNEXCA125,
  targets = c("NB", "RU", "probuseful"),
  targets_per_study = c("NB", "RU"),
  return_ref = TRUE
)
sum


plot_forest_metric_forestploter(
  sum,
  metric = "NB",
#  center = "Median",
  xlim = c(-0.1, 0.7),
  file_png = "forest_nb.png",
  file_pdf = "forest_nb.pdf"
)

plot_forest_metric_forestploter(
  sum,
  metric = "RU",
  xlim = c(-0.1, 1),
  file_png = "forest_ru.png",
  file_pdf = "forest_ru.pdf"
)

# Wishart prior
MA_NBCA125_wishart <- MA_NB_tri(data_ADNEXCA125, t = 0.1, prior = "wishart")

attributes(MA_NBCA125_wishart)

summary(MA_NBCA125_wishart)


# ----------------- Stats Med ------------------------------
# Fever example, weak realistic prior
MA_fever_weak <- MA_NB_tri(data_fever, prior = "weak", t = 0.2,
                               return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                               "NBnew_ref","probharmful","probharmful_ref"))
summary(MA_fever_weak)

# Fever example, wishart prior
MA_fever_wishart <- MA_NB_tri(data_fever, prior = "wishart", t = 0.2,
                                  return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                  "NBnew_ref","probharmful","probharmful_ref"))
summary(MA_fever_wishart)
# Ovarian cancer example, weak realistic prior
MA_ovarian_t005_weak <- MA_NB_tri(data_ovarian_t005, prior = "weak", t = 0.05, prev_ref = 0.15,
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
                                                      "NBnew_ref","probharmful","probharmful_ref"))
summary(MA_ovarian_t005_weak)

MA_ovarian_t01_weak <- MA_NB_tri(data_ovarian_t01, prior = "weak", t = 0.1, prev_ref = 0.15,
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
                                                     "NBnew_ref","probharmful","probharmful_ref"))
summary(MA_ovarian_t01_weak)

# Ovarian cancer example, wishart prior
MA_ovarian_t005_wishart <- MA_NB_tri(data_ovarian_t005, prior = "wishart", t = 0.05, prev_ref = 0.15,
                                         return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                         "NBnew_ref","probharmful","probharmful_ref"))
summary(MA_ovarian_t005_wishart)

MA_ovarian_t01_wishart <- MA_NB_tri(data_ovarian_t01, prior = "wishart", t = 0.1, prev_ref = 0.15,
                                        return_vars = c("pooledsens","pooledspec","pooledNB","pooledNB_ref","pooledNB_TA","pooledNB_TA_ref","NBnew",
                                                        "NBnew_ref","probharmful","probharmful_ref"))
summary(MA_ovarian_t01_wishart)

#################### probuseful definition ###########################
# extract posterior draws into a data frame
draws <- as.data.frame(do.call(rbind, MA_NBCA125_weak))
# define the best alternative per draw
best_alt <- pmax(draws$NBnew_TA, 0)

# at least as good as best alternative (tie-allowed)
as.numeric(draws$NBnew >= best_alt) %>% mean()
# strictly better than best alternative
as.numeric(draws$NBnew > best_alt) %>% mean()
# How often do ties actually happen?
mean(draws$NBnew == best_alt)
