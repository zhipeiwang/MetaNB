# setwd("C:/Users/P70096364/OneDrive - Maastricht University/MetaNB/dev")

line_inits <- list(list("etap"=-1.04,"lambdasens0"=2.31,"lambdaspec0"=2.12,"varprev"=0.82,"varsens"=0.46,"varspec"=0.30,"zss"=-0.23,"corr.sens.prev"=0,"zsp"=-0.74),
                   list("etap"=-0.79,"lambdasens0"=3.59, "lambdaspec0"=1.20, "varprev"=0.55, "varsens"=1.30, "varspec"=0.47, "zss"=-0.42, "corr.sens.prev"=0.54, "zsp"=-0.54))

MA_NBCA125_weak <- MA_NB_tri(data = data_ADNEXCA125, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                             prior_type = "weak", t = 0.1, seed = 123,
                             return_vars = c("NB", "RU", "probuseful", "NBnew", "NBnew_TA", "sens", "spec", "pooledprev", "prevnew"),
                             return_known = TRUE)


summary(MA_NBCA125_weak$samples)

sum_user <- summarize_tri_ma(
  MA_NBCA125_weak,
  data = data_ADNEXCA125,
  label_cols = c("Publication", "Country", "N", "Prev"),
  metrics = c("NB", "RU", "probuseful", "sens", "spec", "pooledprev", "prevnew"),
  per_study_metrics = c("NB", "RU", "sens", "spec"),
  return_known = TRUE,
  include_per_study = FALSE)
sum_user



plot_forest(
  MA_NBCA125_weak$samples,
  data = data_ADNEXCA125,
  label_cols = c("Publication", "Country", "N", "Prev"),
  metric = "NB",
#  center = "Median",
  t = 0.1,
  mark_imputed = FALSE,
  xlim = c(-0.1, 0.7),
  file_png = "forest_nb.png",
  file_pdf = "forest_nb.pdf"
)

plot_forest(
  MA_NBCA125_weak$samples,
  data = data_ADNEXCA125,
  label_cols = c("Publication", "Country", "N", "Prev"),
  metric = "RU",
  mark_imputed = FALSE,
  xlim = c(-0.1, 1),
  file_png = "forest_ru1.png",
  file_pdf = "forest_ru1.pdf"
)

plot_forest(
  MA_NBCA125_weak$samples,
  data = data_ADNEXCA125,
  label_cols = c("Publication", "Country", "N", "Prev"),
  metric = "sens",
  # center = "Median",
  mark_imputed = TRUE,
  reported_est_col = "sens_point",
  reported_low_col = "sens_ci_low",
  reported_high_col = "sens_ci_high",
  xlim = c(0.7, 1.01),
  xticks = seq(0.7, 1, by = 0.05),
  file_png = "forest_sens.png",
  file_pdf = "forest_sens.pdf"
)


plot_forest(
  MA_NBCA125_weak$samples,
  data = data_ADNEXCA125,
  label_cols = c("Publication", "Country", "N", "Prev"),
  metric = "spec",
  # center = "Median",
  mark_imputed = TRUE,
  reported_est_col = "spec_point",
  reported_low_col = "spec_ci_low",
  reported_high_col = "spec_ci_high",
  plot_col_width = 70,
  xlim = c(0.5, 1),
  xticks = seq(0.5, 1, by = 0.1),
  file_png = "forest_spec.png",
  file_pdf = "forest_spec.pdf"
)

# Wishart prior
MA_NBCA125_wishart <- MA_NB_tri(data_ADNEXCA125, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                                t = 0.1, prior_type = "wishart")

attributes(MA_NBCA125_wishart)

summary(MA_NBCA125_wishart)


##################################  VOI  #####################################
MA_NBCA125_weak_forVOI <- MA_NB_tri(data = data_ADNEXCA125, tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
                             prior_type = "weak", t = 0.1, seed = 123,
                             return_vars = c("NBnew", "NBnew_TA", "ENBnew", "ENBnew_TA", "prevnew", "pooledNB", "pooledNB_TA", "probuseful"),
                             compute_EVPI = TRUE)

res <- compute_voi_metrics(MA_NBCA125_weak_forVOI$samples)
res$metrics
res$diagnostics

-------------------------------------------
fit_voi <- sample_voi_draws(
  data = data_ADNEXCA125,
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
  data_ADNEXCA125,
  tp = tp, tn = tn, n_event = n_event, n_nonevent = n_nonevent,
  center_rows = 1:37,
  center_label_cols = c("Publication", "Country"),
  prior_type = "weak",
  t = 0.1,
  inits = line_inits,
  auto_resample = TRUE
)

options(scipen = 999)
glimpse(voi_res$metrics)
c(voi_res$metrics$EVPI_population_perfectinfo,
  voi_res$metrics$EVPPI_cluster_perfectprevalenceinfo,voi_res$metrics$EVPI_cluster_perfectinfo)
glimpse(voi_res$diagnostics)
voi_res$center_metrics


# Make a copy so we don't overwrite our original
dat_test <- data_ADNEXCA125

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

plot_forest(
  MA_NBCA125_weak$samples,
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

# data_ADNEXCA125_full_extended <- data_ADNEXCA125_full %>%
#   mutate(
#     sens_point = adnex_comple$data_MA.Senvalue,
#     sens_ci_low = adnex_comple$data_MA.Senlow_CI,
#     sens_ci_high = adnex_comple$data_MA.Senup_CI,
#     spec_point = adnex_comple$data_MA.Spevalue,
#     spec_ci_low = adnex_comple$data_MA.Spelow_CI,
#     spec_ci_high = adnex_comple$data_MA.Speup_CI
#   )
#
# data_ADNEXCA125_full_extended <-
#   data_ADNEXCA125_full_extended %>%
#   dplyr::mutate(
#     dplyr::across(
#       dplyr::all_of(c(
#         "sens_point", "sens_ci_low", "sens_ci_high",
#         "spec_point", "spec_ci_low", "spec_ci_high"
#       )),
#       ~ suppressWarnings(as.numeric(.x) / 100)
#     )
#   )
#
# data_ADNEXCA125_full <- data_ADNEXCA125_full_extended
# save(data_ADNEXCA125_full, study_info_ADNEXCA125, data_fever, data_ovarian_t005, data_ovarian_t01, file = "data_examples.RData")


bayesdca_nb <- function(
    TP, FP, TN, FN,
    threshold,
    prior_p1 = 1, prior_p2 = 1,
    prior_Se1 = 1, prior_Se2 = 1,
    prior_Sp1 = 1, prior_Sp2 = 1,
    n_draws = 4000,
    summary_probs = c(0.025, 0.975)
) {

  # basic quantities
  N <- TP + FP + TN + FN
  d <- TP + FN
  thresholds <- threshold
  w_t <- thresholds / (1 - thresholds)

  # posterior parameters
  post_shape1_p <- d + prior_p1
  post_shape2_p <- N - d + prior_p2

  post_shape1_Se <- TP + prior_Se1
  post_shape2_Se <- d - TP + prior_Se2   # = FN + prior_Se2

  post_shape1_Sp <- TN + prior_Sp1
  post_shape2_Sp <- N - d - TN + prior_Sp2   # = FP + prior_Sp2

  # posterior draws
  p <- rbeta(n_draws, shape1 = post_shape1_p, shape2 = post_shape2_p)
  se <- rbeta(n_draws, shape1 = post_shape1_Se, shape2 = post_shape2_Se)
  sp <- rbeta(n_draws, shape1 = post_shape1_Sp, shape2 = post_shape2_Sp)

  # net benefit draws
  model <- se * p - w_t * (1 - sp) * (1 - p)
  treat_all <- p - (1 - p) * w_t
  treat_none <- rep(0, n_draws)

  # summaries
  summarize_draws <- function(x, probs = summary_probs) {
    c(
      mean = mean(x),
      median = median(x),
      lower = unname(quantile(x, probs[1])),
      upper = unname(quantile(x, probs[2]))
    )
  }

  means <- c(
    treat_none = mean(treat_none),
    model = mean(model),
    treat_all = mean(treat_all)
  )

  winner_strategy <- names(means)[which.max(means)]

  # EVPI
  enb_current <- max(means)

  enb_perfect <- mean(
    pmax(treat_none, treat_all, model)
  )

  EVPI <- enb_perfect - enb_current

  list(
    inputs = list(
      TP = TP, FP = FP, TN = TN, FN = FN,
      N = N, d = d, threshold = threshold
    ),
    parameters = list(
      prevalence = c(shape1 = post_shape1_p, shape2 = post_shape2_p),
      sensitivity = c(shape1 = post_shape1_Se, shape2 = post_shape2_Se),
      specificity = c(shape1 = post_shape1_Sp, shape2 = post_shape2_Sp)
    ),
    draws = list(
      p = p,
      se = se,
      sp = sp,
      model = model,
      treat_all = treat_all,
      treat_none = treat_none
    ),
    summary = list(
      prevalence = summarize_draws(p),
      sensitivity = summarize_draws(se),
      specificity = summarize_draws(sp),
      model = summarize_draws(model),
      treat_all = summarize_draws(treat_all),
      treat_none = summarize_draws(treat_none),
      winner_strategy = winner_strategy,
      EVPI = EVPI,
      enb_current = enb_current,
      enb_perfect = enb_perfect
    )
  )
}


# single-setting results
evpi_results <- vector("list", nrow(data_ADNEXCA125))

for (i in seq_len(nrow(data_ADNEXCA125))) {

  fit <- bayesdca_nb(
    TP = data_ADNEXCA125$tp[i],
    FP = data_ADNEXCA125$n_nonevent[i] - data_ADNEXCA125$tn[i],
    TN = data_ADNEXCA125$tn[i],
    FN = data_ADNEXCA125$n_event[i] - data_ADNEXCA125$tp[i],
    threshold = 0.1
  )

  evpi_results[[i]] <- data.frame(
    center_row = i,
    Publication = data_ADNEXCA125$Publication[i],
    Country = data_ADNEXCA125$Country[i],
    N = data_ADNEXCA125$N[i],
    Prev = data_ADNEXCA125$Prev[i],
    method = "single_setting",
    NB_center_currentinfo = fit$summary$enb_current,
    NB_center_perfectinfo = fit$summary$enb_perfect,
    center_winner_strategy = fit$summary$winner_strategy,
    EVPI_center_perfectinfo = fit$summary$EVPI
  )
}

evpi_single_df <- do.call(rbind, evpi_results)

# triMA results
evpi_triMA_df <- voi_res$center_metrics %>%
  dplyr::left_join(
    data_ADNEXCA125 %>%
      dplyr::mutate(center_row = dplyr::row_number()) %>%
      dplyr::select(center_row, Publication, Country, N, Prev),
    by = "center_row"
  ) %>%
  dplyr::mutate(method = "triMA") %>%
  dplyr::select(
    center_row, Publication, Country, N, Prev, method,
    NB_center_currentinfo,
    NB_center_perfectinfo,
    center_winner_strategy,
    EVPI_center_perfectinfo
  ) %>%
  data.frame()

# combine
evpi_long <- dplyr::bind_rows(evpi_single_df, evpi_triMA_df)

evpi_wide <- evpi_long %>%
  dplyr::select(
    center_row, Publication, Country, N, Prev, method,
    EVPI_center_perfectinfo,
    center_winner_strategy,
    NB_center_currentinfo,
    NB_center_perfectinfo
  ) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = c(
      EVPI_center_perfectinfo,
      center_winner_strategy,
      NB_center_currentinfo,
      NB_center_perfectinfo
    ),
    names_glue = "{.value}_{method}"
  ) %>%
  dplyr::mutate(
    EVPI_diff = EVPI_center_perfectinfo_triMA - EVPI_center_perfectinfo_single_setting,
    EVPI_diff_perc = 100 * EVPI_diff / EVPI_center_perfectinfo_single_setting,
    winner_same = center_winner_strategy_triMA == center_winner_strategy_single_setting
  )

ggplot(evpi_long, aes(x = Prev, y = EVPI_center_perfectinfo, color = center_winner_strategy)) +
  geom_point() +
  facet_wrap(~ method) +
  theme_minimal()

ggplot(evpi_wide, aes(
  x = EVPI_center_perfectinfo_single_setting,
  y = EVPI_center_perfectinfo_triMA,
  color = winner_same,
  size = N
)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Comparison of center-level EVPI estimates",
    x = "EVPI single-setting",
    y = "EVPI trivariate meta-analysis"
  ) +
  theme_minimal()

ggplot(evpi_wide, aes(
  x = EVPI_center_perfectinfo_single_setting,
  y = EVPI_center_perfectinfo_triMA,
  color = winner_same,
  size = as.numeric(sub("%", "", Prev))
)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Comparison of center-level EVPI estimates",
    x = "EVPI single-setting",
    y = "EVPI trivariate meta-analysis"
  ) +
  theme_minimal()

evpi_wide %>%
  filter(winner_same == FALSE)
