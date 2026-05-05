# compute VOI metrics from MCMC draws
compute_voi_metrics <- function(x,
                                sigma_min = 2,
                                ess_min = 400,
                                compute_evppi_prev = TRUE,
                                evppi_method = "sal",
                                center_rows = NULL) {

  # ---- get draws ----
  draws <- if (inherits(x, "mcmc.list")) {
    mcmc_list_to_draws(x)
  } else if (is.data.frame(x)) {
    tibble::as_tibble(x)
  } else {
    stop("`x` must be an mcmc.list or a data.frame of draws.", call. = FALSE)
  }

  need <- c("NBnew", "NBnew_TA", "ENBnew", "ENBnew_TA", "prevnew")
  missing <- setdiff(need, names(draws))
  if (length(missing) > 0) message("Missing draws: ", paste(missing, collapse=", "))
  # For prevalence EVPPI, prevnew is needed.
  # For population EVPI, ENBnew* are needed.

  # ---- baseline: current info ----
  if (!all(c("NBnew", "NBnew_TA") %in% names(draws))) {
    stop("Draws must contain NBnew and NBnew_TA to compute VOI.", call. = FALSE)
  }
  NBnew    <- draws$NBnew
  NBnew_TA <- draws$NBnew_TA

  NB_current <- max(0, mean(NBnew), mean(NBnew_TA))

  # ---- difference model vs treat-all (if pooled nodes available) ----
  has_pooled <- all(c("pooledNB", "pooledNB_TA") %in% names(draws))

  diff_pooled_mean <- NA_real_
  diff_pooled_low  <- NA_real_
  diff_pooled_high <- NA_real_

  if (has_pooled) {
    d_pool <- draws$pooledNB - draws$pooledNB_TA
    diff_pooled_mean <- mean(d_pool)
    qs <- stats::quantile(d_pool, c(0.025, 0.975), names = FALSE)
    diff_pooled_low  <- qs[1]
    diff_pooled_high <- qs[2]
  }

  # predictive/new-setting difference is always available from NBnew/NBnew_TA
  d_new <- NBnew - NBnew_TA
  diff_new_low  <- stats::quantile(d_new, 0.025, names = FALSE)
  diff_new_high <- stats::quantile(d_new, 0.975, names = FALSE)

  diff_modelvsTA_label <- if (is.finite(diff_pooled_mean)) {
    sprintf(
      "Pooled (mean [2.5%%, 97.5%%]): %.3f [%.3f, %.3f]; New ([2.5%%, 97.5%%]): [%.3f, %.3f]",
      diff_pooled_mean, diff_pooled_low, diff_pooled_high,
      diff_new_low, diff_new_high
    )
  } else {
    sprintf("New ([2.5%%, 97.5%%]): [%.3f, %.3f]", diff_new_low, diff_new_high)
  }

  # ---- determine winner strategy under current information ----
  m_model <- mean(NBnew)
  m_ta    <- mean(NBnew_TA)

  means <- c(treat_none = 0, model = m_model, treat_all = m_ta)
  winner_strategy <- names(means)[which.max(means)]

  winner_NB <- switch(
    winner_strategy,
    treat_none = rep(0, length(NBnew)),
    model      = NBnew,
    treat_all  = NBnew_TA
  )

  # ---- cluster perfect info ----
  tmp_cluster <- pmax(NBnew, NBnew_TA, 0)
  NB_cluster_pi <- mean(tmp_cluster)
  EVPI_cluster  <- max(0, NB_cluster_pi - NB_current) # protection against negative EVPI due to MC jitter

  # ---- population perfect info (if available) ----
  has_pop <- all(c("ENBnew", "ENBnew_TA") %in% names(draws))
  if (has_pop) {
    ENBnew    <- draws$ENBnew
    ENBnew_TA <- draws$ENBnew_TA
    tmp_pop <- pmax(ENBnew, ENBnew_TA, 0)
    NB_pop_pi <- mean(tmp_pop)
    EVPI_pop  <- max(0, NB_pop_pi - NB_current)
  } else {
    NB_pop_pi <- NA_real_
    EVPI_pop  <- NA_real_
  }

  # ---- EVPPI prevalence (if requested and available) ----
  has_prev <- "prevnew" %in% names(draws)
  NB_cluster_partialprev_pi <- NA_real_
  EVPPI_prev <- NA_real_

  if (compute_evppi_prev) {
    if (!has_prev) {
      warning("compute_evppi_prev=TRUE but prevnew not found in draws; skipping EVPPI(prevalence).")
    } else {
      if (!requireNamespace("voi", quietly = TRUE)) {
        warning("Package 'voi' not installed; skipping EVPPI(prevalence).")
      } else {
        NB_mat <- cbind(NBnew, NBnew_TA, rep(0, length(NBnew)))
        evppi_obj <- voi::evppi(NB_mat, inputs = draws$prevnew,
                                method = evppi_method, se = TRUE)
        EVPPI_prev <- max(0, as.numeric(evppi_obj$evppi))  # protection against negative EVPPI due to MC jitter
        NB_cluster_partialprev_pi <- NB_current + EVPPI_prev # for presentation symmetry with other blocks
      }
    }
  }

  # ---- center-specific EVPI for selected study rows (if available) ----
  center_metrics <- NULL

  if (!is.null(center_rows)) {
    center_rows <- unique(as.integer(center_rows))
    center_rows <- center_rows[!is.na(center_rows)]

    pick_col <- function(nm1, nm2) {
      if (nm1 %in% names(draws)) nm1 else if (nm2 %in% names(draws)) nm2 else NA_character_
    }

    center_res <- lapply(center_rows, function(cr) {
      nb_name <- sprintf("NB[%d]", cr)
      ta_name <- sprintf("NB_TA[%d]", cr)

      # allow for possible name mangling
      nb_col <- pick_col(nb_name, make.names(nb_name))
      ta_col <- pick_col(ta_name, make.names(ta_name))

      if (is.na(nb_col) || is.na(ta_col)) {
        return(tibble::tibble(
          center_row = cr,
          NB_center_currentinfo = NA_real_,
          NB_center_perfectinfo = NA_real_,
          EVPI_center_perfectinfo = NA_real_,
          center_winner_strategy = NA_character_,
          sigma_center_strategy = NA_real_,
          sigma_center_evpi = NA_real_,
          center_nodes_found = FALSE
        ))
      }

      NBc <- draws[[nb_col]]
      TAc <- draws[[ta_col]]

      NB_center_currentinfo <- max(0, mean(NBc), mean(TAc))
      NB_center_perfectinfo <- mean(pmax(NBc, TAc, 0))
      EVPI_center <- max(0, NB_center_perfectinfo - NB_center_currentinfo)

      means_c <- c(treat_none = 0, model = mean(NBc), treat_all = mean(TAc))
      center_strategy <- names(means_c)[which.max(means_c)]

      sig0 <- calc_sigma(cbind(NBc, TAc, 0))
      sigma_center_strategy <- if (length(sig0) == 0) NA_real_ else max(sig0, na.rm = TRUE)

      tmp_c <- pmax(NBc, TAc, 0)
      winnerNB_c <- switch(center_strategy,
                           treat_none = rep(0, length(NBc)),
                           model = NBc,
                           treat_all = TAc)

      sig_ev <- calc_sigma(cbind(tmp_c, winnerNB_c))
      sigma_center_evpi <- if (length(sig_ev) == 0) NA_real_ else max(sig_ev, na.rm = TRUE)

      tibble::tibble(
        center_row = cr,
        NB_model_mean = mean(NBc),
        NB_TA_mean = mean(TAc),
        NB_center_currentinfo = NB_center_currentinfo,
        NB_center_perfectinfo = NB_center_perfectinfo,
        EVPI_center_perfectinfo = EVPI_center,
        center_winner_strategy = center_strategy,
        sigma_center_strategy = sigma_center_strategy,
        sigma_center_evpi = sigma_center_evpi,
        center_nodes_found = TRUE
      )
    })

    center_metrics <- dplyr::bind_rows(center_res)

    if (any(!center_metrics$center_nodes_found)) {
      missing_rows <- center_metrics$center_row[!center_metrics$center_nodes_found]
      warning(
        "Some requested center rows were not found in the monitored draws: ",
        paste(missing_rows, collapse = ", "),
        ". Make sure those nodes were monitored in sampling.",
        call. = FALSE
      )
    }
  }

  # ---- perfect info never changes decision guards ----
  if (sum(tmp_cluster > winner_NB) == 0) {
    EVPI_cluster <- 0
  }

  if (has_pop && sum(tmp_pop > ENBnew) == 0) {
    EVPI_pop <- 0
  }

  # ---- MC diagnostics ----
  diagnostics <- compute_voi_diagnostics(
    x              = x,
    draws          = draws,
    NBnew          = NBnew,
    NBnew_TA       = NBnew_TA,
    winner_NB      = winner_NB,
    tmp_cluster    = tmp_cluster,
    tmp_pop        = if (has_pop) tmp_pop else NULL,
    ENBnew         = if (has_pop) ENBnew  else NULL,
    center_metrics = center_metrics,
    center_rows    = center_rows,
    sigma_min      = sigma_min,
    ess_min        = ess_min
  )

  # ---- output metrics ----
  metrics <- tibble::tibble(
    winner_strategy = winner_strategy,
    NB_currentinfo = NB_current,
    NB_cluster_perfectinfo = NB_cluster_pi,
    EVPI_cluster_perfectinfo = EVPI_cluster,
    NB_population_perfectinfo = NB_pop_pi,
    EVPI_population_perfectinfo = EVPI_pop,
    NB_cluster_partialprev_pi = NB_cluster_partialprev_pi,
    EVPPI_cluster_perfectprevalenceinfo = EVPPI_prev,
    diff_modelvsTA = diff_modelvsTA_label
  )

  list(metrics = metrics, diagnostics = diagnostics, center_metrics = center_metrics)
}

compute_voi_diagnostics <- function(
    x,
    draws,
    NBnew,
    NBnew_TA,
    winner_NB,
    tmp_cluster,
    tmp_pop = NULL,
    ENBnew = NULL,
    center_metrics = NULL,
    center_rows = NULL,
    sigma_min = 2,
    ess_min = 400
) {
  has_pop <- !is.null(tmp_pop)

  # ---- sigma calculations ----

  # sigma for differences among strategies at the new setting level
  sigma_strategy <- {
    sig <- calc_sigma(cbind(NBnew, NBnew_TA, 0))
    if (length(sig) == 0) NA_real_ else max(sig, na.rm = TRUE)
  }

  # sigma for the EVPI difference itself: tmp_cluster - winner_NB
  sigma_evpi_cluster <- {
    sig <- calc_sigma(cbind(tmp_cluster, winner_NB))
    if (length(sig) == 0) NA_real_ else max(sig, na.rm = TRUE)
  }

  # sigma for population EVPI difference, if available
  sigma_evpi_pop <- if (has_pop) {
    sig <- calc_sigma(cbind(tmp_pop, winner_NB))
    if (length(sig) == 0) NA_real_ else max(sig, na.rm = TRUE)
  } else NA_real_

  # ---- perfect info never changes decision guards ----
  perfect_info_never_beats_model <- sum(tmp_cluster > winner_NB) == 0
  if (perfect_info_never_beats_model) {
    sigma_evpi_cluster <- Inf
  }

  perfect_info_never_beats_model_pop <- has_pop && sum(tmp_pop > ENBnew) == 0
  if (perfect_info_never_beats_model_pop) {
    sigma_evpi_pop <- Inf
  }

  # ---- flags ----

  # we do more sampling if the monte carlo error on the difference between two strategies is smaller than 2
  # or if the monte carlo error on the random effects summary measures is more than 5% of the standard deviation of their posterior distribution, i.e. MCSE/posterior sd > 0.05
  flag_sigma      <- !is.na(sigma_strategy) && sigma_strategy < sigma_min
  flag_sigma_evpi <- !is.na(sigma_evpi_cluster) && sigma_evpi_cluster < sigma_min
  flag_sigma_evpi_pop <- if (has_pop) {
    !is.na(sigma_evpi_pop) && sigma_evpi_pop < sigma_min
  } else FALSE

  flag_sigma_center_strategy <- FALSE
  flag_sigma_center_evpi     <- FALSE
  if (!is.null(center_metrics)) {
    flag_sigma_center_strategy <- any(center_metrics$sigma_center_strategy < sigma_min, na.rm = TRUE)
    flag_sigma_center_evpi     <- any(center_metrics$sigma_center_evpi < sigma_min, na.rm = TRUE)
  }

  # ---- ESS ----
  # ESS-based precision checks for key summary nodes
  # target relative MCSE (MCSE / posterior SD)
  # MCSE /posterior sd = (sd/sqrt(ESS))/sd = 1/sqrt(ESS) < 0.05 -> ESS > 400
  # ess_min <- (1 / rmcse_max)^2   # rmcse_max <- 0.05, equivalently require ESS > 400

  ess_pooledNB <- ess_pooledNB_TA <- ess_probuseful <- NA_real_
  flag_ess <- FALSE
  if (inherits(x, "mcmc.list")) {
    ess <- coda::effectiveSize(x)
    ess_pooledNB    <- unname(ess["pooledNB"])
    ess_pooledNB_TA <- unname(ess["pooledNB_TA"])
    ess_probuseful  <- unname(ess["probuseful"])
    flag_ess <- any(
      c(ess_pooledNB, ess_pooledNB_TA, ess_probuseful) < ess_min,
      na.rm = TRUE
    )
  }

  # ---- diagnostics tibble ----
  tibble::tibble(
    n_draws             = nrow(draws),
    sigma_strategy      = sigma_strategy,
    sigma_evpi_cluster  = sigma_evpi_cluster,
    sigma_evpi_pop      = sigma_evpi_pop,
    sigma_min           = sigma_min,
    ess_min             = ess_min,
    ess_pooledNB        = ess_pooledNB,
    ess_pooledNB_TA     = ess_pooledNB_TA,
    ess_probuseful      = ess_probuseful,
    min_sigma_center_strategy = if (!is.null(center_metrics)) min(center_metrics$sigma_center_strategy, na.rm = TRUE) else NA_real_,
    min_sigma_center_evpi     = if (!is.null(center_metrics)) min(center_metrics$sigma_center_evpi,     na.rm = TRUE) else NA_real_,
    n_centers_requested = if (is.null(center_rows)) 0L else length(center_rows),
    needs_more_sampling = (
      flag_sigma | flag_sigma_evpi | flag_sigma_evpi_pop |
        flag_sigma_center_strategy | flag_sigma_center_evpi | flag_ess
    )
  )
}

