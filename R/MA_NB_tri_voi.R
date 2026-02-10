# compute VOI metrics from MCMC draws
compute_voi_metrics <- function(x,
                                sigma_min = 2,
                                rmcse_max = 0.05,
                                compute_evppi_prev = TRUE,
                                evppi_method = "sal",
                                center_row = NULL) {

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

  # ---- center-specific EVPI for an existing study row (if available) ----
  sigma_center_strategy <- NA_real_
  sigma_center_evpi     <- NA_real_
  EVPI_center <- NA_real_
  NB_center_currentinfo <- NA_real_
  NB_center_perfectinfo <- NA_real_
  center_strategy <- NA_character_

  if (!is.null(center_row)) {
    nb_name <- sprintf("NB[%d]", center_row)
    ta_name <- sprintf("NB_TA[%d]", center_row)

    # allow for potential name mangling just in case
    nb_name2 <- make.names(nb_name)
    ta_name2 <- make.names(ta_name)

    pick_col <- function(nm1, nm2) {
      if (nm1 %in% names(draws)) nm1 else if (nm2 %in% names(draws)) nm2 else NA_character_
    }

    nb_col <- pick_col(nb_name, nb_name2)
    ta_col <- pick_col(ta_name, ta_name2)

    if (!is.na(nb_col) && !is.na(ta_col)) {
      NBc <- draws[[nb_col]]
      TAc <- draws[[ta_col]]

      NB_center_currentinfo <- max(0, mean(NBc), mean(TAc))
      NB_center_perfectinfo <- mean(pmax(NBc, TAc, 0))
      EVPI_center <- max(0, NB_center_perfectinfo - NB_center_currentinfo)

      means_c <- c(treat_none = 0, model = mean(NBc), treat_all = mean(TAc))
      center_strategy <- names(means_c)[which.max(means_c)]

      # sigma for deciding among strategies at that center
      sig0 <- calc_sigma(cbind(NBc, TAc, 0))
      sigma_center_strategy <- if (length(sig0) == 0) NA_real_ else max(sig0, na.rm = TRUE)

      # sigma for EVPI difference at that center
      tmp_c <- pmax(NBc, TAc, 0)
      winnerNB_c <- switch(center_strategy,
                           treat_none = rep(0, length(NBc)),
                           model = NBc,
                           treat_all = TAc)

      sig_ev <- calc_sigma(cbind(tmp_c, winnerNB_c))
      sigma_center_evpi <- if (length(sig_ev) == 0) NA_real_ else max(sig_ev, na.rm = TRUE)

    } else {
      warning("center_row requested, but draws do not include ", nb_name, " and ", ta_name,
              ". Make sure those nodes were monitored in sampling.", call. = FALSE)
    }
  }



  # ---- MC diagnostics ----

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

  # perfect info never changes decision guard
  perfect_info_never_beats_model <- sum(tmp_cluster > winner_NB) == 0
  if (perfect_info_never_beats_model) {
    EVPI_cluster <- 0
    sigma_evpi_cluster <- Inf
  }

  perfect_info_never_beats_model_pop <- has_pop && sum(tmp_pop > ENBnew) == 0
  if (perfect_info_never_beats_model_pop) {
    EVPI_pop <- 0
    sigma_evpi_pop <- Inf
  }


  # flags
  # we do more sampling if the monte carlo error on the difference between two strategies is smaller than 2
  # or if the monte carlo error on the random effects summary measures is more than 5% of the standard deviation of their posterior distribution, i.e. MCSE/posterior sd > 0.05
  flag_sigma <- !is.na(sigma_strategy) && sigma_strategy < sigma_min
  flag_sigma_evpi <- !is.na(sigma_evpi_cluster) && sigma_evpi_cluster < sigma_min
  flag_sigma_evpi_pop <- if (has_pop) {
    !is.na(sigma_evpi_pop) && sigma_evpi_pop < sigma_min
  } else FALSE
  flag_sigma_center_strategy <- !is.na(sigma_center_strategy) && sigma_center_strategy < sigma_min
  flag_sigma_center_evpi <- !is.na(sigma_center_evpi) && sigma_center_evpi < sigma_min
  # ESS-based precision checks for key summary nodes
  # target relative MCSE (MCSE / posterior SD)
  # MCSE /posterior sd = (sd/sqrt(ESS))/sd = 1/sqrt(ESS) < 0.05 -> ESS > 400
  ess_min <- (1 / rmcse_max)^2   # equivalently: require ESS > 400
  ess_pooledNB <- ess_pooledNB_TA <- ess_probuseful <- NA_real_
  flag_rmcse <- FALSE
  if (inherits(x, "mcmc.list")) {
    ess <- coda::effectiveSize(x)
    ess_pooledNB    <- unname(ess["pooledNB"])
    ess_pooledNB_TA <- unname(ess["pooledNB_TA"])
    ess_probuseful  <- unname(ess["probuseful"])
    flag_rmcse <- any(c(ess_pooledNB, ess_pooledNB_TA, ess_probuseful) < ess_min, na.rm = TRUE)
  }

  diagnostics <- tibble::tibble(
    n_draws = nrow(draws),
    winner_strategy = winner_strategy,
    sigma_strategy = sigma_strategy,
    sigma_evpi_cluster = sigma_evpi_cluster,
    sigma_evpi_pop = sigma_evpi_pop,
    sigma_min = sigma_min,
    rmcse_max = rmcse_max,
    ess_min = ess_min,
    ess_pooledNB = ess_pooledNB,
    ess_pooledNB_TA = ess_pooledNB_TA,
    ess_probuseful = ess_probuseful,
    needs_more_sampling = (
      flag_sigma || flag_sigma_evpi || flag_sigma_evpi_pop ||
        flag_sigma_center_strategy || flag_sigma_center_evpi || flag_rmcse
    )
  )
  if (is.null(center_row)) {
    # center flags should not participate if center not requested
    diagnostics <- dplyr::mutate(
      diagnostics,
      needs_more_sampling = (flag_sigma || flag_sigma_evpi || flag_sigma_evpi_pop || flag_rmcse)
    )
  } else {
    diagnostics <- dplyr::bind_cols(
      diagnostics,
      tibble::tibble(
        sigma_center_strategy = sigma_center_strategy,
        sigma_center_evpi = sigma_center_evpi
      )
    )
  }

  # ---- output metrics ----
  metrics <- tibble::tibble(
    NB_currentinfo = NB_current,
    NB_cluster_perfectinfo = NB_cluster_pi,
    EVPI_cluster_perfectinfo = EVPI_cluster,
    NB_population_perfectinfo = NB_pop_pi,
    EVPI_population_perfectinfo = EVPI_pop,
    NB_cluster_partialprev_pi = NB_cluster_partialprev_pi,
    EVPPI_cluster_perfectprevalenceinfo = EVPPI_prev,
    diff_modelvsTA = diff_modelvsTA_label
  )
  if (!is.null(center_row)) {
    metrics <- dplyr::bind_cols(
      metrics,
      tibble::tibble(
        NB_center_currentinfo = NB_center_currentinfo,
        NB_center_perfectinfo = NB_center_perfectinfo,
        EVPI_center_perfectinfo = EVPI_center,
        center_winner_strategy = center_strategy
      )
    )
  }

  list(metrics = metrics, diagnostics = diagnostics)
}



# Sample from the tri-variate NB meta-analysis model with VOI metrics
sample_voi_draws <- function(data,
                             tp, tn, n_event, n_nonevent,
                             center_row = NULL,
                             center_label_cols = c("Publication", "Country"),
                             prior_type = c("weak", "wishart"),
                             t = NULL,
                             prev_ref = 0.5,
                             J = 1000,
                             inits = NULL,
                             n.chains = 2, n.adapt = 1000,
                             burnin = 3000,
                             sample_by = 2000,
                             thin = 1,
                             auto_resample = TRUE,
                             max_rounds = 10,
                             max_draws = 200000,
                             sigma_min = 2,
                             rmcse_max = 0.05,
                             compute_evppi_prev = TRUE,
                             evppi_method = "sal",
                             weak_priors = list(),
                             wishart_priors = list(),
                             check_convergence = c("trace", "none"),
                             diag_vars = c("pooledsens","pooledprev","pooledspec","pooledNB", "pooledNB_TA")) {

  check_convergence <- match.arg(check_convergence)
  prior_type <- match.arg(prior_type)

  # ---- require t ----
  if (is.null(t) || !is.numeric(t) || length(t) != 1 || t <= 0 || t >= 1) {
    stop("Argument `t` must be a single value in (0, 1) and must be supplied.", call. = FALSE)
  }

  # ---- data must be a data frame ----
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)

  # resolve column names (unquoted or strings)
  tp_col         <- resolve_col(rlang::enquo(tp),         data, "tp")
  tn_col         <- resolve_col(rlang::enquo(tn),         data, "tn")
  n_event_col    <- resolve_col(rlang::enquo(n_event),    data, "n_event")
  n_nonevent_col <- resolve_col(rlang::enquo(n_nonevent), data, "n_nonevent")

  tp_vec         <- as.numeric(data[[tp_col]])
  tn_vec         <- as.numeric(data[[tn_col]])
  n_event_vec    <- as.numeric(data[[n_event_col]])
  n_nonevent_vec <- as.numeric(data[[n_nonevent_col]])

  if (anyNA(tp_vec) || anyNA(tn_vec) || anyNA(n_event_vec) || anyNA(n_nonevent_vec)) {
    stop("Missing values detected in tp/tn/n_event/n_nonevent. Please clean data first.", call. = FALSE)
  }

  # JAGS binomial needs integer counts
  tp_vec         <- round(tp_vec)
  tn_vec         <- round(tn_vec)
  n_event_vec    <- round(n_event_vec)
  n_nonevent_vec <- round(n_nonevent_vec)

  if (any(tp_vec < 0) || any(tn_vec < 0) || any(n_event_vec < 0) || any(n_nonevent_vec < 0)) {
    stop("Counts must be non-negative.", call. = FALSE)
  }
  if (any(tp_vec > n_event_vec)) stop("Some tp values exceed n_event.", call. = FALSE)
  if (any(tn_vec > n_nonevent_vec)) stop("Some tn values exceed n_nonevent.", call. = FALSE)

  n_study <- nrow(data)
  n_vec <- n_event_vec + n_nonevent_vec

  if (!is.numeric(J) || length(J) != 1 || J <= 0) stop("`J` must be a positive integer.", call. = FALSE)
  J <- as.integer(J)

  # ---- extended VOI node set ----
  voi_nodes <- c(
    "NBnew", "NBnew_TA",
    "ENBnew", "ENBnew_TA",
    "prevnew",
    "pooledNB", "pooledNB_TA",
    "probuseful"
  )

  # ---- optional: center-specific EVPI for an existing study row ----
  center_meta <- NULL
  center_nodes <- character(0)

  if (!is.null(center_row)) {
    if (!is.numeric(center_row) || length(center_row) != 1 || is.na(center_row)) {
      stop("`center_row` must be a single row number.", call. = FALSE)
    }
    center_row <- as.integer(center_row)
    if (center_row < 1 || center_row > n_study) {
      stop("`center_row` out of range: got ", center_row,
           ", but data has ", n_study, " rows.", call. = FALSE)
    }

    # monitor existing per-study NB nodes
    center_nodes <- c(
      sprintf("NB[%d]", center_row),
      sprintf("NB_TA[%d]", center_row)
    )

    voi_nodes <- unique(c(voi_nodes, center_nodes))

    # store a human-readable label so output is self-documenting
    keep_cols <- intersect(center_label_cols, names(data))
    if (length(keep_cols) == 0) keep_cols <- names(data)[1]  # fallback
    center_meta <- list(
      center_row = center_row,
      center_label = paste(
        sprintf("%s=%s", keep_cols, as.character(data[center_row, keep_cols, drop = TRUE])),
        collapse = ", "
      )
    )
    message("Computing center-specific EVPI for row ", center_row,
            if (!is.null(center_meta)) paste0(" (", center_meta$center_label, ")") else "")

  }


  # ---- base JAGS data ----
  jags_data <- list(
    n_study    = n_study,
    n          = n_vec,
    n_event    = n_event_vec,
    n_nonevent = n_nonevent_vec,
    tp         = tp_vec,
    tn         = tn_vec,
    t          = t,
    prev_ref   = prev_ref,
    J          = J
  )

  # ---- priors ----
  default_weak_priors <- list(
    mu_etap=0, tau_etap=0.001,
    mu_lambdasens0=0, tau_lambdasens0=0.001,
    mu_lambdaspec0=0, tau_lambdaspec0=0.001,
    mu_zss=-0.2, tau_zss=4,
    mu_zsp=-0.2, tau_zsp=4,
    a_csp=-0.99, b_csp=0.99,
    tau_varprev=0.25, tau_varsens=0.25, tau_varspec=0.25
  )
  default_wishart_priors <- list(
    mn   = c(0,0,0),
    prec = diag(1/1000, 3),
    R    = diag(1/10,   3),
    df   = 3
  )

  if (prior_type == "weak") {
    unknown <- setdiff(names(weak_priors), names(default_weak_priors))
    if (length(unknown)) warning("Unknown weak_priors names ignored: ", paste(unknown, collapse=", "))
    priors_used <- utils::modifyList(default_weak_priors, weak_priors)
    jags_data   <- utils::modifyList(jags_data, priors_used)
    model_string <- get_tri_model_weak(include_EVPI = TRUE, J = J)

  } else {
    unknown <- setdiff(names(wishart_priors), names(default_wishart_priors))
    if (length(unknown)) warning("Unknown wishart_priors names ignored: ", paste(unknown, collapse=", "))
    priors_used <- utils::modifyList(default_wishart_priors, wishart_priors)

    if (length(priors_used$mn) != 3) stop("wishart_priors$mn must be length 3", call. = FALSE)
    if (!all(dim(priors_used$prec) == c(3,3))) stop("wishart_priors$prec must be 3x3", call. = FALSE)
    if (!all(dim(priors_used$R) == c(3,3))) stop("wishart_priors$R must be 3x3", call. = FALSE)

    jags_data <- utils::modifyList(jags_data, priors_used)

    # IMPORTANT: get_tri_model_wishart must accept include_EVPI/J
    model_string <- get_tri_model_wishart(include_EVPI = TRUE, J = J)
  }

  priors_used_df <- tibble::tibble(
    prior = names(priors_used),
    value = vapply(priors_used, function(x) {
      if (is.numeric(x) && length(x) == 1) format(x, digits = 8)
      else if (is.numeric(x)) paste(x, collapse = ", ")
      else if (is.matrix(x)) paste0("matrix(", nrow(x), "x", ncol(x), ")")
      else as.character(x)
    }, character(1))
  )

  # ---- compile + burnin ----
  model <- rjags::jags.model(
    textConnection(model_string),
    data = jags_data,
    inits = inits,
    n.chains = n.chains,
    n.adapt = n.adapt
  )
  update(model, burnin)

  # optional traceplot after burn-in
  if (check_convergence == "trace" && length(diag_vars) > 0) {
    diag_samps <- rjags::coda.samples(model, variable.names = diag_vars, n.iter = 1000, thin = 1)
    coda::traceplot(diag_samps)
  }

  # ---- iterative sampling ----
  all_draws <- NULL        # tibble (for returning)
  all_samps <- NULL        # mcmc.list (for ESS + diagnostics)
  last_voi  <- NULL

  round_id <- 0
  repeat {
    round_id <- round_id + 1

    samps <- rjags::coda.samples(
      model,
      variable.names = voi_nodes,
      n.iter = sample_by,
      thin = thin
    )

    # accumulate mcmc.list for ESS-based diagnostics
    all_samps <- append_mcmc_list(all_samps, samps)

    # accumulate tibble for user-facing draws
    new_draws <- mcmc_list_to_draws(samps)
    all_draws <- if (is.null(all_draws)) new_draws else dplyr::bind_rows(all_draws, new_draws)

    # IMPORTANT: pass mcmc.list so ESS can be computed correctly
    last_voi <- compute_voi_metrics(
      all_samps,
      sigma_min = sigma_min,
      rmcse_max = rmcse_max,
      compute_evppi_prev = compute_evppi_prev,
      evppi_method = evppi_method,
      center_row = center_row
    )

    # stop logic
    if (!auto_resample) break
    if (!isTRUE(last_voi$diagnostics$needs_more_sampling[[1]])) break
    if (round_id >= max_rounds) break
    if (nrow(all_draws) >= max_draws) break
  }


  stop_reason <- dplyr::case_when(
    !auto_resample ~ "auto_resample=FALSE",
    !is.null(last_voi) && !isTRUE(last_voi$diagnostics$needs_more_sampling[[1]]) ~ "diagnostics_passed",
    round_id >= max_rounds ~ "max_rounds_reached",
    nrow(all_draws) >= max_draws ~ "max_draws_reached",
    TRUE ~ "stopped"
  )


  diagnostics <- dplyr::mutate(
    last_voi$diagnostics,
    rounds = round_id,
    total_draws = nrow(all_draws),
    stop_reason = stop_reason
  )

  if (auto_resample && isTRUE(diagnostics$needs_more_sampling[[1]]) &&
      diagnostics$stop_reason[[1]] != "diagnostics_passed") {
    warning("VOI sampling stopped before diagnostics passed (", diagnostics$stop_reason[[1]],
            "). Consider increasing max_rounds or max_draws.")
  }

  list(
    draws = all_draws,
    metrics = last_voi$metrics,
    diagnostics = diagnostics,
    priors_used = priors_used_df,
    center = center_meta,
    meta = list(
      prior_type = prior_type,
      t = t,
      prev_ref = prev_ref,
      J = J,
      n_study = n_study,
      burnin = burnin,
      sample_by = sample_by,
      thin = thin,
      n.chains = n.chains,
      center_row = if (is.null(center_row)) NA_integer_ else center_row
    )
  )
}

# wrapper function for VOI sampling from tri-variate NB meta-analysis
MA_NB_tri_voi <- function(data,
                          tp, tn, n_event, n_nonevent,
                          prior_type = c("weak", "wishart"),
                          t,
                          prev_ref = 0.5,
                          center_row = NULL,
                          center_label_cols = c("Publication", "Country"),
                          J = 1000,
                          inits = NULL,
                          n.chains = 2, n.adapt = 1000,
                          burnin = 3000,
                          sample_by = 2000,
                          thin = 1,
                          auto_resample = TRUE,
                          max_rounds = 10,
                          max_draws = 200000,
                          sigma_min = 2,
                          rmcse_max = 0.05,
                          compute_evppi_prev = TRUE,
                          evppi_method = "sal",
                          weak_priors = list(),
                          wishart_priors = list(),
                          check_convergence = c("trace", "none"),
                          diag_vars = c("pooledsens","pooledprev","pooledspec","pooledNB", "pooledNB_TA"),
                          return_draws = FALSE) {

  fit <- sample_voi_draws(
    data = data,
    tp = {{ tp }}, tn = {{ tn }},
    n_event = {{ n_event }}, n_nonevent = {{ n_nonevent }},
    prior_type = prior_type,
    t = t,
    prev_ref = prev_ref,
    center_row = center_row,
    center_label_cols = center_label_cols,
    J = J,
    inits = inits,
    n.chains = n.chains, n.adapt = n.adapt,
    burnin = burnin,
    sample_by = sample_by,
    thin = thin,
    auto_resample = auto_resample,
    max_rounds = max_rounds,
    max_draws = max_draws,
    sigma_min = sigma_min,
    rmcse_max = rmcse_max,
    compute_evppi_prev = compute_evppi_prev,
    evppi_method = evppi_method,
    weak_priors = weak_priors,
    wishart_priors = wishart_priors,
    check_convergence = check_convergence,
    diag_vars = diag_vars
  )

  out <- list(
    metrics = fit$metrics,
    diagnostics = fit$diagnostics,
    priors_used = fit$priors_used,
    meta = fit$meta,
    center = fit$center
  )

  if (isTRUE(return_draws)) {
    out$draws <- fit$draws
  }

  class(out) <- c("MA_NB_tri_voi", class(out))
  out
}
