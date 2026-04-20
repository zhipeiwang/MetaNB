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
  flag_sigma_center_strategy <- FALSE
  flag_sigma_center_evpi <- FALSE
  if (!is.null(center_metrics)) {
    flag_sigma_center_strategy <- any(center_metrics$sigma_center_strategy < sigma_min, na.rm = TRUE)
    flag_sigma_center_evpi <- any(center_metrics$sigma_center_evpi < sigma_min, na.rm = TRUE)
  }
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
    flag_ess <- any(c(ess_pooledNB, ess_pooledNB_TA, ess_probuseful) < ess_min, na.rm = TRUE)
  }

  diagnostics <- tibble::tibble(
    n_draws = nrow(draws),
    winner_strategy = winner_strategy,
    sigma_strategy = sigma_strategy,
    sigma_evpi_cluster = sigma_evpi_cluster,
    sigma_evpi_pop = sigma_evpi_pop,
    sigma_min = sigma_min,
    ess_min = ess_min,
    ess_pooledNB = ess_pooledNB,
    ess_pooledNB_TA = ess_pooledNB_TA,
    ess_probuseful = ess_probuseful,
    min_sigma_center_strategy = if (!is.null(center_metrics)) min(center_metrics$sigma_center_strategy, na.rm = TRUE) else NA_real_,
    min_sigma_center_evpi = if (!is.null(center_metrics)) min(center_metrics$sigma_center_evpi, na.rm = TRUE) else NA_real_,
    n_centers_requested = if (is.null(center_rows)) 0L else length(center_rows),
    needs_more_sampling = (
      flag_sigma || flag_sigma_evpi || flag_sigma_evpi_pop ||
        flag_sigma_center_strategy || flag_sigma_center_evpi || flag_ess
    )
  )


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

  list(metrics = metrics, diagnostics = diagnostics, center_metrics = center_metrics)
}



# Sample from the tri-variate NB meta-analysis model with VOI metrics
sample_voi_draws <- function(data,
                             tp, tn, n_event, n_nonevent,
                             center_rows = NULL,
                             center_label_cols = c("Publication", "Country"),
                             prior_type = c("weak", "wishart"),
                             t = NULL,
                             prev_known = 0.5,
                             J = 1000,
                             inits = NULL,
                             seed = NULL,
                             rng_name = "base::Wichmann-Hill",
                             n.chains = 2, n.adapt = 1000,
                             burnin = 3000,
                             sample_by = 2000,
                             thin = 1,
                             auto_resample = TRUE,
                             max_rounds = 10,
                             max_draws = 200000,
                             sigma_min = 2,
                             ess_min = 400,
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

  # ---- optional: center-specific EVPI for selected existing study rows ----
  center_meta <- NULL
  center_nodes <- character(0)

  if (!is.null(center_rows)) {
    if (!is.numeric(center_rows) || anyNA(center_rows)) {
      stop("`center_rows` must be a numeric vector of study row numbers.", call. = FALSE)
    }

    center_rows <- unique(as.integer(center_rows))

    if (any(center_rows < 1 | center_rows > n_study)) {
      stop(
        "`center_rows` contains values out of range. Valid rows are 1 to ", n_study, ".",
        call. = FALSE
      )
    }

    # monitor existing per-study NB nodes
    center_nodes <- unlist(lapply(center_rows, function(cr) {
      c(sprintf("NB[%d]", cr), sprintf("NB_TA[%d]", cr))
    }))

    voi_nodes <- unique(c(voi_nodes, center_nodes))

    # store a human-readable label so output is self-documenting
    keep_cols <- intersect(center_label_cols, names(data))
    if (length(keep_cols) == 0) keep_cols <- names(data)[1]

    center_meta <- data.frame(
      center_row = center_rows,
      center_label = vapply(center_rows, function(cr) {
        paste(
          sprintf("%s=%s", keep_cols, as.character(data[cr, keep_cols, drop = TRUE])),
          collapse = ", "
        )
      }, character(1)),
      stringsAsFactors = FALSE
    )

    message(
      "Computing center-specific EVPI for rows: ",
      paste(center_rows, collapse = ", ")
    )
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
    prev_known   = prev_known,
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

  priors_used_tibble <- tibble::tibble(
    prior = names(priors_used),
    value = vapply(priors_used, function(x) {
      if (is.numeric(x) && length(x) == 1) format(x, digits = 8)
      else if (is.numeric(x)) paste(x, collapse = ", ")
      else if (is.matrix(x)) paste0("matrix(", nrow(x), "x", ncol(x), ")")
      else as.character(x)
    }, character(1))
  )

  # ---- compile + burnin ----
  final_inits <- resolve_jags_inits(
    inits = inits,
    n.chains = n.chains,
    seed = seed,
    rng_name = rng_name
  )
  rng_meta <- extract_rng_meta(final_inits)

  model <- rjags::jags.model(
    textConnection(model_string),
    data = jags_data,
    inits = final_inits,
    n.chains = n.chains,
    n.adapt = n.adapt
  )
  stats::update(model, burnin)

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
      ess_min = ess_min,
      compute_evppi_prev = compute_evppi_prev,
      evppi_method = evppi_method,
      center_rows = center_rows
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
    center_metrics = last_voi$center_metrics,
    priors_used = priors_used_tibble,
    center = center_meta,
    meta = list(
      prior_type = prior_type,
      t = t,
      prev_known = prev_known,
      J = J,
      n_study = n_study,
      burnin = burnin,
      sample_by = sample_by,
      thin = thin,
      n.chains = n.chains,
      seed = seed,
      RNG_name_per_chain = rng_meta$rng_name_per_chain,
      RNG_seed_per_chain = rng_meta$rng_seed_per_chain,
      center_rows = if (is.null(center_rows)) NA_integer_ else center_rows
    )
  )
}

#' Value-of-information analysis for Bayesian trivariate meta-analysis
#'
#' @description
#' Performs value-of-information analysis based on the Bayesian trivariate
#' random-effects meta-analysis model for prevalence, sensitivity, and
#' specificity. The function estimates expected value of perfect information
#' (EVPI), expected value of partial perfect information (EVPPI) for prevalence,
#' and optional center-specific EVPI measures.
#'
#' @param data A data frame containing study-level data.
#' @param tp Column name or vector of true positives.
#' @param tn Column name or vector of true negatives.
#' @param n_event Column name or vector of number of events.
#' @param n_nonevent Column name or vector of number of non-events.
#' @param prior_type Type of prior specification ("weak" or "wishart").
#' @param t Threshold probability (harm-to-benefit ratio).
#' @param prev_known Optional known prevalence value.
#' @param center_rows Optional vector of study rows for center-specific EVPI.
#' @param center_label_cols Optional data columns used to label selected centers.
#' @param J Number of simulations for population-level VOI calculations.
#' @param inits Optional list of initial values for JAGS.
#' @param seed Optional seed for reproducibility.
#' @param rng_name Name of the JAGS random number generator.
#' @param n.chains Number of MCMC chains.
#' @param n.adapt Number of adaptation iterations.
#' @param burnin Number of burn-in iterations.
#' @param sample_by Number of iterations sampled per resampling round.
#' @param thin Thinning interval.
#' @param auto_resample Logical; whether to continue sampling until precision
#' diagnostics are met or a stopping limit is reached.
#' @param max_rounds Maximum number of resampling rounds.
#' @param max_draws Maximum total number of sampled draws.
#' @param sigma_min Minimum sigma threshold used in checks.
#' @param ess_min Minimum effective sample size used in checks.
#' @param compute_evppi_prev Logical; whether to compute EVPPI for prevalence.
#' @param evppi_method Method passed to [voi::evppi()].
#' @param weak_priors List specifying weak prior parameters.
#' @param wishart_priors List specifying Wishart prior parameters.
#' @param check_convergence Whether to show trace plots ("trace" or "none").
#' @param diag_vars Variables to monitor for convergence diagnostics.
#' @param return_draws Logical; whether to return posterior draws used for VOI.
#'
#' @returns
#' An object of class `"MA_NB_tri_voi"` containing:
#' \item{metrics}{Estimated VOI metrics.}
#' \item{diagnostics}{Precision and sampling diagnostics.}
#' \item{center_metrics}{Optional center-specific EVPI summaries.}
#' \item{priors_used}{Summary of prior specifications used.}
#' \item{meta}{Metadata including MCMC settings and RNG information.}
#' \item{center}{Optional labels for selected centers.}
#' \item{draws}{Posterior draws used for VOI calculations, if requested.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Minimal example with weak priors
#' voi_res <- MA_NB_tri_voi(data, tp = tp, tn = tn,
#'                  n_event = n_event, n_nonevent = n_nonevent,
#'                  t = 0.1, prior_type = "weak", seed = 123)

#' }
MA_NB_tri_voi <- function(data,
                          tp, tn, n_event, n_nonevent,
                          prior_type = c("weak", "wishart"),
                          t,
                          prev_known = 0.5,
                          center_rows = NULL,
                          center_label_cols = c("Publication", "Country"),
                          J = 1000,
                          inits = NULL,
                          seed = NULL,
                          rng_name = "base::Wichmann-Hill",
                          n.chains = 2, n.adapt = 1000,
                          burnin = 3000,
                          sample_by = 2000,
                          thin = 1,
                          auto_resample = TRUE,
                          max_rounds = 10,
                          max_draws = 200000,
                          sigma_min = 2,
                          ess_min = 400,
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
    prev_known = prev_known,
    center_rows = center_rows,
    center_label_cols = center_label_cols,
    J = J,
    inits = inits,
    seed = seed,
    rng_name = rng_name,
    n.chains = n.chains, n.adapt = n.adapt,
    burnin = burnin,
    sample_by = sample_by,
    thin = thin,
    auto_resample = auto_resample,
    max_rounds = max_rounds,
    max_draws = max_draws,
    sigma_min = sigma_min,
    ess_min = ess_min,
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
    center_metrics = fit$center_metrics,
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
