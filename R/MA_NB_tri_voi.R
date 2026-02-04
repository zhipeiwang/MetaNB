# compute VOI metrics from MCMC draws
compute_voi_metrics <- function(x,
                                sigma_min = 2,
                                se_max = 0.05,
                                compute_evppi_prev = TRUE,
                                evppi_method = "sal") {

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
  # For prevalence EVPPI, prevnew is needed.
  # For population EVPI, ENBnew* are needed.

  # ---- baseline: current info ----
  if (!all(c("NBnew", "NBnew_TA") %in% names(draws))) {
    stop("Draws must contain NBnew and NBnew_TA to compute VOI.", call. = FALSE)
  }
  NBnew    <- draws$NBnew
  NBnew_TA <- draws$NBnew_TA

  NB_current <- max(0, mean(NBnew), mean(NBnew_TA))

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


  # SE checks
  se_pooledNB <- if ("pooledNB" %in% names(draws)) stats::sd(draws$pooledNB) / sqrt(nrow(draws)) else NA_real_
  se_pooledNB_TA <- if ("pooledNB_TA" %in% names(draws)) stats::sd(draws$pooledNB_TA) / sqrt(nrow(draws)) else NA_real_
  se_probuseful <- if ("probuseful" %in% names(draws)) stats::sd(draws$probuseful) / sqrt(nrow(draws)) else NA_real_

  # flags
  flag_sigma <- !is.na(sigma_strategy) && sigma_strategy < sigma_min
  flag_sigma_evpi <- !is.na(sigma_evpi_cluster) && sigma_evpi_cluster < sigma_min
  flag_sigma_evpi_pop <- if (has_pop) {
    !is.na(sigma_evpi_pop) && sigma_evpi_pop < sigma_min
  } else FALSE
  flag_se <- any(c(se_pooledNB, se_pooledNB_TA, se_probuseful) > se_max, na.rm = TRUE)

  diagnostics <- tibble::tibble(
    n_draws = nrow(draws),
    winner_strategy = winner_strategy,
    sigma_strategy = sigma_strategy,
    sigma_evpi_cluster = sigma_evpi_cluster,
    sigma_evpi_pop = sigma_evpi_pop,
    sigma_min = sigma_min,
    se_max = se_max,
    se_pooledNB = se_pooledNB,
    se_pooledNB_TA = se_pooledNB_TA,
    se_probuseful = se_probuseful,
    needs_more_sampling = (flag_sigma || flag_sigma_evpi || flag_sigma_evpi_pop || flag_se)

  )

  # ---- output metrics ----
  metrics <- tibble::tibble(
    NB_currentinfo = NB_current,
    NB_cluster_perfectinfo = NB_cluster_pi,
    EVPI_cluster_perfectinfo = EVPI_cluster,
    NB_population_perfectinfo = NB_pop_pi,
    EVPI_population_perfectinfo = EVPI_pop,
    NB_cluster_partialprev_pi = NB_cluster_partialprev_pi,
    EVPPI_cluster_perfectprevalenceinfo = EVPPI_prev
  )

  list(metrics = metrics, diagnostics = diagnostics)
}



# Sample from the tri-variate NB meta-analysis model with VOI metrics
sample_voi_draws <- function(data,
                             tp, tn, n_event, n_nonevent,
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
                             se_max = 0.05,
                             compute_evppi_prev = TRUE,
                             evppi_method = "sal",
                             weak_priors = list(),
                             wishart_priors = list(),
                             check_convergence = c("none", "trace"),
                             diag_vars = c("pooledNB","pooledNB_TA","probuseful")) {

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
    diag_samps <- coda.samples(model, variable.names = diag_vars, n.iter = 1000, thin = 1)
    traceplot(diag_samps)
  }

  # ---- iterative sampling ----
  all_draws <- NULL
  last_voi <- NULL

  round_id <- 0
  repeat {
    round_id <- round_id + 1

    samps <- coda.samples(
      model,
      variable.names = voi_nodes,
      n.iter = sample_by,
      thin = thin
    )

    new_draws <- mcmc_list_to_draws(samps)
    all_draws <- if (is.null(all_draws)) new_draws else dplyr::bind_rows(all_draws, new_draws)

    last_voi <- compute_voi_metrics(
      all_draws,
      sigma_min = sigma_min,
      se_max = se_max,
      compute_evppi_prev = compute_evppi_prev,
      evppi_method = evppi_method
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

  if (auto_resample && isTRUE(diagnostics$needs_more_sampling[[1]]) &&
      diagnostics$stop_reason[[1]] != "diagnostics_passed") {
    warning("VOI sampling stopped before diagnostics passed (", diagnostics$stop_reason[[1]],
            "). Consider increasing max_rounds or max_draws.")
  }

  diagnostics <- dplyr::mutate(
    last_voi$diagnostics,
    rounds = round_id,
    total_draws = nrow(all_draws),
    stop_reason = stop_reason
  )

  list(
    draws = all_draws,
    metrics = last_voi$metrics,
    diagnostics = diagnostics,
    priors_used = priors_used_df,
    meta = list(
      prior_type = prior_type,
      t = t,
      prev_ref = prev_ref,
      J = J,
      n_study = n_study,
      burnin = burnin,
      sample_by = sample_by,
      thin = thin,
      n.chains = n.chains
    )
  )
}

# wrapper function for VOI sampling from tri-variate NB meta-analysis
MA_NB_tri_voi <- function(data,
                          tp, tn, n_event, n_nonevent,
                          prior_type = c("weak", "wishart"),
                          t,
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
                          se_max = 0.05,
                          compute_evppi_prev = TRUE,
                          evppi_method = "sal",
                          weak_priors = list(),
                          wishart_priors = list(),
                          check_convergence = c("none", "trace"),
                          diag_vars = c("pooledNB","pooledNB_TA","probuseful"),
                          return_draws = FALSE) {

  fit <- sample_voi_draws(
    data = data,
    tp = {{ tp }}, tn = {{ tn }},
    n_event = {{ n_event }}, n_nonevent = {{ n_nonevent }},
    prior_type = prior_type,
    t = t,
    prev_ref = prev_ref,
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
    se_max = se_max,
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
    meta = fit$meta
  )

  if (isTRUE(return_draws)) {
    out$draws <- fit$draws
  }

  class(out) <- c("MA_NB_tri_voi", class(out))
  out
}
