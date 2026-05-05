#' Bayesian trivariate meta-analysis of net benefit
#'
#' @description
#' Fits a Bayesian trivariate random-effects meta-analysis model jointly
#' modelling sensitivity, specificity, and prevalence across studies using
#' study-level data. Returns posterior MCMC samples for Net Benefit (NB),
#' Relative Utility (RU), and related clinical utility metrics at a
#' user-specified decision threshold. When `compute_EVPI = TRUE`,
#' also computes Value of Information (VOI) metrics including the Expected
#' Value of Perfect Information (EVPI) and optionally the Expected Value of
#' Partial Perfect Information (EVPPI) for prevalence, with optional
#' iterative resampling until Monte Carlo precision criteria are met.
#'
#' @param data A data frame containing one row per study.
#' @param tp Column of true positives (unquoted name or string).
#' @param tn Column of true negatives (unquoted name or string).
#' @param n_event Column of total events per study.
#' @param n_nonevent Column of total non-events per study.
#' @param prior_type Character string specifying the prior structure.
#'   `"weak"` uses weakly informative priors. `"wishart"` uses an inverse-Wishart
#'   prior. Default `"weak"`.
#' @param t Numeric scalar in (0, 1). The decision threshold (i.e. harm-to-benefit ratio).
#'   Required; no default.
#' @param prev_known Numeric scalar in (0, 1). An assumed known prevalence
#'   used to compute additional NB/RU estimates at a fixed prevalence. Only
#'   returned if `return_known = TRUE`. Default `0.5`.
#' @param return_known Logical. If `TRUE`, also returns posteriors
#'   evaluated at `prev_known` (e.g. `pooledNB_known`,
#'   `NBnew_known`). Default `FALSE`.
#' @param compute_EVPI Logical. If `TRUE`, enables VOI computation.
#'   Forces the necessary JAGS nodes into the monitored set and activates
#'   the iterative sampling loop if `auto_resample = TRUE`.
#'   Default `FALSE`.
#' @param J Positive integer. Number of Monte Carlo samples used for the
#'   within-JAGS population EVPI approximation. Only used when
#'   `compute_EVPI = TRUE`. Default `1000`.
#' @param inits Optional list of length `n.chains`, where each element
#'   is a named list of initial values for one JAGS chain. If `NULL`
#'   and `seed` is provided, RNG seeds are set automatically per chain.
#' @param n.chains Positive integer. Number of MCMC chains. Default `2`.
#' @param n.adapt Positive integer. Number of JAGS adaptation iterations.
#'   Default `1000`.
#' @param burnin Positive integer. Number of burn-in iterations (discarded).
#'   Default `3000`.
#' @param iter Positive integer. Number of posterior sampling iterations per
#'   chain in the first sampling round. Default `1000`.
#' @param thin Positive integer. Thinning interval for MCMC chains.
#'   Default `1`.
#' @param seed Optional integer. Random seed for reproducibility, passed to
#'   set each JAGS chain's RNG. Ignored if `inits` already contains
#'   `.RNG.seed`.
#' @param rng_name Character string. JAGS RNG type. One of
#'   `"base::Wichmann-Hill"` (default),
#'   `"base::Marsaglia-Multicarry"`, `"base::Super-Duper"`, or
#'   `"base::Mersenne-Twister"`.
#' @param diag_vars Character vector of JAGS node names to monitor for
#'   convergence diagnostics (trace plots). Set to `NULL` or
#'   `character(0)` to suppress trace plots. Default:
#'   `c("pooledsens", "pooledprev", "pooledspec", "pooledNB",
#'   "pooledNB_TA", "pooledRU")`.
#' @param return_vars Character vector of variable families or individual
#'   JAGS node names to monitor and return. Family names (`"NB"`,
#'   `"RU"`, `"probuseful"`, `"sens"`, `"spec"`) are
#'   automatically expanded to include per-study, pooled, and predictive
#'   nodes. When `compute_EVPI = TRUE`, VOI-required nodes are added
#'   automatically regardless of this argument. Default:
#'   `c("NB", "probuseful")`.
#' @param weak_priors Named list of scalar overrides for the weak prior
#'   hyperparameters. Valid names: `mu_etap`, `tau_etap`,
#'   `mu_lambdasens0`, `tau_lambdasens0`, `mu_lambdaspec0`,
#'   `tau_lambdaspec0`, `mu_zss`, `tau_zss`, `mu_zsp`,
#'   `tau_zsp`, `a_csp`, `b_csp`, `tau_varprev`,
#'   `tau_varsens`, `tau_varspec`. Unknown names are ignored with
#'   a warning. Only used when `prior_type = "weak"`. Default values are
#'   based on Wynants et al. (2018).
#' @param wishart_priors Named list of overrides for the Wishart prior
#'   hyperparameters. Valid names: `mn` (length-3 numeric vector),
#'   `prec` (3x3 matrix), `R` (3x3 matrix), `df` (numeric
#'   scalar). Unknown names are ignored with a warning. Only used when
#'   `prior_type = "wishart"`. Default values are based on Wynants et al. (2018).
#' @param auto_resample Logical. If `TRUE` and `compute_EVPI = TRUE`,
#'   iteratively draws additional samples of size `sample_by` until
#'   Monte Carlo precision criteria (`sigma_min`, `ess_min`) are
#'   met or stopping conditions are reached. Default `TRUE`.
#' @param max_rounds Positive integer. Maximum number of sampling rounds
#'   when `auto_resample = TRUE`. Default `10`.
#' @param sample_by Positive integer. Number of additional iterations per
#'   chain drawn in each resampling round after the first. Default `2000`.
#' @param max_draws Positive integer. Maximum total number of posterior draws
#'   (across all chains) before stopping resampling. Default `200000`.
#' @param sigma_min Numeric. Minimum acceptable signal-to-noise ratio for the
#'   difference between competing strategies; resampling continues while
#'   sigma < `sigma_min`. Default `2`.
#' @param ess_min Numeric. Minimum acceptable effective sample size for key
#'   pooled nodes (`pooledNB`, `pooledNB_TA`, `probuseful`);
#'   resampling continues while any ESS < `ess_min`. Corresponds to a
#'   relative Monte Carlo standard error of at most 5%. Default `400`.
#' @param compute_evppi_prev Logical. If `TRUE` and
#'   `compute_EVPI = TRUE`, computes the EVPPI for prevalence using the
#'   \pkg{voi} package. Requires \pkg{voi} to be installed. Default
#'   `TRUE`.
#' @param evppi_method Character string. Method passed to `voi::evppi()`
#'   for EVPPI estimation. Default `"sal"`.
#' @param center_rows Optional integer vector. Row indices of studies in
#'   `data` for which center-specific EVPI should be computed. Only
#'   active when `compute_EVPI = TRUE`. These rows determine which
#'   per-study NB nodes are monitored and which centers are included in the
#'   resampling stopping rule. Default `NULL`.
#' @param center_label_cols Optional character vector of column names in
#'   `data` used to construct human-readable center labels in the
#'   output. If `NULL` or none of the specified columns exist in
#'   `data`, the first column of `data` is used. Default
#'   `NULL`.
#'
#' @return A named list with the following elements:
#'
#'   - `samples`: `mcmc.list` object containing posterior
#'     draws for all monitored nodes.
#'   ` priors_used`: A tibble summarizing the prior hyperparameter
#'     values used after merging defaults with any user overrides.
#'   - `meta`: A named list of run metadata including threshold
#'     (`returned`), MCMC settings, actual draws per chain
#'     (`saved_per_chain`), total draws (`total_saved`), and
#'     per-chain RNG information. When `compute_EVPI = TRUE`, also
#'     includes `rounds`, `stop_reason`.
#'   - `voi_metrics`: A tibble of VOI metrics including
#'     `EVPI_cluster_perfectinfo`, `EVPI_population_perfectinfo`,
#'     and `EVPPI_cluster_perfectprevalenceinfo`. `NULL` when
#'     `compute_EVPI = FALSE`.
#'   - `voi_diagnostics`: A tibble of Monte Carlo precision
#'     diagnostics including sigma values, ESS, and
#'     `needs_more_sampling`. `NULL` when
#'     `compute_EVPI = FALSE`.
#'   - `voi_center_metrics`: A tibble of center-specific EVPI
#'     estimates and diagnostics for each row in `center_rows`.
#'     `NULL` when `center_rows` is `NULL` or
#'     `compute_EVPI = FALSE`.
#'   - `voi_center_meta`: A data frame of human-readable center
#'     labels corresponding to `center_rows`. `NULL` when
#'     `center_rows` is `NULL` or `compute_EVPI = FALSE`.
#'
#' @references
#' Wynants L, Riley R, Timmerman D, Van Calster B. Random-effects
#' meta-analysis of the clinical utility of tests and prediction models.
#' \emph{Stat Med} 2018;37(12):2034--52.
#' \doi{10.1002/sim.7653}
#'
#' Jackson C, Heath A (2024). \emph{voi: Value of Information Analysis in R}
#' \doi{10.32614/CRAN.package.voi}
#'
#' @seealso [summarize_tri_ma()], [plot_forest()]
#'
#' @examples
#' \dontrun{
#' # Tri-MA only
#' fit <- MA_NB_tri(
#'   data       = data,
#'   tp         = tp,
#'   tn         = tn,
#'   n_event    = n_event,
#'   n_nonevent = n_nonevent,
#'   t          = 0.1,
#'   prior_type = "weak",
#'   seed       = 123
#' )
#' summary(fit$samples)
#'
#' # With VOI and center-specific EVPI
#' fit_voi <- MA_NB_tri(
#'   data          = data,
#'   tp            = tp,
#'   tn            = tn,
#'   n_event       = n_event,
#'   n_nonevent    = n_nonevent,
#'   t             = 0.1,
#'   prior_type    = "weak",
#'   compute_EVPI  = TRUE,
#'   center_rows   = 1:10,
#'   center_label_cols = c("Publication", "Country"),
#'   seed          = 123
#' )
#' fit_voi$voi_metrics
#' fit_voi$voi_center_metrics
#' }
#'
#' @export
MA_NB_tri <- function(data,
                      tp, tn, n_event, n_nonevent,
                      prior_type = c("weak", "wishart"),
                      t = NULL,
                      prev_known = 0.5,
                      return_known = FALSE,
                      compute_EVPI = FALSE,
                      J = 1000,
                      inits = NULL,
                      n.chains = 2, n.adapt = 1000, burnin = 3000, iter = 1000, thin = 1,
                      seed = NULL, rng_name = "base::Wichmann-Hill",
                      diag_vars = c("pooledsens","pooledprev","pooledspec","pooledNB", "pooledNB_TA", "pooledRU"),
                      return_vars = c("NB", "probuseful"),
                      # user may pass PARTIAL overrides here:
                      weak_priors = list(),
                      wishart_priors = list(),
                      # VOI settings (only active when compute_EVPI = TRUE)
                      auto_resample = TRUE,
                      max_rounds = 10,
                      sample_by = 2000,
                      max_draws = 200000,
                      sigma_min = 2,
                      ess_min = 400,
                      compute_evppi_prev = TRUE,
                      evppi_method = "sal",
                      center_rows = NULL,
                      center_label_cols = NULL
){
  # 1. Resolve and validate study-level input columns (tp, tn, n_event, n_nonevent)
  # 2. Expand requested `return_vars` into actual JAGS node names
  # 3. Construct JAGS data list and prior settings (weak or Wishart)
  # 4. Compile JAGS model and run MCMC (adaptation, burn-in, sampling)
  # 5. Return posterior samples and metadata (e.g. call, settings)

  prior_type <- match.arg(prior_type)

  # require t
  if (is.null(t) || !is.numeric(t) || length(t) != 1 || t <= 0 || t >= 1) {
    stop("Argument `t` must be a single value in (0, 1) and must be supplied.")
  }

  # ---- data must be a data frame ----
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  # resolve column names (unquoted or strings)
  tp_col         <- resolve_col(rlang::enquo(tp),        data, "tp")
  tn_col         <- resolve_col(rlang::enquo(tn),        data, "tn")
  n_event_col    <- resolve_col(rlang::enquo(n_event),   data, "n_event")
  n_nonevent_col <- resolve_col(rlang::enquo(n_nonevent),data, "n_nonevent")

  # extract vectors
  tp_vec         <- data[[tp_col]]
  tn_vec         <- data[[tn_col]]
  n_event_vec    <- data[[n_event_col]]
  n_nonevent_vec <- data[[n_nonevent_col]]

  n_study <- nrow(data)

  # coerce & validate
  tp_vec         <- as.numeric(tp_vec)
  tn_vec         <- as.numeric(tn_vec)
  n_event_vec    <- as.numeric(n_event_vec)
  n_nonevent_vec <- as.numeric(n_nonevent_vec)

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
  if (any(tp_vec > n_event_vec)) {
    stop("Some tp values exceed n_event.", call. = FALSE)
  }
  if (any(tn_vec > n_nonevent_vec)) {
    stop("Some tn values exceed n_nonevent.", call. = FALSE)
  }

  n_vec <- n_event_vec + n_nonevent_vec

  # build the JAGS list (same names the model expects)
  line_data <- list(
    n_study = n_study,
    n = n_vec,
    n_event = n_event_vec,
    n_nonevent = n_nonevent_vec,
    tp = tp_vec,
    tn = tn_vec
  )


  # ---- required data checks ----
  if (compute_EVPI) {
    if (!is.numeric(J) || length(J) != 1 || J <= 0)
      stop("When compute_EVPI = TRUE, J must be a positive integer")
    J <- as.integer(J)
  }

  # ---- expansion logic of returned parameters ----
  # Expand high-level variable families (e.g. "NB", "RU") into the corresponding JAGS node names.
  # Atomic names like "pooledNB", "NBnew", "probuseful_known" pass through unchanged.
  expand_target <- function(x, return_known) {

    if (x == "NB") {
      base <- c(
        "NB",              # per-study NB
        "pooledNB",        # pooled NB
        "pooledNB_TA",     # pooled treat-all NB
        "NBnew",           # future-study NB
        "NBnew_TA"         # future-study treat-all NB
      )
      if (return_known) {
        base <- c(base,
                  "pooledNB_known",
                  "pooledNB_TA_known",
                  "NBnew_known",
                  "NBnew_TA_known")
      }
      return(base)
    }

    if (x == "RU") {
      base <- c(
        "RU",        # per-study RU
        "pooledRU",  # pooled RU
        "RUnew"      # RU in new setting
      )
      if (return_known) {
        base <- c(base,
                  "pooledRU_known",
                  "RUnew_known")
      }
      return(base)
    }

    if (x == "probuseful") {
      # probuseful is scalar
      base <- "probuseful"
      if (return_known) {
        base <- c(base, "probuseful_known")
      }
      return(base)
    }

    if (x == "sens") {
      base <- c(
        "sens",       # per-study sens[i]
        "pooledsens", # pooled
        "sensnew"     # predictive (new setting)
      )
      return(base)
    }

    if (x == "spec") {
      base <- c(
        "spec",       # per-study spec[i]
        "pooledspec", # pooled
        "specnew"     # predictive (new setting)
      )
      return(base)
    }

    # everything else is treated as a literal JAGS node name
    return(x)
  }

  expanded <- unlist(lapply(return_vars, expand_target, return_known = return_known))
  return_vars <- unique(expanded)

  # force VOI nodes when compute_EVPI = TRUE
  if (compute_EVPI) {
    voi_nodes <- c(
      "NBnew", "NBnew_TA",
      "ENBnew", "ENBnew_TA",
      "prevnew",
      "pooledNB", "pooledNB_TA",
      "probuseful"
    )
    return_vars <- unique(c(return_vars, voi_nodes))
  }

  if (return_known) {
    message("Including estimates/predictions at assumed know prevalence: prev_known = ",
            prev_known)
  }

  # ---- optional: center-specific EVPI for selected existing study rows ----
  center_meta  <- NULL
  center_nodes <- character(0)

  if (compute_EVPI && !is.null(center_rows)) {
    # validation (copied from sample_voi_draws)
    if (!is.numeric(center_rows) || anyNA(center_rows))
      stop("`center_rows` must be a numeric vector of study row numbers.", call. = FALSE)

    center_rows <- unique(as.integer(center_rows))

    if (any(center_rows < 1 | center_rows > n_study))
      stop("`center_rows` contains values out of range. Valid rows are 1 to ",
           n_study, ".", call. = FALSE)

    # monitor existing per-study NB nodes for the selected center rows
    center_nodes <- unlist(lapply(center_rows, function(cr)
      c(sprintf("NB[%d]", cr), sprintf("NB_TA[%d]", cr))
    ))

    return_vars <- unique(c(return_vars, center_nodes))

    # store a human-readable label so output is self-documenting
    keep_cols <- intersect(center_label_cols, names(data))
    if (length(keep_cols) == 0) keep_cols <- names(data)[1]

    center_meta <- data.frame(
      center_row = center_rows,
      center_label = vapply(center_rows, function(cr) {
        paste(
          sprintf("%s=%s", keep_cols,
                  as.character(data[cr, keep_cols, drop = TRUE])),
          collapse = ", "
        )
      }, character(1)),
      stringsAsFactors = FALSE
    )

    message("Computing center-specific EVPI for rows: ", paste(center_rows, collapse = ", "))
  }

  # ---- full defaults (do NOT edit by user call) from Wynants et al. 2018 ----
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
    mn  = c(0,0,0),
    prec = diag(1/1000, 3),
    R    = diag(1/10,   3),
    df   = 3
  )

  # ---- assemble data passed to JAGS ----
  jags_data <- line_data
  jags_data$t <- t
  jags_data$prev_known <- prev_known
  if (compute_EVPI) {
    jags_data$J <- J
  }

  # Merge user-specified priors with defaults depending on `prior_type`
  if (prior_type == "weak") {
    # merge defaults <- user overrides
    unknown <- setdiff(names(weak_priors), names(default_weak_priors))
    if (length(unknown)) warning("Unknown weak_priors names ignored: ", paste(unknown, collapse=", "))
    priors_used <- utils::modifyList(default_weak_priors, weak_priors)
    jags_data   <- utils::modifyList(jags_data, priors_used)
    model_string  <- get_tri_model_weak(include_EVPI = compute_EVPI, J = J)

    # sanity check: ensure all required weak prior names present
    need <- names(default_weak_priors)
    miss <- setdiff(need, names(jags_data))
    if (length(miss)) stop("Missing weak prior values for: ", paste(miss, collapse=", "))

  } else {
    unknown <- setdiff(names(wishart_priors), names(default_wishart_priors))
    if (length(unknown)) warning("Unknown wishart_priors names ignored: ", paste(unknown, collapse=", "))
    priors_used <- utils::modifyList(default_wishart_priors, wishart_priors)
    # quick shape checks
    if (length(priors_used$mn) != 3) stop("wishart_priors$mn must be length 3")
    if (!all(dim(priors_used$prec) == c(3,3))) stop("wishart_priors$prec must be 3x3")
    if (!all(dim(priors_used$R)   == c(3,3))) stop("wishart_priors$R must be 3x3")

    jags_data  <- utils::modifyList(jags_data, priors_used)
    model_string <- get_tri_model_wishart(include_EVPI = compute_EVPI, J = J)
  }

  # ---- compile, adapt, burn-in ----
  final_inits <- resolve_jags_inits(inits = inits, n.chains = n.chains, seed = seed, rng_name = rng_name)
  rng_meta <- extract_rng_meta(final_inits)
  model <- rjags::jags.model(textConnection(model_string), data=jags_data, inits=final_inits,
                             n.chains=n.chains, n.adapt=n.adapt)
  stats::update(model, burnin)

  # ---- diagnostics ----
  # trace plots for diagnostic variables
  if (!is.null(diag_vars) && length(diag_vars) > 0) {
    diag_samples <- rjags::coda.samples(model, variable.names=diag_vars, n.iter=1000, thin=1)
    coda::traceplot(diag_samples)
  }

  # ---- sampling ----
  if (!compute_EVPI) {
    samples <- rjags::coda.samples(model, variable.names = return_vars,
                                   n.iter = iter, thin = thin)
    all_samps <- samples
    voi_out   <- NULL

  } else {
    # iterative sampling loop
    all_samps <- NULL
    voi_out   <- NULL
    round_id  <- 0

    repeat {
      round_id  <- round_id + 1
      n_iter_this_round <- if (round_id == 1) iter else sample_by

      samps     <- rjags::coda.samples(model, variable.names = return_vars,
                                       n.iter = n_iter_this_round, thin = thin)
      all_samps <- append_mcmc_list(all_samps, samps)

      # compute_voi_metrics() is called here only to check the stopping rule.
      # The outputs are intentionally not returned; the user calls
      # compute_voi_metrics() on fit$samples after MA_NB_tri() completes.
      voi_out <- compute_voi_metrics(
        all_samps,
        sigma_min          = sigma_min,
        ess_min            = ess_min,
        compute_evppi_prev = compute_evppi_prev,
        evppi_method       = evppi_method,
        center_rows        = center_rows
      )

      if (!auto_resample)                                             break
      if (!isTRUE(voi_out$diagnostics$needs_more_sampling[[1]]))     break
      if (round_id >= max_rounds)                                     break
      if (coda::niter(all_samps) * n.chains >= max_draws)            break
    }

    stop_reason <- dplyr::case_when(
      !auto_resample ~ "auto_resample=FALSE",
      !isTRUE(voi_out$diagnostics$needs_more_sampling[[1]]) ~ "diagnostics_passed",
      round_id >= max_rounds ~ "max_rounds_reached",
      coda::niter(all_samps) * n.chains >= max_draws ~ "max_draws_reached",
      TRUE ~ "stopped"
    )

    voi_out$diagnostics <- dplyr::mutate(
      voi_out$diagnostics,
      rounds      = round_id,
      stop_reason = stop_reason
    )

    if (auto_resample &&
        isTRUE(voi_out$diagnostics$needs_more_sampling[[1]]) &&
        voi_out$diagnostics$stop_reason[[1]] != "diagnostics_passed") {
      warning("VOI sampling stopped before diagnostics passed (",
              voi_out$diagnostics$stop_reason[[1]],
              "). Consider increasing max_rounds or max_draws.")
    }
  }

  # ---- return ----
  # format priors for easy printing
  priors_used_df <- data.frame(
    prior = names(priors_used),
    value = vapply(priors_used, function(x) {
      if (is.numeric(x) && length(x) == 1) {
        format(x, digits = 6)
      } else if (is.matrix(x)) {
        paste0(
          "diag=c(", paste(round(diag(x), 6), collapse = ", "), ")",
          if (any(x[lower.tri(x)] != 0))
            paste0(", off-diag=", round(x[lower.tri(x)][1], 6))
        )
      } else if (is.numeric(x)) {
        paste(x, collapse = ", ")
      } else {
        as.character(x)
      }
    }, character(1)),
    stringsAsFactors = FALSE
  )
  priors_used_tibble <- tibble::as_tibble(priors_used_df)

  meta = list(
    t = t,
    prev_known = prev_known,
    prior_type = prior_type,
    returned = return_vars,
    n.chains = n.chains,
    n.adapt = n.adapt,
    burnin = burnin,
    iter = iter,
    thin = thin,
    saved_per_chain = coda::niter(all_samps),
    total_saved = coda::niter(all_samps) * n.chains,
    seed = seed,
    RNG_name_per_chain = rng_meta$rng_name_per_chain,
    RNG_seed_per_chain = rng_meta$rng_seed_per_chain
  )

  # append VOI sampling metadata if applicable
  if (compute_EVPI && !is.null(voi_out)) {
    meta$rounds      <- voi_out$diagnostics$rounds[[1]]
    meta$stop_reason <- voi_out$diagnostics$stop_reason[[1]]
  }

  # return a list with samples and metadata
  # Return raw MCMC samples; downstream functions handle summarization
  return(list(
    samples = all_samps,
    priors_used = priors_used_tibble,
    meta = meta,
    # only present when compute_EVPI = TRUE:
    voi_metrics     = if (compute_EVPI) voi_out$metrics     else NULL,
    voi_diagnostics = if (compute_EVPI) voi_out$diagnostics else NULL,
    voi_center_metrics  = if (compute_EVPI && !is.null(center_rows) ) voi_out$center_metrics else NULL,
    voi_center_meta          = if (compute_EVPI && !is.null(center_rows) ) center_meta else NULL
  ))

}
