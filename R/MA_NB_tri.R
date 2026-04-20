#' Bayesian trivariate meta-analysis of net benefit
#'
#' @description
#' Fits a Bayesian trivariate random-effects meta-analysis model for prevalence,
#' sensitivity, and specificity using study-level data. Returns
#' posterior MCMC samples for Net Benefit (NB), Relative Utility (RU),
#' and related clinical utility metrics at a user-specified decision threshold.
#'
#' @param data A data frame containing study-level data.
#' @param tp Column name or vector of true positives.
#' @param tn Column name or vector of true negatives.
#' @param n_event Column name or vector of number of events.
#' @param n_nonevent Column name or vector of number of non-events.
#' @param prior_type Type of prior specification ("weak" or "wishart").
#' @param t Threshold probability (harm-to-benefit ratio).
#' @param prev_known Optional known prevalence value.
#' @param return_known Logical; whether to return results for known prevalence.
#' @param compute_EVPI Logical; whether to compute EVPI.
#' @param J Number of simulations for VOI calculations.
#' @param inits Optional list of initial values for JAGS.
#' @param n.chains Number of MCMC chains.
#' @param n.adapt Number of adaptation iterations.
#' @param burnin Number of burn-in iterations.
#' @param iter Number of sampling iterations.
#' @param thin Thinning interval.
#' @param seed Optional seed for reproducibility.
#' @param rng_name Name of the JAGS random number generator.
#' @param diag_vars Variables to monitor for diagnostics.
#' @param return_vars Variables to return from posterior samples.
#' @param weak_priors List specifying weak prior parameters.
#' @param wishart_priors List specifying Wishart prior parameters.
#'
#' @returns
#' A list containing:
#' \item{samples}{Posterior samples of model parameters and derived quantities.}
#' \item{priors_used}{Details of prior specifications used.}
#' \item{meta}{Metadata including MCMC settings and RNG information.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Minimal example with weak priors
#' fit <- MA_NB_tri(data, tp = tp, tn = tn,
#'                  n_event = n_event, n_nonevent = n_nonevent,
#'                  t = 0.1, prior_type = "weak", seed = 123)

#' }
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
                      wishart_priors = list()
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

  if (!compute_EVPI && any(c("ENBnew","ENBnew_TA") %in% return_vars)) {
    warning("ENBnew / ENBnew_TA requested but compute_EVPI = FALSE; values will be unavailable")
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

  if (return_known) {
    message("Including estimates/predictions at assumed know prevalence: prev_known = ",
            prev_known)
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

  # ---- final sample ----
  samples <- rjags::coda.samples(model, variable.names=return_vars, n.iter=iter, thin=thin)

  # ---- return ----
  # format priors for easy printing
  priors_used_df <- data.frame(
    prior = names(priors_used),
    value = vapply(priors_used, function(x) {
      if (is.numeric(x) && length(x) == 1) format(x, digits = 6)
      else if (is.numeric(x)) paste(x, collapse = ", ")
      else if (is.matrix(x)) paste0("matrix(", nrow(x), "x", ncol(x), ")")
      else as.character(x)
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
    saved_per_chain = iter / thin,
    total_saved = n.chains * (iter / thin),
    seed = seed,
    RNG_name_per_chain = rng_meta$rng_name_per_chain,
    RNG_seed_per_chain = rng_meta$rng_seed_per_chain
  )

  # return a list with samples and metadata
  # Return raw MCMC samples; downstream functions handle summarization
  return(list(
    samples = samples,
    priors_used = priors_used_tibble,
    meta = meta
  ))

}
