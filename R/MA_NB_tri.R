MA_NB_tri <- function(line_data,
                      prior_type = c("weak", "wishart"),
                      t = NULL,
                      prev_ref = 0.5,
                      return_ref = FALSE,
                      compute_EVPI = FALSE,
                      J = 1000,
                      inits = NULL,
                      n.chains = 2, n.adapt = 1000, burnin = 3000, iter = 1000, thin = 1,
                      diag_vars = c("pooledsens","pooledprev","pooledspec","pooledNB"),
                      return_vars = c("NB", "probharmful"),
                      # user may pass PARTIAL overrides here:
                      weak_priors = list(),
                      wishart_priors = list()
){
  prior_type <- match.arg(prior_type)

  # require t
  if (is.null(t) || !is.numeric(t) || length(t) != 1 || t <= 0 || t >= 1) {
    stop("Argument `t` must be a single value in (0, 1) and must be supplied.")
  }

  # ---- required data checks ----
  missing <- setdiff(c("N","n","n_event","n_nonevent","tp","tn"), names(line_data))
  if (length(missing)) stop("line_data is missing: ", paste(missing, collapse=", "))


  if (compute_EVPI) {
    if (!is.numeric(J) || length(J) != 1 || J <= 0)
      stop("When compute_EVPI = TRUE, J must be a positive integer")
    J <- as.integer(J)
  }

  if (!compute_EVPI && any(c("ENBnew","ENBnew_TA") %in% return_vars)) {
    warning("ENBnew / ENBnew_TA requested but compute_EVPI = FALSE; values will be unavailable")
  }

  # ---- expansion logic of returned parameters ----
  # Only the *family names* "NB", "RU", "probharmful" trigger expansion.
  # Atomic names like "pooledNB", "NBnew", "probharmful_ref" pass through unchanged.
  expand_target <- function(x, return_ref) {

    if (x == "NB") {
      base <- c(
        "NB",              # per-study NB
        "pooledNB",        # pooled NB
        "pooledNB_TA",     # pooled treat-all NB
        "NBnew",           # future-study NB
        "NBnew_TA"         # future-study treat-all NB
      )
      if (return_ref) {
        base <- c(base,
                  "pooledNB_ref",
                  "pooledNB_TA_ref",
                  "NBnew_ref",
                  "NBnew_TA_ref")
      }
      return(base)
    }

    if (x == "RU") {
      base <- c(
        "RU",        # per-study RU
        "pooledRU",  # pooled RU
        "RUnew"      # RU in new setting
      )
      if (return_ref) {
        base <- c(base,
                  "pooledRU_ref",
                  "RUnew_ref")
      }
      return(base)
    }

    if (x == "probharmful") {
      # probharmful is scalar, no "new" / "new_ref" in your model
      base <- "probharmful"
      if (return_ref) {
        base <- c(base, "probharmful_ref")
      }
      return(base)
    }

    # everything else is treated as a literal JAGS node name
    return(x)
  }

  expanded <- unlist(lapply(return_vars, expand_target, return_ref = return_ref))
  return_vars <- unique(expanded)

  if (return_ref) {
    message("Including estimates/predictions at assumed reference prevalence: prev_ref = ",
            prev_ref)
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
  jags_data$prev_ref <- prev_ref
  if (compute_EVPI) {
    jags_data$J <- J
  }

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
    model_string <- get_tri_model_wishart()
  }

  # ---- compile, adapt, burn-in ----
  model <- rjags::jags.model(textConnection(model_string), data=jags_data, inits=inits,
                             n.chains=n.chains, n.adapt=n.adapt)
  update(model, burnin)

  # ---- diagnostics ----
  # trace plots for diagnostic variables
  if (!is.null(diag_vars) && length(diag_vars) > 0) {
    diag_samples <- coda.samples(model, variable.names=diag_vars, n.iter=1000, thin=1)
    traceplot(diag_samples)
  }

  # ---- final sample ----
  samples <- coda.samples(model, variable.names=return_vars, n.iter=iter, thin=thin)

  attr(samples,"t") <- t
  attr(samples,"prev_ref") <- prev_ref
  attr(samples,"prior_type") <- prior_type
  attr(samples,"returned") <- return_vars
  attr(samples,"saved_per_chain") <- iter/thin
  attr(samples,"total_saved") <- n.chains * (iter/thin)
  samples
}
