# for the MA_NB_tri() function to accept df instead of a list of lists
# turn an argument into a column name
resolve_col <- function(arg, data, arg_name) {
  # arg is a quosure: enquo(tp)
  if (rlang::quo_is_missing(arg)) {
    stop("You must supply `", arg_name, "` when `data` is a data frame.", call. = FALSE)
  }

  expr <- rlang::get_expr(arg)

  # Case 1: user passed a string, e.g. tp = "TP"
  if (is.character(expr)) {
    if (length(expr) != 1) {
      stop("`", arg_name, "` must be a single column name.", call. = FALSE)
    }
    col <- expr

    # Case 2: user passed an unquoted name, e.g. tp = TP
  } else if (is.symbol(expr)) {
    col <- rlang::as_name(expr)

  } else {
    stop("`", arg_name, "` must be a column name (unquoted) or a single string.",
         call. = FALSE)
  }

  if (!col %in% names(data)) {
    stop(
      "Column `", col, "` (from `", arg_name, "`) not found in `data`.\n",
      "Available columns: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }

  col
}

# recursively drop NULLs from a list
# Remove empty components so the returned object stays compact and easier to use.
# Helper for .summarize_for_forestplot() and summarize_tri_ma()
drop_nulls_recursive <- function(x) {
  # Keep data frames intact
  if (is.data.frame(x)) return(x)

  # Keep atomic vectors (including named numeric vectors) intact
  if (!is.list(x)) return(x)

  # Recurse only through true lists
  x <- lapply(x, drop_nulls_recursive)

  # Drop NULL elements from this list level
  x <- x[!vapply(x, is.null, logical(1))]

  # Optional: drop empty lists (but NOT data frames; handled above)
  x <- x[!vapply(x, function(z) is.list(z) && length(z) == 0, logical(1))]

  x
}

# extract mcmc.list from fit object or pass through if already mcmc.list
get_samples <- function(x) {
  if (inherits(x, "mcmc.list")) return(x)
  if (is.list(x) && !is.null(x$samples) && inherits(x$samples, "mcmc.list")) return(x$samples)
  stop("Expected an mcmc.list or a fit object with $samples (mcmc.list).", call. = FALSE)
}


####################-------------------summarizer for forest plot function----------------------######################
.summarize_for_forestplot <- function(
    samples,
    data = NULL,
    label_cols = NULL,
    targets = c("NB", "probuseful"),
    targets_per_study = c("NB"),
    return_known = FALSE
) {
  # 1. Extract posterior summaries from the fitted samples
  # 2. Map user-facing metric names (e.g. NB, RU, sens) to JAGS node names
  # 3. Build per-study summaries when study-level nodes are available
  # 4. Build pooled / predictive / known-prevalence summaries
  # 5. Return a structured list used by user-facing summaries and forest plots

  study_info <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)

    # Keep only requested label columns that actually exist in `data`
    missing_labels <- setdiff(label_cols, names(data))
    if (length(missing_labels)) {
      warning("Missing label_cols ignored: ", paste(missing_labels, collapse = ", "))
    }
    label_cols <- intersect(label_cols, names(data))
    if (is.null(label_cols)) label_cols <- names(data)[1]
    study_info <- data[, label_cols, drop = FALSE]
  }

  # Accept either a fit object or a raw mcmc.list, then compute standard posterior summaries
  samples <- get_samples(samples)
  ss <- summary(samples)
  stats  <- ss$statistics
  quants <- ss$quantiles

  # ---- helpers ----
  # Posterior summaries are built from node names in the MCMC output.
  # `get_summary()` returns mean/median/95% interval,
  # and `get_mean_only()` is used for scalar probabilities such as probuseful.
  get_summary <- function(name) {
    if (!name %in% rownames(stats)) return(NULL)
    c(
      Mean = stats[name, "Mean"],
      Median = quants[name, "50%"],
      Low  = quants[name, "2.5%"],
      High = quants[name, "97.5%"]
    )
  }

  get_mean_only <- function(name) {
    if (!name %in% rownames(stats)) return(NULL)
    c(Mean = stats[name, "Mean"])
  }

  # ---- family specification ----
  # Central mapping from user-facing metric families to JAGS node names.
  family_spec <- list(
    NB = list(
      per_study_base = "NB",

      pooled = list(
        model     = "pooledNB",
        treat_all = "pooledNB_TA"
      ),

      predictive = list(
        model     = "NBnew",
        treat_all = "NBnew_TA"
      ),

      pooled_known = list(
        model     = "pooledNB_known",
        treat_all = "pooledNB_TA_known"
      ),

      pred_known = list(
        model     = "NBnew_known",
        treat_all = "NBnew_TA_known"
      )
    ),

    RU = list(
      per_study_base = "RU",
      pooled = list(model = "pooledRU"),
      predictive = list(model = "RUnew"),
      pooled_known = list(model = "pooledRU_known"),
      pred_known = list(model = "RUnew_known")
    ),

    probuseful = list(
      per_study_base = NULL,
      pooled = list(model = "probuseful"),
      predictive = NULL,
      pooled_known = list(model = "probuseful_known"),
      pred_known = NULL
    ),
    sens = list(
      per_study_base = "sens",
      pooled = list(model = "pooledsens"),
      predictive = list(model = "sensnew"),
      pooled_known = NULL,
      pred_known = NULL
    ),

    spec = list(
      per_study_base = "spec",
      pooled = list(model = "pooledspec"),
      predictive = list(model = "specnew"),
      pooled_known = NULL,
      pred_known = NULL
    )

  )

  out <- list()

  for (tar in targets) {

    # -------------------------------
    # CASE 1: family target (NB / RU / probuseful / sens / spec)
    # -------------------------------
    if (tar %in% names(family_spec)) {

      metric_def <- family_spec[[tar]]

      # ---- per-study ----
      per_study <- NULL
      if (!is.null(metric_def$per_study_base) && tar %in% targets_per_study) {

        # Per-study nodes are stored in the MCMC output as e.g. NB[1], NB[2], ...
        pat  <- paste0("^", metric_def$per_study_base, "\\[\\d+\\]$")
        rows <- grep(pat, rownames(stats), value = TRUE)

        if (length(rows) > 0) {
          if (is.null(study_info)) {
            stop("data required for per-study summaries of ", tar)
          }

          # Parse study indices from the node names and check that they align
          # exactly with the rows of `data`. This prevents accidental relabeling.
          idx <- as.integer(sub(paste0("^", metric_def$per_study_base, "\\[(\\d+)\\]$"),
                                "\\1", rows))
          ord <- order(idx)

          n_data <- nrow(study_info)

          # Basic sanity
          if (anyNA(idx)) {
            stop("Could not parse study indices from MCMC row names for ", tar, ".", call. = FALSE)
          }
          if (any(idx < 1) || any(idx > n_data)) {
            stop(
              "Per-study indices for ", tar, " are out of range.\n",
              "Found indices: ", paste(sort(unique(idx)), collapse = ", "), "\n",
              "But data has ", n_data, " rows.",
              call. = FALSE
            )
          }

          # Strong alignment check to prevent mislabeling
          if (length(idx) != n_data || !setequal(idx, seq_len(n_data))) {
            stop(
              "Mismatch between per-study draws and `data` rows for ", tar, ".\n",
              "MCMC has indices: ", paste(sort(unique(idx)), collapse = ", "), "\n",
              "But data expects: ", paste(seq_len(n_data), collapse = ", "), "\n",
              "This would misalign study labels with estimates, so I am stopping.",
              call. = FALSE
            )
          }

          # add observed rates for sens / spec in study_info
          if (tar %in% c("sens", "spec")) {

            # compute observed rates if the needed columns exist
            has_needed <- all(c("tp", "tn", "n_event", "n_nonevent") %in% names(data))

            if (has_needed) {
              if (tar == "sens") {
                study_info$Observed <- with(data, tp / n_event)
              } else {
                study_info$Observed <- with(data, tn / n_nonevent)
              }
            }
          }

          per_study <- study_info %>%
            dplyr::mutate(.idx = dplyr::row_number()) %>%
            dplyr::arrange(.idx) %>%
            dplyr::mutate(
              Mean = stats[rows, "Mean"][ord],
              Median = quants[rows, "50%"][ord],
              Low  = quants[rows, "2.5%"][ord],
              High = quants[rows, "97.5%"][ord]
            ) %>%
            dplyr::select(-.idx)
        }
      }

      # ---- pooled (CrI) ----
      pooled <- list()
      if (!is.null(metric_def$pooled)) {
        for (k in names(metric_def$pooled)) {
          pooled[[k]] <- get_summary(metric_def$pooled[[k]])
        }
        pooled <- pooled[!sapply(pooled, is.null)]
        if (length(pooled) == 0) pooled <- NULL
      }

      # ---- predictive (PI) ----
      predictive <- list()
      if (!is.null(metric_def$predictive)) {
        for (k in names(metric_def$predictive)) {
          predictive[[k]] <- get_summary(metric_def$predictive[[k]])
        }
        predictive <- predictive[!sapply(predictive, is.null)]
        if (length(predictive) == 0) predictive <- NULL
      }

      # ---- pooled_known (CrI) ----
      pooled_known <- NULL
      if (return_known && !is.null(metric_def$pooled_known)) {
        pooled_known <- list()
        for (k in names(metric_def$pooled_known)) {
          pooled_known[[k]] <- get_summary(metric_def$pooled_known[[k]])
        }
        pooled_known <- pooled_known[!sapply(pooled_known, is.null)]
        if (length(pooled_known) == 0) pooled_known <- NULL
      }

      # ---- pred_known ----
      pred_known <- NULL
      if (return_known && !is.null(metric_def$pred_known)) {
        pred_known <- list()
        for (k in names(metric_def$pred_known)) {
          pred_known[[k]] <- get_summary(metric_def$pred_known[[k]])
        }
        pred_known <- pred_known[!sapply(pred_known, is.null)]
        if (length(pred_known) == 0) pred_known <- NULL
      }

      # ---- probuseful special case (mean only) ----
      if (tar == "probuseful") {
        pooled <- list(model = get_mean_only("probuseful"))
        if (return_known) {
          pooled_known <- list(model = get_mean_only("probuseful_known"))
        }
        predictive <- NULL
        pred_known   <- NULL
      }

      out[[tar]] <- list(
        per_study  = per_study,
        pooled     = pooled,
        predictive = predictive,
        pooled_known = pooled_known,
        pred_known   = pred_known
      )

    } else {

      # -------------------------------
      # CASE 2: atomic variable name
      # -------------------------------
      # Atomic targets are interpreted as direct JAGS node names rather than
      # metric families. They are returned as scalar summaries when available.
      scalar <- NULL
      if (tar %in% rownames(stats)) {
        if (tar %in% rownames(quants)) {
          scalar <- get_summary(tar)
        } else {
          scalar <- get_mean_only(tar)
        }
      }

      out[[tar]] <- list(scalar = scalar)
    }
  }

  # ---- drop NULLs recursively ----
  # Remove empty components so the returned object stays compact and easier to use.
  drop_nulls_recursive(out)


}


####################------------------- user-facing summarizer ----------------------######################

#' Summarize posterior samples from a Bayesian trivariate meta-analysis
#'
#' @description
#' Extracts and structures posterior summaries from a fitted Bayesian
#' trivariate meta-analysis. For each requested metric, returns pooled
#' estimates (posterior mean, median, and 95% credible interval),
#' predictive summaries for a new center (mean, median, and 95% prediction
#' interval), and optionally per-study summaries and summaries at a fixed
#' known prevalence.
#'
#' @param samples A fitted object returned by [MA_NB_tri()], or a
#'   \code{coda::mcmc.list} of posterior samples directly.
#' @param data Optional data frame containing study-level information. Required
#'   when per-study summaries are requested (i.e. when \code{include_per_study
#'   = TRUE} and a metric in \code{per_study_metrics} is requested).
#' @param label_cols Optional character vector of column names in \code{data}
#'   to retain in per-study output (e.g. study labels, country, sample size).
#' @param metrics Character vector of metrics to summarize. Family names
#'   (\code{"NB"}, \code{"RU"}, \code{"probuseful"}) and atomic JAGS node
#'   names (e.g. \code{"pooledprev"}, \code{"NBnew_known"}) are both accepted.
#'   Default: \code{c("NB", "probuseful")}.
#' @param per_study_metrics Character vector of metrics for which per-study
#'   summaries should be returned. Only used when \code{include_per_study =
#'   TRUE}. Default: \code{c("NB")}.
#' @param return_known Logical. If \code{TRUE}, also returns summaries
#'   evaluated at the fixed known prevalence (e.g. \code{pooledNB_known},
#'   \code{NBnew_known}). Requires that \code{return_known = TRUE} was set
#'   when fitting with [MA_NB_tri()]. Default \code{FALSE}.
#' @param include_per_study Logical. If \code{TRUE}, per-study posterior
#'   summaries are included for metrics listed in \code{per_study_metrics}.
#'   Requires \code{data} to be supplied. Default \code{TRUE}.
#'
#' @return A named list with one element per requested metric. For family
#'   metrics (\code{"NB"}, \code{"RU"}), each element contains:
#' \describe{
#'   \item{\code{per_study}}{A data frame of per-study posterior summaries
#'     (mean, median, 95\% CrI), merged with \code{label_cols} from
#'     \code{data}. \code{NULL} if \code{include_per_study = FALSE}.}
#'   \item{\code{pooled}}{Named list of pooled posterior summaries (mean,
#'     median, 95\% CrI) for the model and treat-all strategies.}
#'   \item{\code{predictive}}{Named list of posterior predictive summaries
#'     (mean, median, 95\% PI) for a new center for the model and treat-all
#'     strategies.}
#'   \item{\code{pooled_known}}{As \code{pooled} but at fixed known
#'     prevalence. \code{NULL} if \code{return_known = FALSE}.}
#'   \item{\code{pred_known}}{As \code{predictive} but at fixed known
#'     prevalence. \code{NULL} if \code{return_known = FALSE}.}
#' }
#' For \code{"probuseful"}, a list with \code{value} and
#' optionally \code{value_known}. For atomic node names, a list with
#' \code{scalar} containing the posterior summary.
#'
#' @seealso [MA_NB_tri()], [plot_forest()]
#'
#' @examples
#' \dontrun{
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
#'
#' # Default summaries
#' summary <- summarize_tri_ma(fit, data = data,
#'                         label_cols = c("Publication", "Country"))
#' summary$NB$pooled$model       # pooled NB: mean, median, 95% CrI
#' summary$NB$predictive$model   # predictive NB for new center
#' summary$probuseful$value      # P(useful)
#'
#' # Include known-prevalence summaries
#' fit2 <- MA_NB_tri(..., prev_known = 0.5, return_known = TRUE)
#' summary2 <- summarize_tri_ma(fit2, data = data, return_known = TRUE)
#' summary2$NB$pooled_known$model
#' }
#'
#' @export
summarize_tri_ma <- function(
    samples,
    data = NULL,
    label_cols = NULL,
    metrics = c("NB", "probuseful"),
    per_study_metrics = c("NB"),
    return_known = FALSE,
    include_per_study = TRUE
) {
  core <- .summarize_for_forestplot(
    samples = get_samples(samples),
    data = data,
    label_cols = label_cols,
    targets = metrics,
    targets_per_study = if (include_per_study) per_study_metrics else character(0),
    return_known = return_known
  )

  out <- list()

  # --- NB / RU: keep structure; named vectors (no 1-row dfs) ---
  for (m in intersect(c("NB","RU"), names(core))) {
    obj <- core[[m]]

    out[[m]] <- drop_nulls_recursive(list(
      per_study  = if (include_per_study) obj$per_study else NULL,
      pooled     = obj$pooled,        # named vectors inside
      predictive = obj$predictive,    # named vectors inside
      pooled_known = if (return_known) obj$pooled_known else NULL,
      pred_known   = if (return_known) obj$pred_known else NULL
    ))
  }

  # --- probuseful: scalar value, not pooled/model ---
  if ("probuseful" %in% names(core)) {
    pu <- core$probuseful

    value <- NULL
    if (!is.null(pu$pooled$model) && "Mean" %in% names(pu$pooled$model)) {
      value <- unname(pu$pooled$model["Mean"])
    }

    value_known <- NULL
    if (return_known && !is.null(pu$pooled_known$model) && "Mean" %in% names(pu$pooled_known$model)) {
      value_known <- unname(pu$pooled_known$model["Mean"])
    }

    out$probuseful <- drop_nulls_recursive(list(
      value = value,
      value_known = value_known
    ))
  }

  # atomic scalars (if you allow them in metrics)
  other <- setdiff(metrics, c("NB","RU","probuseful"))
  for (nm in intersect(other, names(core))) {
    out[[nm]] <- core[[nm]]
  }

  drop_nulls_recursive(out)
}

####################------------------- helpers for forest plotter ----------------------######################
# fallback CI for sensitivity and specificity based on Wilson’s method via mada::madad()
compute_sens_spec_ci_mada <- function(tp, tn, n_event, n_nonevent) {
  out <- mada::madad(
    TP = tp,
    FN = n_event - tp,
    FP = n_nonevent - tn,
    TN = tn
  )
  list(
    sens_low = out$sens$sens.ci[, 1],
    sens_high = out$sens$sens.ci[, 2],
    spec_low = out$spec$spec.ci[, 1],
    spec_high = out$spec$spec.ci[, 2]
  )
}

# analytical solution for confidence interval of NB
compute_nb_ci_delta <- function(tp, tn, n_event, n_nonevent, t,
                                conf.level = 0.95) {
  n <- n_event + n_nonevent

  prev_hat <- n_event / n
  sens_hat <- tp / n_event
  spec_hat <- tn / n_nonevent

  nb_hat <- prev_hat * sens_hat -
    (1 - prev_hat) * (1 - spec_hat) * t / (1 - t)

  var_hat <- (1 / n) * (
    prev_hat * sens_hat * (1 - sens_hat) +
      (t^2 / (1 - t)^2) * (1 - prev_hat) * spec_hat * (1 - spec_hat) +
      prev_hat * (1 - prev_hat) *
      (sens_hat + t * (1 - spec_hat) / (1 - t))^2
  )

  se_hat <- sqrt(var_hat)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  tibble::tibble(
    est = nb_hat,
    se = se_hat,
    var = var_hat,
    low = nb_hat - z * se_hat,
    high = nb_hat + z * se_hat
  )
}

# a helper that decides which point estimate and interval to display for each study in the forest plot, based on user-specified columns and fallback rules
build_per_study_display <- function(
    metric,
    per,
    data,
    center = c("Median", "Mean"),
    t = NULL,

    reported_est_col = NULL,
    reported_low_col = NULL,
    reported_high_col = NULL,

    interval_fallback = c("none", "frequentist", "analytic", "model")
) {
  center <- match.arg(center)
  interval_fallback <- match.arg(interval_fallback)

  n <- nrow(per)

  display_est  <- per[[center]]
  display_low  <- per$Low
  display_high <- per$High

  point_source <- rep("model", n)
  interval_source <- rep("model", n)

  ## ---- 1. point estimates ----

  # reported point estimate has priority
  if (!is.null(reported_est_col) && reported_est_col %in% names(data)) {
    rep_est <- suppressWarnings(as.numeric(data[[reported_est_col]]))
    has_rep_est <- !is.na(rep_est)

    # for sens/spec, convert percentages to proportions if needed
    if (metric %in% c("sens", "spec") && any(rep_est > 1, na.rm = TRUE)) {
      rep_est <- rep_est / 100
    }

    display_est[has_rep_est] <- rep_est[has_rep_est]
    point_source[has_rep_est] <- "reported"
  }

  # observed fallback for sens/spec
  if (metric %in% c("sens", "spec") &&
      all(c("tp", "tn", "n_event", "n_nonevent") %in% names(data))) {
    obs_est <- if (metric == "sens") data$tp / data$n_event else data$tn / data$n_nonevent
    use_obs <- point_source != "reported" & !is.na(obs_est)

    display_est[use_obs] <- obs_est[use_obs]
    point_source[use_obs] <- "observed"
  }

  # observed fallback for NB
  if (metric == "NB" &&
      all(c("tp", "tn", "n_event", "n_nonevent") %in% names(data))) {
    if (is.null(t)) {
      stop("`t` must be supplied when building per-study NB from observed study data.",
           call. = FALSE)
    }

    nb_obs <- compute_nb_ci_delta(
      tp = data$tp,
      tn = data$tn,
      n_event = data$n_event,
      n_nonevent = data$n_nonevent,
      t = t
    )

    use_obs <- point_source != "reported" & !is.na(nb_obs$est)
    display_est[use_obs] <- nb_obs$est[use_obs]
    point_source[use_obs] <- "observed"
  }

  # RU stays model-based unless user reported values exist
  # no extra block needed

  ## ---- 2. reported intervals ----

  has_rep_ci <- rep(FALSE, n)

  if (!is.null(reported_low_col) && !is.null(reported_high_col) &&
      reported_low_col %in% names(data) && reported_high_col %in% names(data)) {

    rep_low  <- suppressWarnings(as.numeric(data[[reported_low_col]]))
    rep_high <- suppressWarnings(as.numeric(data[[reported_high_col]]))

    # for sens/spec, convert percentages to proportions if needed
    if (metric %in% c("sens", "spec")) {
      if (any(rep_low > 1, na.rm = TRUE))  rep_low  <- rep_low / 100
      if (any(rep_high > 1, na.rm = TRUE)) rep_high <- rep_high / 100
    }

    has_rep_ci <- !is.na(rep_low) & !is.na(rep_high)

    display_low[has_rep_ci] <- rep_low[has_rep_ci]
    display_high[has_rep_ci] <- rep_high[has_rep_ci]
    interval_source[has_rep_ci] <- "reported"
  }

  ## ---- 3. fallback intervals for rows with missing reported CI ----

  need_ci <- !has_rep_ci

  if (any(need_ci) && interval_fallback != "none") {

    # sens/spec frequentist fallback via mada
    if (metric %in% c("sens", "spec") && interval_fallback == "frequentist") {
      if (all(c("tp", "tn", "n_event", "n_nonevent") %in% names(data))) {
        ci <- compute_sens_spec_ci_mada(
          tp = data$tp,
          tn = data$tn,
          n_event = data$n_event,
          n_nonevent = data$n_nonevent
        )

        if (metric == "sens") {
          display_low[need_ci]  <- ci$sens_low[need_ci]
          display_high[need_ci] <- ci$sens_high[need_ci]
        } else {
          display_low[need_ci]  <- ci$spec_low[need_ci]
          display_high[need_ci] <- ci$spec_high[need_ci]
        }

        interval_source[need_ci] <- "frequentist"
      }
    }

    # NB analytical fallback
    if (metric == "NB" && interval_fallback == "analytic") {
      if (is.null(t)) {
        stop("`t` must be supplied when using analytical CI fallback for NB.",
             call. = FALSE)
      }

      if (all(c("tp", "tn", "n_event", "n_nonevent") %in% names(data))) {
        nb_ci <- compute_nb_ci_delta(
          tp = data$tp,
          tn = data$tn,
          n_event = data$n_event,
          n_nonevent = data$n_nonevent,
          t = t
        )

        display_low[need_ci]  <- nb_ci$low[need_ci]
        display_high[need_ci] <- nb_ci$high[need_ci]
        interval_source[need_ci] <- "analytic"
      }
    }

    # model-based fallback
    if (interval_fallback == "model") {
      display_low[need_ci]  <- per$Low[need_ci]
      display_high[need_ci] <- per$High[need_ci]
      interval_source[need_ci] <- "model"
    }
  }

  tibble::tibble(
    display_est = display_est,
    display_low = display_low,
    display_high = display_high,
    point_source = point_source,
    interval_source = interval_source
  )
}



# choose interval column label for the per-study display column
infer_interval_label <- function(metric, interval_source) {
  # RU is always CrI in current methodology
  if (identical(metric, "RU")) {
    return("95% CrI")
  }

  # keep only non-missing unique values
  src <- unique(stats::na.omit(interval_source))

  # all displayed study-row intervals are model-based
  if (length(src) > 0 && all(src == "model")) {
    return("95% CrI")
  }

  # otherwise default to CI
  "95% CI"
}


# print short messages explaining fallback sources actually used
emit_forest_messages <- function(metric, point_source, interval_source,
                                 interval_fallback, mark_imputed = FALSE) {

  # keep non-missing row-level sources
  point_src <- point_source[!is.na(point_source)]
  int_src   <- interval_source[!is.na(interval_source)]

  # non-reported sources actually used
  point_nonrep <- unique(point_src[point_src != "reported"])
  int_nonrep   <- unique(int_src[int_src != "reported"])

  # counts
  n_point_total  <- length(point_src)
  n_point_nonrep <- sum(point_src != "reported")

  n_int_total  <- length(int_src)
  n_int_nonrep <- sum(int_src != "reported")

  # choose lead-in
  point_leadin <- if (n_point_nonrep == n_point_total) "All" else "Some"
  int_leadin   <- if (n_int_nonrep == n_int_total) "All" else "Some"

  # ---- point estimate message ----
  if (length(point_nonrep) > 0) {

    point_text <- NULL

    if (setequal(point_nonrep, "observed")) {
      point_text <- "study-level observed or calculated values"
    } else if (setequal(point_nonrep, "model")) {
      point_text <- "Bayesian trivariate model-based values"
    } else if (setequal(point_nonrep, c("observed", "model"))) {
      point_text <- paste(
        "study-level observed or calculated values,",
        "or Bayesian trivariate model-based values"
      )
    }

    if (!is.null(point_text)) {
      if (isTRUE(mark_imputed)) {
        message("\u2020 indicates point estimates replaced by ", point_text,
                " when reported estimates were unavailable.")
      } else {
        message(point_leadin, " point estimates were replaced by ", point_text,
                " because reported estimates were unavailable.")
      }
    }
  }

  # ---- interval message ----
  if (length(int_nonrep) > 0) {

    interval_text <- NULL

    if (setequal(int_nonrep, "frequentist")) {
      interval_text <- "confidence interval fallback based on Wilson's method"
    } else if (setequal(int_nonrep, "analytic")) {
      interval_text <- "analytical confidence interval fallback"
    } else if (setequal(int_nonrep, "model")) {
      interval_text <- "model-based credible interval fallback from the Bayesian trivariate meta-analysis"
    } else {
      interval_text <- "fallback intervals"
    }

    if (isTRUE(mark_imputed)) {
      message("* indicates ", interval_text,
              " when reported study intervals were unavailable.")
    } else {
      message(int_leadin, " intervals were replaced by ", interval_text,
              " because reported study intervals were unavailable.")
    }
  }
}


# helper for reproducibility by building inits
resolve_jags_inits <- function(inits = NULL,
                               n.chains,
                               seed = NULL,
                               rng_name = "base::Wichmann-Hill") {

  valid_rng <- c(
    "base::Wichmann-Hill",
    "base::Marsaglia-Multicarry",
    "base::Super-Duper",
    "base::Mersenne-Twister"
  )

  if (!is.character(rng_name) || length(rng_name) != 1 || !rng_name %in% valid_rng) {
    stop(
      "`rng_name` must be one of: ",
      paste(valid_rng, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(n.chains) || length(n.chains) != 1 || n.chains < 1) {
    stop("`n.chains` must be a positive integer.", call. = FALSE)
  }
  n.chains <- as.integer(n.chains)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
      stop("`seed` must be a single non-missing numeric value.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  # ---------- case 1: no inits supplied ----------
  if (is.null(inits)) {
    if (is.null(seed)) {
      return(NULL)
    }

    # reproducible per-chain RNG settings
    out <- vector("list", n.chains)
    for (k in seq_len(n.chains)) {
      out[[k]] <- list(
        .RNG.name = rng_name,
        .RNG.seed = seed + (k - 1L)
      )
    }
    return(out)
  }

  # ---------- validate inits ----------
  if (!is.list(inits) || length(inits) != n.chains) {
    stop(
      "`inits` must be a list of length `n.chains` (", n.chains, "), ",
      "with one sub-list per chain.",
      call. = FALSE
    )
  }

  if (!all(vapply(inits, is.list, logical(1)))) {
    stop("Each element of `inits` must itself be a list.", call. = FALSE)
  }

  out <- inits

  has_rng_name <- vapply(out, function(x) ".RNG.name" %in% names(x), logical(1))
  has_rng_seed <- vapply(out, function(x) ".RNG.seed" %in% names(x), logical(1))
  has_full_rng <- has_rng_name & has_rng_seed

  # ---------- case 2: inits supplied, no seed ----------
  if (is.null(seed)) {
    if (!all(has_full_rng)) {
      warning(
        "Some chains in `inits` do not contain both `.RNG.name` and `.RNG.seed`. ",
        "Results may not be exactly reproducible.",
        call. = FALSE
      )
    }
    return(out)
  }

  # ---------- case 3: inits supplied, seed supplied ----------
  # Fill only missing RNG fields; never overwrite user-supplied parameter initials.
  # Also never overwrite user-supplied RNG fields.
  all_full_before <- all(has_full_rng)

  for (k in seq_len(n.chains)) {
    if (!has_rng_name[k]) {
      out[[k]]$.RNG.name <- rng_name
    }
    if (!has_rng_seed[k]) {
      out[[k]]$.RNG.seed <- seed + (k - 1L)
    }
  }

  if (all_full_before) {
    message(
      "`seed` ignored because `.RNG.name` and `.RNG.seed` were already supplied ",
      "in `inits` for all chains."
    )
  } else {
    message(
      "Filled missing JAGS RNG settings in `inits` using `seed = ", seed, "`."
    )
  }

  out
}

# extract per-chain JAGS RNG metadata from resolved init lists
extract_rng_meta <- function(final_inits) {
  get_rng_name <- function(x) {
    if (".RNG.name" %in% names(x)) x$.RNG.name else NA_character_
  }

  get_rng_seed <- function(x) {
    if (".RNG.seed" %in% names(x)) as.integer(x$.RNG.seed) else NA_integer_
  }

  if (is.null(final_inits)) {
    return(list(
      rng_name_per_chain = NULL,
      rng_seed_per_chain = NULL
    ))
  }

  list(
    rng_name_per_chain = vapply(final_inits, get_rng_name, character(1)),
    rng_seed_per_chain = vapply(final_inits, get_rng_seed, integer(1))
  )
}

####################------------------- helpers for MA_NB_tri_voi.R ----------------------######################
# convert mcmc.list to tibble of draws, a helper for compute_voi_metrics()
mcmc_list_to_draws <- function(samples) {
  if (!inherits(samples, "mcmc.list")) {
    stop("`samples` must be a coda::mcmc.list.", call. = FALSE)
  }
  # row-bind all chains; keep parameter names as columns
  mats <- lapply(samples, as.matrix)
  draws <- do.call(rbind, mats)

  # deduplicate columns defensively before coercing to data frame
  # This is needed because expand_target in MA_NB_tri() expands to per-study NB[]s, but center_rows() also requests specific per-study NB[]s
  dups <- colnames(draws)[duplicated(colnames(draws))]
  if (length(dups) > 0) {
    message("Removing duplicate columns from mcmc matrix: ",
            paste(unique(dups), collapse = ", "))
    draws <- draws[, !duplicated(colnames(draws)), drop = FALSE]
  }

  # safer as tibble/data.frame; preserve numeric columns
  draws <- as.data.frame(draws, check.names = FALSE)
  tibble::as_tibble(draws)
}

# compute sigma values for each pair of columns in df, MCMC diagnostic helper
calc_sigma <- function(df) {
  df <- as.data.frame(df)
  n <- nrow(df)
  m <- ncol(df)
  if (n < 2 || m < 2) return(numeric(0))

  out <- numeric(m * (m - 1) / 2)
  idx <- 1
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      d <- df[[i]] - df[[j]]
      se <- stats::sd(d) / sqrt(n)

      if (is.na(se) || se == 0) {
        out[idx] <- Inf
      } else {
        out[idx] <- abs(mean(d)) / se
      }

      idx <- idx + 1
    }
  }
  out
}

# for the precision/ESS stopping rule
append_mcmc_list <- function(old, new) {
  stopifnot(inherits(new, "mcmc.list"))
  if (is.null(old)) return(new)
  stopifnot(inherits(old, "mcmc.list"))
  if (length(old) != length(new)) stop("Different number of chains.")

  out <- vector("list", length(old))
  for (k in seq_along(old)) {
    m_old <- as.matrix(old[[k]])
    m_new <- as.matrix(new[[k]])

    # Safety: align columns in case ordering differs
    common <- intersect(colnames(m_old), colnames(m_new))
    if (length(common) == 0) stop("No common variables between old and new samples.")
    m_old <- m_old[, common, drop = FALSE]
    m_new <- m_new[, common, drop = FALSE]

    out[[k]] <- coda::mcmc(rbind(m_old, m_new))
  }
  coda::mcmc.list(out)
}


