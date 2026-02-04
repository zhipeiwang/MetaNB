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
summarize_for_forestplot <- function(
    samples,
    data = NULL,
    label_cols = c("Publication", "Country", "N", "Prev"),
    targets = c("NB", "probuseful"),
    targets_per_study = c("NB"),
    return_ref = FALSE
) {
  study_info <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)

    missing_labels <- setdiff(label_cols, names(data))
    if (length(missing_labels)) {
      warning("Missing label_cols ignored: ", paste(missing_labels, collapse = ", "))
    }
    label_cols <- intersect(label_cols, names(data))
    study_info <- data[, label_cols, drop = FALSE]
  }

  samples <- get_samples(samples)
  ss <- summary(samples)
  stats  <- ss$statistics
  quants <- ss$quantiles

  # ---- helpers ----
  get_cri <- function(name) {
    if (!name %in% rownames(stats)) return(NULL)
    c(
      Mean = stats[name, "Mean"],
      Median = quants[name, "50%"],
      Low  = quants[name, "2.5%"],
      High = quants[name, "97.5%"]
    )
  }

  get_pi <- function(name) {
    if (!name %in% rownames(stats)) return(NULL)
    c(
      Low  = quants[name, "2.5%"],
      High = quants[name, "97.5%"]
    )
  }

  get_mean_only <- function(name) {
    if (!name %in% rownames(stats)) return(NULL)
    c(Mean = stats[name, "Mean"])
  }

  # ---- family specification ----
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

      pooled_ref = list(
        model     = "pooledNB_ref",
        treat_all = "pooledNB_TA_ref"
      ),

      pred_ref = list(
        model     = "NBnew_ref",
        treat_all = "NBnew_TA_ref"
      )
    ),

    RU = list(
      per_study_base = "RU",
      pooled = list(model = "pooledRU"),
      predictive = list(model = "RUnew"),
      pooled_ref = list(model = "pooledRU_ref"),
      pred_ref = list(model = "RUnew_ref")
    ),

    probuseful = list(
      per_study_base = NULL,
      pooled = list(model = "probuseful"),
      predictive = NULL,
      pooled_ref = list(model = "probuseful_ref"),
      pred_ref = NULL
    ),
    sens = list(
      per_study_base = "sens",
      pooled = list(model = "pooledsens"),
      predictive = list(model = "sensnew"),
      pooled_ref = NULL,
      pred_ref = NULL
    ),

    spec = list(
      per_study_base = "spec",
      pooled = list(model = "pooledspec"),
      predictive = list(model = "specnew"),
      pooled_ref = NULL,
      pred_ref = NULL
    )

  )

  out <- list()

  for (tar in targets) {

    # -------------------------------
    # CASE 1: family target (NB / RU / probuseful)
    # -------------------------------
    if (tar %in% names(family_spec)) {

      metric_def <- family_spec[[tar]]

      # ---- per-study ----
      per_study <- NULL
      if (!is.null(metric_def$per_study_base) && tar %in% targets_per_study) {

        pat  <- paste0("^", metric_def$per_study_base, "\\[\\d+\\]$")
        rows <- grep(pat, rownames(stats), value = TRUE)

        if (length(rows) > 0) {
          if (is.null(study_info)) {
            stop("data required for per-study summaries of ", tar)
          }

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
          pooled[[k]] <- get_cri(metric_def$pooled[[k]])
        }
        pooled <- pooled[!sapply(pooled, is.null)]
        if (length(pooled) == 0) pooled <- NULL
      }

      # ---- predictive (PI) ----
      predictive <- list()
      if (!is.null(metric_def$predictive)) {
        for (k in names(metric_def$predictive)) {
          predictive[[k]] <- get_pi(metric_def$predictive[[k]])
        }
        predictive <- predictive[!sapply(predictive, is.null)]
        if (length(predictive) == 0) predictive <- NULL
      }

      # ---- pooled_ref (CrI) ----
      pooled_ref <- NULL
      if (return_ref && !is.null(metric_def$pooled_ref)) {
        pooled_ref <- list()
        for (k in names(metric_def$pooled_ref)) {
          pooled_ref[[k]] <- get_cri(metric_def$pooled_ref[[k]])
        }
        pooled_ref <- pooled_ref[!sapply(pooled_ref, is.null)]
        if (length(pooled_ref) == 0) pooled_ref <- NULL
      }

      # ---- pred_ref (PI or mean-only) ----
      pred_ref <- NULL
      if (return_ref && !is.null(metric_def$pred_ref)) {
        pred_ref <- list()
        for (k in names(metric_def$pred_ref)) {
          pred_ref[[k]] <- get_pi(metric_def$pred_ref[[k]])
        }
        pred_ref <- pred_ref[!sapply(pred_ref, is.null)]
        if (length(pred_ref) == 0) pred_ref <- NULL
      }

      # ---- probuseful special case (mean only) ----
      if (tar == "probuseful") {
        pooled <- list(model = get_mean_only("probuseful"))
        if (return_ref) {
          pooled_ref <- list(model = get_mean_only("probuseful_ref"))
        }
        predictive <- NULL
        pred_ref   <- NULL
      }

      out[[tar]] <- list(
        per_study  = per_study,
        pooled     = pooled,
        predictive = predictive,
        pooled_ref = pooled_ref,
        pred_ref   = pred_ref
      )

    } else {

      # -------------------------------
      # CASE 2: atomic variable name
      # -------------------------------
      scalar <- NULL
      if (tar %in% rownames(stats)) {
        if (tar %in% rownames(quants)) {
          scalar <- get_cri(tar)
        } else {
          scalar <- get_mean_only(tar)
        }
      }

      out[[tar]] <- list(scalar = scalar)
    }
  }

  # ---- drop NULLs recursively ----
  drop_nulls_recursive(out)


}


####################------------------- user-facing summarizer ----------------------######################
summarize_for_users <- function(
    samples,
    data = NULL,
    label_cols = c("Publication", "Country", "N", "Prev"),
    metrics = c("NB", "RU", "probuseful"),
    per_study_metrics = c("NB", "RU"),
    return_ref = FALSE,
    include_per_study = TRUE
) {
  core <- summarize_for_forestplot(
    samples = get_samples(samples),
    data = data,
    label_cols = label_cols,
    targets = metrics,
    targets_per_study = if (include_per_study) per_study_metrics else character(0),
    return_ref = return_ref
  )

  out <- list()

  # --- NB / RU: keep structure; named vectors (no 1-row dfs) ---
  for (m in intersect(c("NB","RU"), names(core))) {
    obj <- core[[m]]

    out[[m]] <- drop_nulls_recursive(list(
      per_study  = if (include_per_study) obj$per_study else NULL,
      pooled     = obj$pooled,        # named vectors inside
      predictive = obj$predictive,    # named vectors inside
      pooled_ref = if (return_ref) obj$pooled_ref else NULL,
      pred_ref   = if (return_ref) obj$pred_ref else NULL
    ))
  }

  # --- probuseful: scalar value, not pooled/model ---
  if ("probuseful" %in% names(core)) {
    pu <- core$probuseful

    value <- NULL
    if (!is.null(pu$pooled$model) && "Mean" %in% names(pu$pooled$model)) {
      value <- unname(pu$pooled$model["Mean"])
    }

    value_ref <- NULL
    if (return_ref && !is.null(pu$pooled_ref$model) && "Mean" %in% names(pu$pooled_ref$model)) {
      value_ref <- unname(pu$pooled_ref$model["Mean"])
    }

    out$probuseful <- drop_nulls_recursive(list(
      value = value,
      value_ref = value_ref
    ))
  }

  # atomic scalars (if you allow them in metrics)
  other <- setdiff(metrics, c("NB","RU","probuseful"))
  for (nm in intersect(other, names(core))) {
    out[[nm]] <- core[[nm]]
  }

  drop_nulls_recursive(out)
}


# convert mcmc.list to tibble of draws, a helper for compute_voi_metrics()
mcmc_list_to_draws <- function(samples) {
  if (!inherits(samples, "mcmc.list")) {
    stop("`samples` must be a coda::mcmc.list.", call. = FALSE)
  }
  # row-bind all chains; keep parameter names as columns
  mats <- lapply(samples, as.matrix)
  draws <- do.call(rbind, mats)

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
      out[idx] <- abs(mean(d)) / se
      idx <- idx + 1
    }
  }
  out
}


