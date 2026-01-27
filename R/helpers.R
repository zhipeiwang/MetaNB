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



# summarizer for forest plot function
summarize_mcmc_outputs <- function(
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
    )
  )

  out <- list()

  for (tar in targets) {

    # -------------------------------
    # CASE 1: family target (NB / RU / probuseful)
    # -------------------------------
    if (tar %in% names(family_spec)) {

      spec <- family_spec[[tar]]

      # ---- per-study ----
      per_study <- NULL
      if (!is.null(spec$per_study_base) && tar %in% targets_per_study) {

        pat  <- paste0("^", spec$per_study_base, "\\[\\d+\\]$")
        rows <- grep(pat, rownames(stats), value = TRUE)

        if (length(rows) > 0) {
          if (is.null(study_info)) {
            stop("data required for per-study summaries of ", tar)
          }

          idx <- as.integer(sub(paste0("^", spec$per_study_base, "\\[(\\d+)\\]$"),
                                "\\1", rows))
          ord <- order(idx)

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
      if (!is.null(spec$pooled)) {
        for (k in names(spec$pooled)) {
          pooled[[k]] <- get_cri(spec$pooled[[k]])
        }
        pooled <- pooled[!sapply(pooled, is.null)]
        if (length(pooled) == 0) pooled <- NULL
      }

      # ---- predictive (PI) ----
      predictive <- list()
      if (!is.null(spec$predictive)) {
        for (k in names(spec$predictive)) {
          predictive[[k]] <- get_pi(spec$predictive[[k]])
        }
        predictive <- predictive[!sapply(predictive, is.null)]
        if (length(predictive) == 0) predictive <- NULL
      }

      # ---- pooled_ref (CrI) ----
      pooled_ref <- NULL
      if (return_ref && !is.null(spec$pooled_ref)) {
        pooled_ref <- list()
        for (k in names(spec$pooled_ref)) {
          pooled_ref[[k]] <- get_cri(spec$pooled_ref[[k]])
        }
        pooled_ref <- pooled_ref[!sapply(pooled_ref, is.null)]
        if (length(pooled_ref) == 0) pooled_ref <- NULL
      }

      # ---- pred_ref (PI or mean-only) ----
      pred_ref <- NULL
      if (return_ref && !is.null(spec$pred_ref)) {
        pred_ref <- list()
        for (k in names(spec$pred_ref)) {
          pred_ref[[k]] <- get_pi(spec$pred_ref[[k]])
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

  out
}
