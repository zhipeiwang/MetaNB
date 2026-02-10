plot_forest_metric_forestploter <- function(
    sum,
    label_cols = c("Publication", "Country", "N", "Prev"),
    study_label_col = NULL,
    prev_col = NULL,
    n_col = NULL,
    metric = c("NB", "RU", "sens", "spec"),
    center = c("Median", "Mean"),
    xlim = NULL,
    xticks = NULL,
    plot_col_width = 50,
    data = NULL,
    use_reported_ci = FALSE,
    reported_low_col = NULL,
    reported_high_col = NULL,
    reported_est_col = NULL,
    file_png = NULL,
    file_pdf = NULL
) {
  xlim_user <- !is.null(xlim)

  metric <- match.arg(metric)
  sum_m  <- sum[[metric]]

  if (is.null(sum_m)) {
    stop("Metric ", metric, " not found in summary object.")
  }
  if (is.null(sum_m$per_study)) {
    stop("No per-study results found for metric ", metric,
         ". Did you return it in MA_NB_tri() and summarize it?")
  }


  center_missing <- missing(center)
  center <- match.arg(center)

  if (!(metric %in% c("sens","spec") && isTRUE(use_reported_ci))) {
    if (center_missing) {
      message("Using center = '", center, "' (default). Set center = 'Mean' if desired.")
    } else {
      message("Using center = '", center, "'.")
    }
  } else {
    message("Using observed point estimates for ", metric, " (center ignored).")
  }


  # P(useful) (NB only)
  prob_useful <- NULL
  if (metric == "NB" && "probuseful" %in% names(sum)) {
    prob_useful <- sum$probuseful$pooled$model["Mean"]
  }
  show_pu <- (metric == "NB") && !is.null(prob_useful)

  # ---------- LABELS ----------
  label_map <- list(
    NB   = list(xlab = "Net Benefit",        col = "NB (95% CrI)"),
    RU   = list(xlab = "Relative Utility",   col = "RU (95% CrI)"),
    sens = list(xlab = "Sensitivity",        col = NULL),
    spec = list(xlab = "Specificity",        col = NULL)
  )

  if (!metric %in% names(label_map)) {
    stop("Unknown metric: ", metric)
  }

  xlab <- label_map[[metric]]$xlab

  # interval label for sens/spec depends on use_reported_ci
  if (metric %in% c("sens", "spec")) {
    interval_label <- if (isTRUE(use_reported_ci)) "CI" else "CrI"
    col_label <- if (metric == "sens") {
      paste0("Sensitivity (95% ", interval_label, ")")
    } else {
      paste0("Specificity (95% ", interval_label, ")")
    }
  } else {
    col_label <- label_map[[metric]]$col
  }



  # keep only label columns that exist
  label_cols <- intersect(label_cols, names(sum_m$per_study))
  if (length(label_cols) == 0) stop("None of `label_cols` exist in per-study results.", call. = FALSE)

  # which column should receive "Pooled estimate", etc.
  if (is.null(study_label_col)) {
    study_label_col <- label_cols[1]
  } else {
    if (!study_label_col %in% names(sum_m$per_study)) {
      stop("`study_label_col` = '", study_label_col, "' not found in per-study results.", call. = FALSE)
    }
    # ensure it's included in the table
    if (!study_label_col %in% label_cols) label_cols <- c(study_label_col, label_cols)
  }


  # ---------- 1. Per-study table (sorted by point estimate) ----------
  per <- sum_m$per_study

  # ---------- 1B. Optional: use observed sens/spec and reported CIs ----------
  display_est <- per[[center]]
  display_low <- per$Low
  display_high <- per$High
  used_reported_ci <- rep(FALSE, nrow(per))

  if (isTRUE(use_reported_ci) && metric %in% c("sens", "spec")) {

    # Use observed point estimate if available
    if ("Observed" %in% names(per)) {
      display_est <- per$Observed
    } else if (!is.null(data) &&
               all(c("tp","tn","n_event","n_nonevent") %in% names(data))) {
      if (metric == "sens") display_est <- data$tp / data$n_event
      if (metric == "spec") display_est <- data$tn / data$n_nonevent
    }

    # Use reported CI if provided; fallback remains Bayesian CrI
    if (!is.null(data) &&
        !is.null(reported_low_col) && !is.null(reported_high_col) &&
        reported_low_col %in% names(data) && reported_high_col %in% names(data)) {

      rep_low  <- suppressWarnings(as.numeric(data[[reported_low_col]]))
      rep_high <- suppressWarnings(as.numeric(data[[reported_high_col]]))

      has_rep <- !is.na(rep_low) & !is.na(rep_high)

      display_low[has_rep] <- rep_low[has_rep]
      display_high[has_rep] <- rep_high[has_rep]
      used_reported_ci[has_rep] <- TRUE
    }

    # Optional: use reported point estimate if user supplies it
    if (!is.null(data) && !is.null(reported_est_col) &&
        reported_est_col %in% names(data)) {
      rep_est <- suppressWarnings(as.numeric(data[[reported_est_col]]))
      has_est <- !is.na(rep_est)
      display_est[has_est] <- rep_est[has_est]
    }

    # Sort by what you're actually plotting
    ord <- order(display_est)
    per <- per[ord, , drop = FALSE]
    display_est <- display_est[ord]
    display_low <- display_low[ord]
    display_high <- display_high[ord]
    used_reported_ci <- used_reported_ci[ord]
  } else {
    # original sorting for NB/RU or Bayesian-only sens/spec
    per <- per %>% dplyr::arrange(.data[[center]])
  }


  # Add a flag column to indicate which rows used reported CIs (for potential styling later)
  df <- per %>%
    dplyr::mutate(
      .display_est = display_est,
      .display_low = display_low,
      .display_high = display_high,
      ` `      = strrep(" ", plot_col_width),
      value_ci = sprintf("%5.2f  [%.2f; %.2f]", .display_est, .display_low, .display_high)
    ) %>%
    dplyr::select(dplyr::all_of(label_cols), ` `, value_ci)


  # pretty formatting (only if user tells us which columns these are)
  if (!is.null(prev_col) && prev_col %in% names(df)) {
    df <- df %>% dplyr::mutate(!!prev_col := scales::percent(as.numeric(.data[[prev_col]]), accuracy = 1))
  }

  if (!is.null(n_col) && n_col %in% names(df)) {
    df <- df %>% dplyr::mutate(!!n_col := format(as.numeric(.data[[n_col]]),
                                                 big.mark = ",", scientific = FALSE, trim = TRUE))
  }

  # IMPORTANT: make label columns character so bind_rows works with blank rows
  df <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(label_cols), as.character))


  names(df)[ncol(df)] <- col_label

  # build a "blank row" template
  make_row <- function(title, value_ci) {
    row <- as.list(rep("", length(label_cols)))
    names(row) <- label_cols
    row[[study_label_col]] <- title

    row <- tibble::as_tibble(row) %>%
      dplyr::mutate(
        ` ` = strrep(" ", plot_col_width),
        value_ci = value_ci
      )

    names(row)[ncol(row)] <- col_label
    row
  }

  # Pooled and predictive objects
  pooled <- sum_m$pooled$model
  pred   <- sum_m$predictive$model

  if (is.null(pooled)) stop("No pooled results found for metric ", metric, call. = FALSE)
  if (is.null(pred))   stop("No predictive results found for metric ", metric, call. = FALSE)

  # ---------- 2. Pooled row (model strategy only) ----------
  pooled_row <- make_row(
    "Pooled estimate",
    sprintf("%5.2f  [%.2f; %.2f]", pooled[center], pooled["Low"], pooled["High"])
  )

  # ---------- 3. Prediction interval row ----------
  pred_row <- make_row(
    "Prediction interval",
    sprintf("[%.2f; %.2f]", pred["Low"], pred["High"])
  )

  # ---------- 3B. P(useful) row (text only) ----------
  pu_row <- if (show_pu) make_row("P(useful)", sprintf("%.2f", prob_useful)) else NULL


  # ---------- 4. Combine table ----------
  tbl <- dplyr::bind_rows(df, pooled_row, pred_row, pu_row) %>%
    dplyr::mutate(dplyr::across(everything(), as.character))

  # ---------- 5. Numeric vectors ----------
  # base vectors (studies + pooled + prediction interval)
  est   <- c(display_est,  pooled[center], NA)
  lower <- c(display_low,  pooled["Low"],  NA)
  upper <- c(display_high, pooled["High"], NA)


  is_sum <- c(
    rep(FALSE, nrow(per)),
    TRUE,
    FALSE
  )

  # add one more entry ONLY if P(useful) row exists
  if (show_pu) {
    est   <- c(est,   NA)
    lower <- c(lower, NA)
    upper <- c(upper, NA)
    is_sum <- c(is_sum, FALSE)
  }
  ci_col <- which(names(tbl) == " ")

  # ---------- 5B. X-limits (avoid truncation) ----------
  if (is.null(xlim)) {
    # Include study CIs + pooled CI + prediction interval
    all_limits <- c(lower, upper)

    # pred is a named vector with Low/High in your summarizer
    all_limits <- c(all_limits, pred["Low"], pred["High"])

    xlim <- range(all_limits, na.rm = TRUE)

    # add a bit of padding
    pad <- 0.05 * diff(xlim)
    if (is.finite(pad) && pad > 0) xlim <- xlim + c(-pad, pad)
  }

  if (!xlim_user) {
    # ensure 0 is visible (but do NOT force axis to start at 0)
    xlim[1] <- min(xlim[1], 0, na.rm = TRUE)
    xlim[2] <- max(xlim[2], 0, na.rm = TRUE)
  }


  # ---------- 6. Theme ----------
  tm <- forestploter::forest_theme(
    base_size    = 10,
    # refline_gp   = grid::gpar(lty = "solid"),
    summary_fill = "grey20",
    summary_col  = "grey20",
    point_size   = 2
  )

  # ---------- 7. Draw forest ----------
  p <- forestploter::forest(
    tbl,
    est        = est,
    lower      = lower,
    upper      = upper,
    is_summary = is_sum,
    ci_column  = ci_col,
    xlab       = xlab,
    xlim       = xlim,
    ticks_at   = xticks,
    ref_line   = pooled[center],
    # vert_line = pooled[center],
    theme      = tm
  )

  # ---------- 7B. Prediction interval thick bar ----------
  npc_x <- function(x) (x - xlim[1]) / diff(xlim)

  p <- forestploter::add_grob(
    p,
    row = nrow(per) + 2, # the prediction interval row is always after the pooled row
    col = ci_col,
    gb_fn = grid::segmentsGrob,
    x0 = grid::unit(npc_x(pred["Low"]),  "npc"),
    x1 = grid::unit(npc_x(pred["High"]), "npc"),
    y0 = grid::unit(0.5, "npc"),
    y1 = grid::unit(0.5, "npc"),
    gp = grid::gpar(lwd = 3)
  )

  # ---------- Display ----------
  plot(p)

  # ---------- Save ----------
  wh <- forestploter::get_wh(p, unit = "in")

  if (!is.null(file_png)) {
    png(file_png, width = wh[1], height = wh[2],
        units = "in", res = 300, bg = "white")
    plot(p)
    dev.off()
  }

  if (!is.null(file_pdf)) {
    grDevices::cairo_pdf(file_pdf, width = wh[1], height = wh[2], bg = "white")
    plot(p)
    dev.off()
  }

  invisible(p)
}
