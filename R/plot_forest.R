plot_forest_metric_forestploter <- function(
    sum,
    metric = c("NB", "RU"),
    center = c("Median", "Mean"),
    xlim = NULL,
    plot_col_width = 50,
    file_png = NULL,
    file_pdf = NULL
) {
  metric <- match.arg(metric)
  sum_m  <- sum[[metric]]

  if (is.null(sum_m$per_study)) {
    stop("No per-study results found for metric ", metric,
         ". Did you return RU[i] in MA_NB_tri() and summarize it?")
  }

  if (is.null(sum_m)) {
    stop("Metric ", metric, " not found in summary object.")
  }

  center_missing <- missing(center)
  center <- match.arg(center)

  if (center_missing) {
    message("Using center = '", center, "' (default). Set center = 'Mean' if desired.")
  } else {
    message("Using center = '", center, "'.")
  }

  # P(useful) (NB only)
  prob_useful <- NULL
  if (metric == "NB" && "probuseful" %in% names(sum)) {
    prob_useful <- sum$probuseful$pooled$model["Mean"]
  }
  show_pu <- (metric == "NB") && !is.null(prob_useful)

  # ---------- LABELS ----------
  xlab <- if (metric == "NB") "Net Benefit" else "Relative Utility"
  col_label <- if (metric == "NB") "NB (95% CrI)" else "RU (95% CrI)"

  # ---------- CHECK IF COUNTRY COLUMN IS PRESENT ----------
  has_country <- "Country" %in% names(sum_m$per_study)

  # ---------- 1. Per-study table (sorted by mean) ----------
  per <- sum_m$per_study %>%
    dplyr::arrange(.data[[center]])

  df <- per %>%
    dplyr::mutate(
      Publication = Publication,
      Country     = if (has_country) Country else NULL,
      N           = as.character(N),
      Prev        = scales::percent(Prev, accuracy = 1),
      ` `         = strrep(" ", plot_col_width),
      value_ci    = sprintf("%5.2f  [%.2f; %.2f]", .data[[center]], Low, High)
    )

  df <- if (has_country) {
    df %>% dplyr::select(Publication, Country, N, Prev, ` `, value_ci)
  } else {
    df %>% dplyr::select(Publication, N, Prev, ` `, value_ci)
  }

  names(df)[ncol(df)] <- col_label

  # ---------- 2. Pooled row (model strategy only) ----------
  pooled <- sum_m$pooled$model

  pooled_row <- tibble::tibble(
    Publication = "Pooled estimate",
    Country     = if (has_country) "" else NULL,
    N = "", Prev = "",
    ` ` = strrep(" ", plot_col_width),
    value_ci = sprintf("%5.2f  [%.2f; %.2f]",
                       pooled[center], pooled["Low"], pooled["High"])
  )

  if (!has_country) pooled_row$Country <- NULL
  names(pooled_row)[ncol(pooled_row)] <- col_label

  # ---------- 3. Prediction interval row ----------
  pred <- sum_m$predictive$model

  pred_row <- tibble::tibble(
    Publication = "Prediction interval",
    Country     = if (has_country) "" else NULL,
    N = "", Prev = "",
    ` ` = strrep(" ", plot_col_width),
    value_ci = sprintf("[%.2f; %.2f]", pred["Low"], pred["High"])
  )

  if (!has_country) pred_row$Country <- NULL
  names(pred_row)[ncol(pred_row)] <- col_label

  # ---------- 3B. P(useful) row (text only) ----------
  pu_row <- NULL
  if (show_pu) {
    pu_row <- tibble::tibble(
      Publication = "P(useful)",
      Country     = if (has_country) "" else NULL,
      N = "", Prev = "",
      ` ` = strrep(" ", plot_col_width),
      value_ci = sprintf("%5.2f", prob_useful)
    )
    if (!has_country) pu_row$Country <- NULL
    names(pu_row)[ncol(pu_row)] <- col_label
  }



  # ---------- 4. Combine table ----------
  tbl <- dplyr::bind_rows(df, pooled_row, pred_row, pu_row) %>%
    dplyr::mutate(dplyr::across(everything(), as.character))

  # ---------- 5. Numeric vectors ----------
  # base vectors (studies + pooled + prediction interval)
  est   <- c(per[[center]], pooled[center], NA)
  lower <- c(per$Low,       pooled["Low"],  NA)
  upper <- c(per$High,      pooled["High"], NA)

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

  # ensure 0 is visible (but do NOT force axis to start at 0)
  xlim[1] <- min(xlim[1], 0, na.rm = TRUE)
  xlim[2] <- max(xlim[2], 0, na.rm = TRUE)

  # ---------- 6. Theme ----------
  tm <- forestploter::forest_theme(
    base_size    = 10,
    refline_gp   = grid::gpar(lty = "solid"),
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
    ref_line   = 0,
    vert_line = pooled[center],
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
