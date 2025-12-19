plot_nb_forest_forestploter <- function(sum_nb,
                                        xlim = c(-0.1, 0.7),
                                        file_png = NULL,
                                        file_pdf = NULL) {

  # ---------- CHECK IF COUNTRY COLUMN IS PRESENT ----------
  has_country <- "Country" %in% names(sum_nb$per_study)

  # ---------- 1. Per-study table (sorted by NB) ----------
  per <- sum_nb$per_study %>%
    arrange(NB)

  df <- per %>%
    mutate(
      Publication = Publication,
      Country     = if (has_country) Country else NULL,
      N           = as.character(N),
      Prev        = percent(Prev, accuracy = 1),
      ` `         = strrep(" ", 40),  # blank column for CI panel
      `NB (95% CrI)` = sprintf("%.2f  [%.2f; %.2f]", NB, NB_low, NB_high)
    )

  # dynamically choose columns based on availability of Country
  df <- if (has_country) {
    df %>% select(Publication, Country, N, Prev, ` `, `NB (95% CrI)`)
  } else {
    df %>% select(Publication, N, Prev, ` `, `NB (95% CrI)`)
  }

  # ---------- 2. Pooled row ----------
  pooled_mean <- unname(sum_nb$pooled["mean"])
  pooled_low  <- unname(sum_nb$pooled["low"])
  pooled_high <- unname(sum_nb$pooled["high"])

  pooled_row <- tibble(
    Publication = "Random effects model",
    Country     = if (has_country) "" else NULL,
    N = "", Prev = "",
    ` ` = strrep(" ", 40),
    `NB (95% CrI)` =
      sprintf("%.2f  [%.2f; %.2f]", pooled_mean, pooled_low, pooled_high)
  )

  if (!has_country) pooled_row$Country <- NULL

  # ---------- 3. Prediction interval row ----------
  pred_low  <- unname(sum_nb$pred["low"])
  pred_high <- unname(sum_nb$pred["high"])

  pred_row <- tibble(
    Publication = "Prediction interval",
    Country     = if (has_country) "" else NULL,
    N = "", Prev = "",
    ` ` = strrep(" ", 40),
    `NB (95% CrI)` = sprintf("[%.2f; %.2f]", pred_low, pred_high)
  )

  if (!has_country) pred_row$Country <- NULL

  # ---------- 4. Combine table ----------
  tbl <- bind_rows(df, pooled_row, pred_row) %>%
    mutate(across(everything(), as.character))

  # ---------- 5. Numeric vectors ----------
  est   <- c(per$NB, pooled_mean, NA)    # NA = hide pred interval point
  lower <- c(per$NB_low, pooled_low, NA)
  upper <- c(per$NB_high, pooled_high, NA)

  is_sum <- c(rep(FALSE, nrow(per)), TRUE, FALSE)

  ci_col <- which(names(tbl) == " ")

  # ---------- 6. Theme ----------
  tm <- forest_theme(
    base_size    = 10,
    summary_fill = "grey20",
    summary_col  = "grey20",
    point_size   = 2
  )

  # ---------- 7. Draw main forest ----------
  p <- forestploter::forest(
    tbl,
    est        = est,
    lower      = lower,
    upper      = upper,
    is_summary = is_sum,
    ci_column  = ci_col,
    xlab       = "Net Benefit",
    xlim       = xlim,
    ref_line   = 0,
    theme      = tm
  )

  # ---------- 7B. Prediction interval thick bar ----------
  npc_x <- function(x) (x - xlim[1]) / diff(xlim)

  p <- add_grob(
    p,
    row = nrow(tbl), col = ci_col,
    gb_fn = segmentsGrob,
    x0 = unit(npc_x(pred_low),  "npc"),
    x1 = unit(npc_x(pred_high), "npc"),
    y0 = unit(0.5, "npc"),
    y1 = unit(0.5, "npc"),
    gp = gpar(lwd = 6)
  )

  # ---------- Display ----------
  plot(p)

  # ---------- 8. Save ----------
  wh <- get_wh(p, unit = "in")

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


plot_forest_metric_forestploter <- function(
    sum,
    metric = c("NB", "RU"),
    xlim = NULL,
    file_png = NULL,
    file_pdf = NULL
) {
  metric <- match.arg(metric)
  sum_m  <- sum[[metric]]

  if (is.null(sum_m)) {
    stop("Metric ", metric, " not found in summary object.")
  }

  # ---------- LABELS ----------
  xlab <- if (metric == "NB") "Net Benefit" else "Relative Utility"
  col_label <- if (metric == "NB") "NB (95% CrI)" else "RU (95% CrI)"

  # ---------- CHECK IF COUNTRY COLUMN IS PRESENT ----------
  has_country <- "Country" %in% names(sum_m$per_study)

  # ---------- 1. Per-study table (sorted by mean) ----------
  per <- sum_m$per_study %>%
    dplyr::arrange(mean)

  df <- per %>%
    dplyr::mutate(
      Publication = Publication,
      Country     = if (has_country) Country else NULL,
      N           = as.character(N),
      Prev        = scales::percent(Prev, accuracy = 1),
      ` `         = strrep(" ", 40),
      value_ci    = sprintf("%.2f  [%.2f; %.2f]", mean, low, high)
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
    Publication = "Random effects model",
    Country     = if (has_country) "" else NULL,
    N = "", Prev = "",
    ` ` = strrep(" ", 40),
    value_ci = sprintf("%.2f  [%.2f; %.2f]",
                       pooled["mean"], pooled["low"], pooled["high"])
  )

  if (!has_country) pooled_row$Country <- NULL
  names(pooled_row)[ncol(pooled_row)] <- col_label

  # ---------- 3. Prediction interval row ----------
  pred <- sum_m$predictive$model

  pred_row <- tibble::tibble(
    Publication = "Prediction interval",
    Country     = if (has_country) "" else NULL,
    N = "", Prev = "",
    ` ` = strrep(" ", 40),
    value_ci = sprintf("[%.2f; %.2f]", pred["low"], pred["high"])
  )

  if (!has_country) pred_row$Country <- NULL
  names(pred_row)[ncol(pred_row)] <- col_label

  # ---------- 4. Combine table ----------
  tbl <- dplyr::bind_rows(df, pooled_row, pred_row) %>%
    dplyr::mutate(dplyr::across(everything(), as.character))

  # ---------- 5. Numeric vectors ----------
  est   <- c(per$mean, pooled["mean"], NA)
  lower <- c(per$low,  pooled["low"],  NA)
  upper <- c(per$high, pooled["high"], NA)

  is_sum <- c(rep(FALSE, nrow(per)), TRUE, FALSE)
  ci_col <- which(names(tbl) == " ")

  # ---------- 6. Theme ----------
  tm <- forestploter::forest_theme(
    base_size    = 10,
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
    theme      = tm
  )

  # ---------- 7B. Prediction interval thick bar ----------
  npc_x <- function(x) (x - xlim[1]) / diff(xlim)

  p <- forestploter::add_grob(
    p,
    row = nrow(tbl), col = ci_col,
    gb_fn = grid::segmentsGrob,
    x0 = grid::unit(npc_x(pred["low"]),  "npc"),
    x1 = grid::unit(npc_x(pred["high"]), "npc"),
    y0 = grid::unit(0.5, "npc"),
    y1 = grid::unit(0.5, "npc"),
    gp = grid::gpar(lwd = 6)
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
