#' Forest plot for trivariate meta-analysis summaries
#'
#' @description
#' Creates a forest plot for net benefit, relative utility, sensitivity, or
#' specificity using posterior samples from a Bayesian trivariate meta-analysis
#' and the corresponding study-level data.
#'
#' @param samples A `coda::mcmc.list` of posterior samples, typically
#'   `fit$samples` from [MA_NB_tri()].
#' @param data A data frame containing study-level data and labels.
#' @param label_cols Character vector of columns in `data` to display in the
#'   plot table.
#' @param study_label_col Optional column in `data` used for study labels and
#'   summary row labels. If `NULL`, the first entry of `label_cols` is used.
#' @param metric Metric to plot. One of `"NB"`, `"RU"`, `"sens"`, or `"spec"`.
#' @param center Summary measure used for pooled point estimates. One of
#'   `"Median"` or `"Mean"`.
#' @param t Threshold probability. Required when observed-study net benefit or
#'   analytical net benefit confidence intervals are used.
#' @param xlim Optional x-axis limits.
#' @param xticks Optional tick marks for the x-axis.
#' @param plot_col_width Width of the blank plotting column used internally by
#'   `forestploter`.
#' @param reported_low_col Optional column in `data` containing reported lower
#'   interval bounds.
#' @param reported_high_col Optional column in `data` containing reported upper
#'   interval bounds.
#' @param reported_est_col Optional column in `data` containing reported point
#'   estimates.
#' @param interval_fallback Optional fallback interval method. If `NULL`, a
#'   default is chosen based on `metric`.
#' @param mark_imputed Logical; whether to mark study rows where displayed point
#'   estimates or intervals were filled in rather than fully taken from
#'   reported study values.
#' @param file_png Optional file path for saving the plot as PNG.
#' @param file_pdf Optional file path for saving the plot as PDF.
#'
#' @returns
#' A forest plot object produced by the `forestploter` package. The plot is also
#' drawn as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- MA_NB_tri(data, tp = tp, tn = tn,
#'                  n_event = n_event, n_nonevent = n_nonevent,
#'                  t = 0.1, prior_type = "weak", seed = 123)
#' plot_forest(samples = fit$samples, data = data, metric = "NB",
#'             label_cols = c("Publication", "Country", "N", "Prev"),
#'             t = 0.1)
#' }
plot_forest <- function(
    samples,
    data,
    label_cols = NULL,
    study_label_col = NULL,
    metric = c("NB", "RU", "sens", "spec"),
    center = c("Median", "Mean"),
    t = NULL,
    xlim = NULL,
    xticks = NULL,
    plot_col_width = 50,
    reported_low_col = NULL,
    reported_high_col = NULL,
    reported_est_col = NULL,
    interval_fallback = NULL,
    mark_imputed = TRUE,
    file_png = NULL,
    file_pdf = NULL
) {
  xlim_user <- !is.null(xlim)
  metric <- match.arg(metric)

  # ---- samples -> plot summary ----
  if (!inherits(samples, "mcmc.list")) {
    stop("`samples` must be a coda::mcmc.list (e.g., fit$samples).", call. = FALSE)
  }
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    stop("`data` must be provided as a data.frame for study labels.", call. = FALSE)
  }

  # Decide which families to summarize based on the metric we are plotting
  targets <- if (metric == "NB") c("NB", "probuseful") else metric
  targets_per_study <- metric

  # Build the plot-ready summary object internally
  sum <- .summarize_for_forestplot(
    samples = samples,
    data = data,
    label_cols = label_cols,
    targets = targets,
    targets_per_study = targets_per_study,
    return_known = FALSE
  )
  sum_m  <- sum[[metric]]

  if (is.null(sum_m)) {
    stop("Metric ", metric, " not found in plot summary.", call. = FALSE)
  }
  if (is.null(sum_m$per_study)) {
    stop("No per-study results found for metric ", metric, ".", call. = FALSE)
  }
  if (metric == "NB" && is.null(t) &&
      (is.null(reported_est_col) || is.null(reported_low_col) || is.null(reported_high_col) ||
       identical(interval_fallback, "analytic"))) {
    message("`t` is required for observed-study NB calculations or analytical NB CI fallback.")
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

  # point estimate labels
  point_col_label <- switch(
    metric,
    NB = "NB",
    RU = "RU",
    sens = "Sens",
    spec = "Spec"
  )


  # keep only label columns that exist in the data user provided
  label_cols <- intersect(label_cols, names(data))
  if (length(label_cols) == 0) {
    stop("None of `label_cols` exist in `data`.", call. = FALSE)
  }

  # which column should receive "Pooled estimate", etc.
  if (is.null(study_label_col)) {
    study_label_col <- label_cols[1]
  } else {
    if (!study_label_col %in% names(data)) {
      stop("`study_label_col` = '", study_label_col, "' not found in `data`.", call. = FALSE)
    }
    # ensure it's included in the table
    if (!study_label_col %in% label_cols) label_cols <- c(study_label_col, label_cols)
  }


  # ---------- 1. Per-study table (sorted by point estimate) ----------
  per <- sum_m$per_study

  # ---------- 1B. Build displayed per-study values ----------
  if (is.null(interval_fallback)) {
    interval_fallback <- switch(
      metric,
      sens = "frequentist",
      spec = "frequentist",
      NB   = "analytic",
      RU   = "model",
      "none"
    )
  }

  disp <- build_per_study_display(
    metric = metric,
    per = per,
    data = data,
    center = center,
    t = t,
    reported_est_col = reported_est_col,
    reported_low_col = reported_low_col,
    reported_high_col = reported_high_col,
    interval_fallback = interval_fallback
  )

  display_est  <- disp$display_est
  display_low  <- disp$display_low
  display_high <- disp$display_high

  # interval labels
  interval_col_label <- infer_interval_label(
    metric = metric,
    interval_source = disp$interval_source
  )

  emit_forest_messages(
    metric = metric,
    point_source = disp$point_source,
    interval_source = disp$interval_source,
    mark_imputed = mark_imputed
  )

  # sort by what is actually displayed
  ord <- order(display_est)
  per <- per[ord, , drop = FALSE]
  disp <- disp[ord, , drop = FALSE]
  display_est  <- display_est[ord]
  display_low  <- display_low[ord]
  display_high <- display_high[ord]


  # Mark rows where displayed values were filled in rather than fully taken from reported study values
  point_flag <- ifelse(
    mark_imputed & disp$point_source != "reported",
    "\u2020",   # dagger
    ""
  )

  interval_flag <- ifelse(
    mark_imputed & disp$interval_source != "reported",
    "*",
    ""
  )

  df <- per %>%
    dplyr::mutate(
      .display_est = display_est,
      .display_low = display_low,
      .display_high = display_high,
      ` ` = strrep(" ", plot_col_width),
      value_est = sprintf("%.2f%s", .display_est, point_flag),
      value_int = sprintf("[%.2f; %.2f]%s", .display_low, .display_high, interval_flag)
    ) %>%
    dplyr::select(dplyr::all_of(label_cols), ` `, value_est, value_int)


  # IMPORTANT: make label columns character so bind_rows works with blank rows
  df <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(label_cols), as.character))


  names(df)[ncol(df) - 1] <- point_col_label
  names(df)[ncol(df)] <- interval_col_label

  # build a "blank row" template
  make_row <- function(title, value_est = "", value_int = "") {
    row <- as.list(rep("", length(label_cols)))
    names(row) <- label_cols
    row[[study_label_col]] <- title

    row <- tibble::as_tibble(row) %>%
      dplyr::mutate(
        ` ` = strrep(" ", plot_col_width),
        value_est = value_est,
        value_int = value_int
      )

    names(row)[ncol(row) - 1] <- point_col_label
    names(row)[ncol(row)] <- interval_col_label
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
    value_est = sprintf("%.2f", pooled[center]),
    value_int = sprintf("[%.2f; %.2f]", pooled["Low"], pooled["High"])
  )

  # ---------- 3. Prediction interval row ----------
  pred_row <- make_row(
    "Prediction interval",
    value_est = "",
    value_int = sprintf("[%.2f; %.2f]", pred["Low"], pred["High"])
  )

  # ---------- 3B. P(useful) row (text only) ----------
  pu_row <- if (show_pu) {
    make_row("P(useful)", value_est = sprintf("%.2f", prob_useful), value_int = "")
  } else NULL

  # ---------- 4. Combine table ----------
  tbl <- dplyr::bind_rows(df, pooled_row, pred_row, pu_row) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

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
    grDevices::png(file_png, width = wh[1], height = wh[2],
        units = "in", res = 300, bg = "white")
    plot(p)
    grDevices::dev.off()
  }

  if (!is.null(file_pdf)) {
    grDevices::cairo_pdf(file_pdf, width = wh[1], height = wh[2], bg = "white")
    plot(p)
    grDevices::dev.off()
  }

  invisible(p)
}

