# Forest plot for trivariate meta-analysis summaries

Produces a forest plot for net benefit (NB), relative utility (RU),
sensitivity, or specificity from a fitted Bayesian trivariate
meta-analysis. Per-study point estimates and intervals are displayed
alongside the pooled posterior estimate and the 95% prediction interval
for a new center. For the forest plot for net benefit, the posterior
probability \\P(\mathrm{useful})\\ that the model is useful at a new
center is also shown.

Per-study point estimates and intervals are determined by a priority
hierarchy: reported values from `data` take precedence, followed by
observed or analytically derived values, followed by model-based
posterior summaries. Imputed values are flagged with \\\dagger\\ (point
estimates) and \\\ast\\ (intervals) when `mark_imputed = TRUE`.

## Usage

``` r
plot_forest(
  samples,
  data,
  tp = NULL,
  tn = NULL,
  n_event = NULL,
  n_nonevent = NULL,
  label_cols = NULL,
  study_label_col = NULL,
  metric = c("NB", "RU", "sens", "spec"),
  center = c("Mean", "Median"),
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
)
```

## Arguments

- samples:

  A [`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) of
  posterior samples, typically `fit$samples` from
  [`MA_NB_tri()`](https://zhipeiwang.github.io/MetaNB/reference/MA_NB_tri.md).

- data:

  A data frame containing study-level data and display labels. Must have
  one row per study in the same order as the fitted model.

- tp:

  Column of true positives (unquoted name or string).

- tn:

  Column of true negatives (unquoted name or string).

- n_event:

  Column of total events per study.

- n_nonevent:

  Column of total non-events per study.

- label_cols:

  Character vector of column names in `data` to display as table columns
  in the plot.

- study_label_col:

  Optional column name in `data` that receives the summary row labels
  ("Pooled estimate", "Prediction interval", "P(useful)"). If `NULL`,
  the first entry of `label_cols` is used.

- metric:

  Metric to plot. One of `"NB"` (net benefit), `"RU"` (relative
  utility), `"sens"` (sensitivity), or `"spec"` (specificity). Default
  `"NB"`.

- center:

  Posterior summary used as the point estimate per study and for the
  pooled row, and for sorting studies. One of `"Mean"` (default) or
  `"Median"`.

- t:

  Numeric scalar in (0, 1). The decision threshold. Required when
  `metric = "NB"` and observed or analytically derived per-study
  estimates are needed.

- xlim:

  Optional numeric vector of length 2 giving x-axis limits. If `NULL`,
  limits are derived automatically from the data with padding, and 0 is
  always included.

- xticks:

  Optional numeric vector of x-axis tick positions.

- plot_col_width:

  Integer. Width in characters of the blank column used internally by
  the package `forestploter` to render the plot panel. Increase if the
  plot appears cramped. Default `50`.

- reported_low_col:

  Optional column name in `data` containing reported lower interval
  bounds.

- reported_high_col:

  Optional column name in `data` containing reported upper interval
  bounds.

- reported_est_col:

  Optional column name in `data` containing reported per-study point
  estimates. Takes priority over observed and model-based values.

- interval_fallback:

  Optional character string specifying the fallback interval method for
  studies lacking reported intervals. If `NULL`, defaults are chosen by
  metric: `"frequentist"` for NB (Sande et al. 2020), and also
  `"frequentist"` for sensitivity and specificity (Wilson's method via
  the `madad` function in the package `mada`), and `"model"` (posterior
  CrI) for RU.

- mark_imputed:

  Logical. If `TRUE`, per-study rows where the displayed point estimate
  or interval was not taken directly from reported values are flagged
  with \\\dagger\\ and \\\ast\\, respectively. Default `TRUE`.

- file_png:

  Optional file path. If supplied, the plot is saved as a PNG at 300 dpi
  in addition to being displayed.

- file_pdf:

  Optional file path. If supplied, the plot is saved as a PDF in
  addition to being displayed.

## Value

The `forestploter` plot object, returned invisibly. The plot is also
rendered as a side effect. If `file_png` or `file_pdf` are supplied, the
plot is additionally saved to those paths.

## Details

Studies are sorted in ascending order of the displayed point estimate.
The pooled row uses a filled diamond and reflects the posterior mean or
median of the pooled estimate across the observed studies. The
prediction interval row reflects the 95% interval for a hypothetical new
center and is rendered as a horizontal bar with no point estimate.

Interval column headers are labelled `95% CrI` when all displayed
intervals are model-based, and `95% CI` otherwise.

## References

Wynants L, Riley R, Timmerman D, Van Calster B. Random-effects
meta-analysis of the clinical utility of tests and prediction models.
*Stat Med* 2018;37(12):2034–52.
[doi:10.1002/sim.7653](https://doi.org/10.1002/sim.7653)

Agresti A, Coull BA. Approximate is better than "exact" for interval
estimation of binomial proportions. *Am Stat* 1998;52(2):119–26.
[doi:10.2307/2685469](https://doi.org/10.2307/2685469)

Brown LD, Cai TT, DasGupta A. Interval estimation for a binomial
proportion. *Stat Sci* 2001;16(2):101–17.

Sande SZ, Li J, D'Agostino R, Wong TY, Cheng CY. Statistical inference
for decision curve analysis, with applications to cataract diagnosis.
*Stat Med* 2020;39(22):2980–3002.
[doi:10.1002/sim.8588](https://doi.org/10.1002/sim.8588)

Dayimu A (2026). *forestploter: Create a Flexible Forest Plot*.
[doi:10.32614/CRAN.package.forestploter](https://doi.org/10.32614/CRAN.package.forestploter)

Doebler P (2025). *mada: Meta-Analysis of Diagnostic Accuracy*.
[doi:10.32614/CRAN.package.mada](https://doi.org/10.32614/CRAN.package.mada)

## See also

[`MA_NB_tri()`](https://zhipeiwang.github.io/MetaNB/reference/MA_NB_tri.md),
[`summarize_tri_ma()`](https://zhipeiwang.github.io/MetaNB/reference/summarize_tri_ma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- MA_NB_tri(
  data       = data,
  tp         = tp,
  tn         = tn,
  n_event    = n_event,
  n_nonevent = n_nonevent,
  t          = 0.1,
  prior_type = "weak",
  seed       = 123
)

# Basic NB forest plot
plot_forest(
  samples   = fit$samples,
  data      = data,
  label_cols = c("Publication", "Country"),
  metric    = "NB",
  t         = 0.1
)

# Sensitivity forest plot saved to file
plot_forest(
  samples   = fit$samples,
  data      = data,
  tp        = tp,
  tn        = tn,
  label_cols = c("Publication", "Country"),
  metric    = "sens",
  file_png  = "forest_sens.png"
)
} # }
```
