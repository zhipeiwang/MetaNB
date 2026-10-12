# Summarize posterior samples from a Bayesian trivariate meta-analysis

Extracts and structures posterior summaries from a fitted Bayesian
trivariate meta-analysis. For each requested metric, returns pooled
estimates (posterior mean, median, and 95% credible interval),
predictive summaries for a new center (mean, median, and 95% prediction
interval), and optionally per-study summaries and summaries at a fixed
known prevalence.

## Usage

``` r
summarize_tri_ma(
  samples,
  data = NULL,
  label_cols = NULL,
  metrics = c("NB", "probuseful"),
  per_study_metrics = c("NB"),
  return_known = FALSE,
  include_per_study = TRUE
)
```

## Arguments

- samples:

  A fitted object returned by
  [`MA_NB_tri()`](https://zhipeiwang.github.io/MetaNB/reference/MA_NB_tri.md),
  or a [`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)
  of posterior samples directly.

- data:

  Optional data frame containing study-level information. Required when
  per-study summaries are requested (i.e. when
  `include_per_study = TRUE` and a metric in `per_study_metrics` is
  requested).

- label_cols:

  Optional character vector of column names in `data` to retain in
  per-study output (e.g. study labels, country, sample size).

- metrics:

  Character vector of metrics to summarize. Family names (`"NB"`,
  `"RU"`, `"probuseful"`) and atomic JAGS node names (e.g.
  `"pooledprev"`, `"NBnew_known"`) are both accepted. Default:
  `c("NB", "probuseful")`.

- per_study_metrics:

  Character vector of metrics for which per-study summaries should be
  returned. Only used when `include_per_study = TRUE`. Default:
  `c("NB")`.

- return_known:

  Logical. If `TRUE`, also returns summaries evaluated at the fixed
  known prevalence (e.g. `pooledNB_known`, `NBnew_known`). Requires that
  `return_known = TRUE` was set when fitting with
  [`MA_NB_tri()`](https://zhipeiwang.github.io/MetaNB/reference/MA_NB_tri.md).
  Default `FALSE`.

- include_per_study:

  Logical. If `TRUE`, per-study posterior summaries are included for
  metrics listed in `per_study_metrics`. Requires `data` to be supplied.
  Default `TRUE`.

## Value

A named list with one element per requested metric. For family metrics
(`"NB"`, `"RU"`), each element contains:

- `per_study`: A data frame of per-study posterior summaries (mean,
  median, 95% CrI), merged with `label_cols` from `data`. `NULL` if
  `include_per_study = FALSE`.

- `pooled`: Named list of pooled posterior summaries (mean, median, 95%
  CrI) for the model and treat-all strategies.

- `predictive`: Named list of posterior predictive summaries (mean,
  median, 95% PI) for a new center for the model and treat-all
  strategies.

- `pooled_known`: As `pooled` but at fixed known prevalence. `NULL` if
  `return_known = FALSE`.

- `pred_known`: As `predictive` but at fixed known prevalence. `NULL` if
  `return_known = FALSE`.

  For `"probuseful"`, a list with `value` and optionally `value_known`.
  For atomic node names, a list with scalar\` containing the posterior
  summary.

## See also

[`MA_NB_tri()`](https://zhipeiwang.github.io/MetaNB/reference/MA_NB_tri.md),
[`plot_forest()`](https://zhipeiwang.github.io/MetaNB/reference/plot_forest.md)

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

# Default summaries
summary <- summarize_tri_ma(fit, data = data,
                        label_cols = c("Publication", "Country"))
summary$NB$pooled$model       # pooled NB: mean, median, 95% CrI
summary$NB$predictive$model   # predictive NB for new center
summary$probuseful$value      # P(useful)

# Include known-prevalence summaries
fit2 <- MA_NB_tri(..., prev_known = 0.5, return_known = TRUE)
summary2 <- summarize_tri_ma(fit2, data = data, return_known = TRUE)
summary2$NB$pooled_known$model
} # }
```
