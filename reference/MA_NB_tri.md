# Bayesian trivariate meta-analysis of net benefit

Fits a Bayesian trivariate random-effects meta-analysis model jointly
modelling sensitivity, specificity, and prevalence across studies using
study-level count data. Returns posterior MCMC samples for Net Benefit
(NB), Relative Utility (RU), sensitivity, specificity and related
clinical utility metrics at a user-specified decision threshold.

When `compute_EVPI = TRUE`, also computes Value of Information (VOI)
metrics including the Expected Value of Perfect Information (EVPI) and
optionally the Expected Value of Partial Perfect Information (EVPPI) for
prevalence, with optional iterative resampling until Monte Carlo
precision criteria are met.

## Usage

``` r
MA_NB_tri(
  data,
  tp,
  tn,
  n_event,
  n_nonevent,
  prior_type = c("weak", "wishart"),
  t = NULL,
  prev_known = 0.5,
  return_known = FALSE,
  compute_EVPI = FALSE,
  J = 1000,
  inits = NULL,
  n.chains = 2,
  n.adapt = 1000,
  burnin = 3000,
  iter = 1000,
  thin = 1,
  seed = NULL,
  rng_name = "base::Wichmann-Hill",
  diag_vars = c("pooledsens", "pooledprev", "pooledspec", "pooledNB", "pooledNB_TA",
    "pooledRU"),
  return_vars = c("NB", "probuseful"),
  weak_priors = list(),
  wishart_priors = list(),
  auto_resample = TRUE,
  max_rounds = 10,
  sample_by = 2000,
  max_draws = 2e+05,
  sigma_min = 2,
  ess_min = 400,
  compute_evppi_prev = TRUE,
  evppi_method = "sal",
  center_rows = NULL,
  center_label_cols = NULL
)
```

## Arguments

- data:

  A data frame containing one row per study.

- tp:

  Column of true positives (unquoted name or string).

- tn:

  Column of true negatives (unquoted name or string).

- n_event:

  Column of total events per study.

- n_nonevent:

  Column of total non-events per study.

- prior_type:

  Character string specifying the prior structure. `"weak"` uses weakly
  informative priors. `"wishart"` uses an inverse-Wishart prior. Default
  `"weak"`.

- t:

  Numeric scalar in (0, 1). The decision threshold (i.e. harm-to-benefit
  ratio). Required; no default.

- prev_known:

  Numeric scalar in (0, 1). An assumed known prevalence used to compute
  additional NB/RU estimates at a fixed prevalence. Only returned if
  `return_known = TRUE`. Default `0.5`.

- return_known:

  Logical. If `TRUE`, also returns posteriors evaluated at `prev_known`
  (e.g. `pooledNB_known`, `NBnew_known`, `pooledNB_TA_known`,
  `NBnew_TA_known`). Default `FALSE`.

- compute_EVPI:

  Logical. If `TRUE`, enables VOI computation. Forces the necessary JAGS
  nodes into the monitored set and activates the iterative sampling loop
  if `auto_resample = TRUE`. Default `FALSE`.

- J:

  Positive integer. Number of Monte Carlo samples used for the
  within-JAGS population EVPI approximation. Only used when
  `compute_EVPI = TRUE`. Default `1000`.

- inits:

  Optional list of length `n.chains`, where each element is a named list
  of initial values for one JAGS chain. If `NULL` and `seed` is
  provided, RNG seeds are set automatically per chain.

- n.chains:

  Positive integer. Number of MCMC chains. Default `2`.

- n.adapt:

  Positive integer. Number of JAGS adaptation iterations. Default
  `1000`.

- burnin:

  Positive integer. Number of burn-in iterations (discarded). Default
  `3000`.

- iter:

  Positive integer. Number of posterior sampling iterations per chain in
  the first sampling round. Default `1000`.

- thin:

  Positive integer. Thinning interval for MCMC chains. Default `1`.

- seed:

  Optional integer. Random seed for reproducibility, passed to set each
  JAGS chain's RNG. Ignored if `inits` already contains `.RNG.seed`.

- rng_name:

  Character string. JAGS RNG type. One of `"base::Wichmann-Hill"`
  (default), `"base::Marsaglia-Multicarry"`, `"base::Super-Duper"`, or
  `"base::Mersenne-Twister"`.

- diag_vars:

  Character vector of JAGS node names to monitor for convergence
  diagnostics (trace plots). Set to `NULL` or `character(0)` to suppress
  trace plots. Default:
  `c("pooledsens", "pooledprev", "pooledspec", "pooledNB", "pooledNB_TA", "pooledRU")`.

- return_vars:

  Character vector of variable families or individual JAGS node names to
  monitor and return. Family names (`"NB"`, `"RU"`, `"sens"`, `"spec"`)
  are automatically expanded to include per-study, pooled, and
  predictive nodes. When `compute_EVPI = TRUE`, VOI-required nodes are
  added automatically regardless of this argument. Default:
  `c("NB", "probuseful")`. Advanced users may request additional
  parameters. Please refer to model_text_weak.R and model_text_wishart.R
  on GitHub for the JAGS code. Common options include:

  - `NB`: Net benefit related parameters, including per-study (`NB[1]`
    to `NB[N]`, where N is the total number of studies/centers included
    in the meta-analysis), pooled (`pooledNB`, `pooledNB_TA`) and
    predictive (`NBnew`, `NBnew_TA`) posteriors, for both the model and
    treat all strategies.

  - `RU`: Relative utility related parameters, including per-study
    (`RU[1]` to `RU[N]`, where N is the total number of studies/centers
    included in the meta-analysis), pooled (`pooledRU`) and predictive
    (`RUnew`) posteriors.

  - `probuseful`: Posterior probability that the model is clinically
    useful.

  - `sens`: Sensitivity related parameters, including per-study
    (`sens[1]` to `sens[N]`, where N is the total number of
    studies/centers included in the meta-analysis), pooled
    (`pooledsens`) and predictive (`sensnew`) posteriors.

  - `spec`: Specificity related parameters, including per-study
    (`spec[1]` to `spec[N]`, where N is the total number of
    studies/centers included in the meta-analysis), pooled
    (`pooledspec`) and predictive (`specnew`) posteriors.

- weak_priors:

  Named list of scalar overrides for the weak prior hyperparameters.
  Valid names: `mu_etap`, `tau_etap`, `mu_lambdasens0`,
  `tau_lambdasens0`, `mu_lambdaspec0`, `tau_lambdaspec0`, `mu_zss`,
  `tau_zss`, `mu_zsp`, `tau_zsp`, `a_csp`, `b_csp`, `tau_varprev`,
  `tau_varsens`, `tau_varspec`. Unknown names are ignored with a
  warning. Only used when `prior_type = "weak"`. Default values are
  based on Wynants et al. (2018).

- wishart_priors:

  Named list of overrides for the Wishart prior hyperparameters. Valid
  names: `mn` (length-3 numeric vector), `prec` (3x3 matrix), `R` (3x3
  matrix), `df` (numeric scalar). Unknown names are ignored with a
  warning. Only used when `prior_type = "wishart"`. Default values are
  based on Wynants et al. (2018).

- auto_resample:

  Logical. If `TRUE` and `compute_EVPI = TRUE`, iteratively draws
  additional samples of size `sample_by` until Monte Carlo precision
  criteria (`sigma_min`, `ess_min`) are met or stopping conditions are
  reached. Default `TRUE`.

- max_rounds:

  Positive integer. Maximum number of sampling rounds when
  `auto_resample = TRUE`. Default `10`.

- sample_by:

  Positive integer. Number of additional iterations per chain drawn in
  each resampling round after the first. Default `2000`.

- max_draws:

  Positive integer. Maximum total number of posterior draws (across all
  chains) before stopping resampling. Default `200000`.

- sigma_min:

  Numeric. Minimum acceptable signal-to-noise ratio for the difference
  between competing strategies; resampling continues while sigma \<
  `sigma_min`. Default `2`.

- ess_min:

  Numeric. Minimum acceptable effective sample size for key pooled nodes
  (`pooledNB`, `pooledNB_TA`, `probuseful`); resampling continues while
  any ESS \< `ess_min`. Corresponds to a relative Monte Carlo standard
  error of at most 5%. Default `400`.

- compute_evppi_prev:

  Logical. If `TRUE` and `compute_EVPI = TRUE`, computes the EVPPI for
  prevalence using the voi package. Requires voi to be installed.
  Default `TRUE`.

- evppi_method:

  Character string. Method passed to
  [`voi::evppi()`](https://chjackson.github.io/voi/reference/evppi.html)
  for EVPPI estimation. Default `"sal"`.

- center_rows:

  Optional integer vector. Row indices of studies in `data` for which
  center-specific EVPI should be computed. Only active when
  `compute_EVPI = TRUE`. These rows determine which per-study NB nodes
  are monitored and which centers are included in the resampling
  stopping rule. Default `NULL`.

- center_label_cols:

  Optional character vector of column names in `data` used to construct
  human-readable center labels in the output. If `NULL` or none of the
  specified columns exist in `data`, the first column of `data` is used.
  Default `NULL`.

## Value

A named list with the following elements:

- `samples`: `mcmc.list` object containing posterior draws for all
  monitored nodes. ` priors_used`: A tibble summarizing the prior
  hyperparameter values used after merging defaults with any user
  overrides.

- `meta`: A named list of run metadata including threshold (`returned`),
  MCMC settings, actual draws per chain (`saved_per_chain`), total draws
  (`total_saved`), and per-chain RNG information. When
  `compute_EVPI = TRUE`, also includes `rounds`, `stop_reason`.

- `voi_metrics`: A tibble of VOI metrics including
  `EVPI_cluster_perfectinfo`, `EVPI_population_perfectinfo`, and
  `EVPPI_cluster_perfectprevalenceinfo`. `NULL` when
  `compute_EVPI = FALSE`.

- `voi_diagnostics`: A tibble of Monte Carlo precision diagnostics
  including sigma values, ESS, and `needs_more_sampling`. `NULL` when
  `compute_EVPI = FALSE`.

- `voi_center_metrics`: A tibble of center-specific EVPI estimates and
  diagnostics for each row in `center_rows`. `NULL` when `center_rows`
  is `NULL` or `compute_EVPI = FALSE`.

- `voi_center_meta`: A data frame of human-readable center labels
  corresponding to `center_rows`. `NULL` when `center_rows` is `NULL` or
  `compute_EVPI = FALSE`.

## References

Wynants L, Riley R, Timmerman D, Van Calster B. Random-effects
meta-analysis of the clinical utility of tests and prediction models.
*Stat Med* 2018;37(12):2034–52.
[doi:10.1002/sim.7653](https://doi.org/10.1002/sim.7653)

Jackson C, Heath A (2024). *voi: Value of Information Analysis in R*
[doi:10.32614/CRAN.package.voi](https://doi.org/10.32614/CRAN.package.voi)

## See also

[`summarize_tri_ma()`](https://zhipeiwang.github.io/MetaNB/reference/summarize_tri_ma.md),
[`plot_forest()`](https://zhipeiwang.github.io/MetaNB/reference/plot_forest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Tri-MA only
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
summary(fit$samples)

# With VOI and center-specific EVPI
fit_voi <- MA_NB_tri(
  data          = data,
  tp            = tp,
  tn            = tn,
  n_event       = n_event,
  n_nonevent    = n_nonevent,
  t             = 0.1,
  prior_type    = "weak",
  compute_EVPI  = TRUE,
  center_rows   = 1:10,
  center_label_cols = c("Publication", "Country"),
  seed          = 123
)
fit_voi$voi_metrics
fit_voi$voi_center_metrics
} # }
```
