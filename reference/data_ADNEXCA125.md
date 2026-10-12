# Example dataset for ovarian cancer prediction model validation

Study-level data used to illustrate the trivariate meta-analysis of net
benefit and value-of-information analysis.

## Usage

``` r
data_ADNEXCA125
```

## Format

A data frame with variables:

- Publication:

  Study identifier.

- Country:

  Country of study.

- N:

  Study sample size.

- Prev:

  Study prevalence, stored as a character percentage.

- n_nonevent:

  Number of non-events.

- n_event:

  Number of events.

- tp:

  Number of true positives.

- tn:

  Number of true negatives.

- sens_point:

  Reported study sensitivity point estimate, if available.

- sens_ci_low:

  Reported lower confidence bound for sensitivity, if available.

- sens_ci_high:

  Reported upper confidence bound for sensitivity, if available.

- spec_point:

  Reported study specificity point estimate, if available.

- spec_ci_low:

  Reported lower confidence bound for specificity, if available.

- spec_ci_high:

  Reported upper confidence bound for specificity, if available.

## Source

<https://doi.org/10.1136/bmjmed-2023-000817>
