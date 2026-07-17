# MetaNB

MetaNB is an R package for Bayesian trivariate random-effects meta-analysis of prevalence, sensitivity, and specificity to indirectly synthesize net benefit across validation settings. The package also provides value-of-information (VOI) measures to quantify uncertainty around implementation decisions.

## Vignette

A worked introduction to the package is available here:

https://zhipeiwang.github.io/MetaNB/Introduction_to_MetaNB.html

The vignette demonstrates a typical workflow, including:

* fitting the Bayesian trivariate meta-analysis model;
* assessing convergence;
* summarizing posterior draws;
* visualizing results using forest plots; and
* obtaining value-of-information metrics.

## Installation

MetaNB uses JAGS through the `rjags` package. JAGS must therefore be installed on the system before MetaNB can be used. You can download JAGS from https://sourceforge.net/projects/mcmc-jags/files/

After installing JAGS, install MetaNB with:

```r
remotes::install_github("zhipeiwang/MetaNB")
```

## Main functionality
The key outputs currently provided by MetaNB include:

* Posterior distributions obtained from MCMC sampling, including (but not limited to) per-setting, pooled, and predictive net benefit for model, treat-all and treat-none strategies;
* a posterior probability that a prediction model is clinically useful in a new setting;
* forest plots for net benefit, relative utility, sensitivity, and specificity;
* value-of-information metrics, including:
  * expected value of perfect information (EVPI) for the overall implementation decision;
  * cluster-level EVPI, where each unobserved setting is allowed to choose its own optimal strategy;
  * expected value of perfect partial information (EVPPI) when the decision is to use the optimal strategy for a cluster with a given prevalence; and
  * EVPI in a specific observed cluster.

## Development status

MetaNB is currently under active development. Function names, arguments, and outputs may change as functionality is expanded and refined.

Feedback, bug reports, and suggestions are welcome.
