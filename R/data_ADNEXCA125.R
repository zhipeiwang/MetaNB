#' Example dataset for ovarian cancer prediction model validation
#'
#' Study-level data used to illustrate the trivariate meta-analysis
#' of net benefit and value-of-information analysis.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{Publication}{Study identifier.}
#'   \item{Country}{Country of study.}
#'   \item{N}{Study sample size.}
#'   \item{Prev}{Study prevalence, stored as a character percentage.}
#'   \item{n_nonevent}{Number of non-events.}
#'   \item{n_event}{Number of events.}
#'   \item{tp}{Number of true positives.}
#'   \item{tn}{Number of true negatives.}
#'   \item{sens_point}{Reported study sensitivity point estimate, if available.}
#'   \item{sens_ci_low}{Reported lower confidence bound for sensitivity, if available.}
#'   \item{sens_ci_high}{Reported upper confidence bound for sensitivity, if available.}
#'   \item{spec_point}{Reported study specificity point estimate, if available.}
#'   \item{spec_ci_low}{Reported lower confidence bound for specificity, if available.}
#'   \item{spec_ci_high}{Reported upper confidence bound for specificity, if available.}
#' }
#' @source \url{https://doi.org/10.1136/bmjmed-2023-000817}
"data_ADNEXCA125"
