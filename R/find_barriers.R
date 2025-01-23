#' Calculate Decision Boundaries for Continuous Sweep given pdelta
#'
#' @param pdelta Minimum difference between true and false positive rate.
#' @param tr List containing training parameters:
#' \itemize{
#'   \item mup - Mean of positive class distribution
#'   \item omegap - Variance of positive class distribution
#'   \item mun - Mean of negative class distribution  
#'   \item omegan - Variance of negative class distribution
#'   \item ntrain - Number of training instances
#'   \item prevtrain - Prevalence of positive class in training data
#' }
#' 
#' @return Numeric vector of length 2 containing lower and upper decision barriers
#' @export
#' @importFrom rootSolve uniroot.all
#' @importFrom stats pnorm qnorm
#'
#' @examples
#' tr <- list(mup = 1, omegap = 1, mun = 0, omegan = 0.5, ntrain = 1000, prevtrain = 0.7)
#' find_barriers(pdelta = 0.25, tr = tr)
#' 

find_barriers <- function(pdelta, tr) {
  # Functions for the positive (Fp) and negative (Fn) class. 
  Fp <- function(x) {
    stats::pnorm(x, tr$mup, sqrt(tr$omegap), lower.tail = FALSE)
  }
  Fn <- function(x) {
    stats::pnorm(x, tr$mun, sqrt(tr$omegan), lower.tail = FALSE)
  }
  
  # Target function: Difference between rates and pdelta
  Fdelta <- function(x) {
    Fp(x) - Fn(x) - pdelta
  }
  
  # Calculate search interval with safety margin to ensure root finding
  err_fct <- 0.1  # 10% margin to extend search boundaries
  lft <- tr$mun + (sqrt(tr$omegan) * (stats::qnorm(pdelta) - err_fct))
  rgt <- tr$mup + (sqrt(tr$omegap)* (stats::qnorm(1 - pdelta) + err_fct)) 

  # Find roots where Fdelta crosses zero (decision boundaries)
  out <- 
    rootSolve::uniroot.all(f = function(x) Fdelta(x), interval = c(lft, rgt)) |>
    unique()

  # Warn if unexpected number of barriers found
  if(length(out) != 2) {warning("Length barriers is not 2 - pdelta too low/high?")}
  
  return(out)
}
