#' Optimize the variance of Continuous Sweep
#'
#' @param tr List containing training parameters:
#' \itemize{
#'   \item mup - Mean of positive class distribution
#'   \item omegap - Variance of positive class distribution
#'   \item mun - Mean of negative class distribution  
#'   \item omegan - Variance of negative class distribution
#'   \item ntrain - Number of training instances
#'   \item prevtrain - Prevalence of positive class in training data
#' }
#' @param te Numeric vector of test scores (target population)
#' @param var_type The type of variance estimate:
#' \itemise{
#'   \item basic - assumes training parameters Fp and Fn as known.
#'   \item advanced - assumes training parameters Fp and Fn as estimated.
#' } 
#' 
#' @return A list containing:
#' \itemize{
#'   \item variance - The estimated variance given optimal pdelta.
#'   \item pdelta - Optimal pdelta estimate.
#'   \item barriers - Two decision boundaries.
#'   \item ac - An initial adjusted count estimate using AC-MAX.
#'   \item var_type - The type of variance estimate (basic or advanced).
#' }
#' @export
#' @importFrom rootSolve uniroot.all
#' @importFrom pracma integral2
#' @importFrom stats pnorm dnorm integrate
#' @importFrom utils head tail
#' @importFrom purrr pmap
#'
#' @examples
#' # Example training parameters
#' tr <- list(mup = 1, omegap = 1, mun = 0, omegan = 0.5, ntrain = 1000, prevtrain = 0.7)
#' # Generate test scores
#' set.seed(123)
#' te <- c(stats::rnorm(500, mean = tr$mup, sd = sqrt(tr$omegap)), 
#'         stats::rnorm(500, mean = tr$mun, sd = sqrt(tr$omegan)))
#' # Apply the function
#' optimize_variance(tr, te, var_type = "advanced")
#' 
optimize_variance <- function(tr, te, var_type = "advanced") {

  stopifnot(var_type == "basic" || var_type ==  "advanced")

  Fp <- function(x) {
    stats::pnorm(x, tr$mup, sqrt(tr$omegap), lower.tail = FALSE)
  }
  Fn <- function(x) {
    stats::pnorm(x, tr$mun, sqrt(tr$omegan), lower.tail = FALSE)
  }
  Fdiff <- function(x) {
    Fp(x) - Fn(x)
  }

  Vbasic <- function(x) {
    find_barriers(x, tr) |> find_variance(tr, te, ac_pmax, var_type = "basic")
  }

  Vadv <- function(x) {
    find_barriers(x, tr) |> find_variance(tr, te, ac_pmax, var_type = "advanced")
  }

  pmin = 0.001
  prevtrain <- tr$prevtrain
  ntrain <- tr$ntrain
  ntest <- length(te)

  lft <- tr$mun + sqrt(tr$omegap) * stats::qnorm(pmin)
  rgt <- tr$mup + sqrt(tr$omegan) * stats::qnorm(1 - pmin)

  ppoint <- stats::optimize(Fdiff, lower = lft, upper = rgt, maximum = TRUE)
  pmax <- ppoint$objective
  thmax <- ppoint$maximum

  ccmax <- mean(te >= thmax)
  ac_pmax <- (ccmax - Fn(thmax)) / (Fp(thmax) - Fn(thmax))
  ac_pmax <- ac_pmax |> pmax(0) |> pmin(1)

  eps <- 0.001
  if(var_type == "basic") {
    pdelta_obj <- stats::optimize(Vbasic, lower = pmin, upper = pmax - eps)
  }
  else if(var_type == "advanced") {
    pdelta_obj <- stats::optimize(Vadv, lower = pmin, upper = pmax - eps)
  }

  out <- list(pdelta = pdelta_obj$minimum, barriers = find_barriers(pdelta_obj$minimum, tr),
              variance = pdelta_obj$objective, ac_max = ac_pmax, var_type = var_type)
  return(out)
}
