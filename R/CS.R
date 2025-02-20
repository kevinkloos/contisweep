#' Continuous Sweep Prevalence Estimation
#'
#' Estimates prevalence values using continuous sweep based on:
#' \itemize{
#' \item Training parameters from positive and negative distributions
#' \item Test scores from target population
#' \item Decision boundaries (barriers) for the area to integrate
#' }
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
#' @param barriers Two numeric values indicating the left and right point of the area
#' to be integrated OR a character "basic" or "advanced" to use optimized barriers.
#' @param conf_alpha Either FALSE or a numeric value between 0 and 1 indicating the
#' alpha parameter for the confidence interval.
#' @param print_output Boolean indicating whether to print nicely formatted output.
#' @return A prevalence estimate or a list containing:
#' \itemize{
#'   \item prevalence - Prevalence estimate
#'   \item variance - Variance estimate
#'   \item conf_interval - Lower and upper bounds of the confidence interval
#'   \item conf_alpha - Alpha parameter for the confidence interval
#'   \item pdelta - Optimal pdelta estimate
#' }
#' @examples
#' # Example training parameters
#' tr <- list(mup = 1, omegap = 1, mun = 0, omegan = 0.5, ntrain = 1000, prevtrain = 0.7)
#' # Generate test scores
#' set.seed(123)
#' te <- c(
#'   rnorm(500, mean = tr$mup, sd = sqrt(tr$omegap)),
#'   rnorm(500, mean = tr$mun, sd = sqrt(tr$omegan))
#' )
#'
#' # Calculate prevalence estimate
#' CS(tr, te, barriers = c(-1, 1))
#' CS(tr, te, barriers = "advanced", conf_alpha = 0.95)
#'
#' @references
#' Method implemented in Advanced Variance Estimation for Continuous Sweep.
#'
#' @export
#' @importFrom purrr map_dbl pmap
#' @importFrom dplyr near
#' @importFrom utils head tail
#' @importFrom stats integrate
#' @importFrom purrr pmap_dbl
CS <- function(tr, te, barriers = "advanced", conf_alpha = FALSE, print_output = TRUE) {
  .calculate_one_area <- function(cac, thetal, thetar, mup, sdp, mun, sdn) {
    Fp <- function(q) stats::pnorm(q, mup, sdp, lower.tail = FALSE)
    Fn <- function(q) stats::pnorm(q, mun, sdn, lower.tail = FALSE)

    if (dplyr::near(thetal, thetar)) {
      return(0)
    }
    if (dplyr::near(Fp(thetal), Fn(thetal))) {
      return(Inf)
    }
    if (dplyr::near(Fp(thetar), Fn(thetar))) {
      return(Inf)
    }
    if (sign(Fp(thetal) - Fn(thetal)) != sign(Fp(thetar) - Fn(thetar))) {
      return(Inf)
    }

    return(stats::integrate(function(x) (cac - Fn(x)) / (Fp(x) - Fn(x)), thetal, thetar)$value)
  }
  n_obs <- length(te)
  bar_char <- is.character(barriers)

  if (!isFALSE(conf_alpha) && (conf_alpha <= 0 || conf_alpha >= 1)) {
    stop("conf_alpha must be FALSE or a numeric value between 0 and 1.")
  }
  if (bar_char) {
    opt_var <- optimize_variance(tr, te, var_type = barriers)
    barriers <- opt_var$barriers
  }

  stopifnot(length(barriers) == 2)
  stopifnot(all(!is.na(barriers)))
  stopifnot(is.numeric(barriers))

  thetal <- barriers[1]
  thetar <- barriers[2]

  sorted_te <- sort(te)
  all_theta_l <- utils::head(sorted_te, -1)
  all_theta_r <- utils::tail(sorted_te, -1)
  all_cac <- seq(1 - 1 / n_obs, 1 / n_obs, length.out = n_obs - 1)

  bool_crit <- all_theta_r >= thetal & all_theta_l <= thetar
  n_trunc <- sum(bool_crit)
  if (n_trunc == 0) {
    warning("No values of theta_l and theta_r between the barriers.") # or some other warnings, maybe even an error?
    return(NA_real_)
  }

  # Subset the values of all_theta_l and all_theta_r between theta_l and theta_r
  # Replace the values slightly below and sightly above barriers with theta_l and theta_r
  all_theta_l <-
    subset(all_theta_l, bool_crit) |>
    replace(1, thetal)
  all_theta_r <-
    subset(all_theta_r, bool_crit) |>
    replace(n_trunc, thetar)
  all_cac <- subset(all_cac, bool_crit)

  all_areas <-
    purrr::pmap_dbl(
      list(all_cac, all_theta_l, all_theta_r),
      ~ .calculate_one_area(..1, ..2, ..3, tr$mup, sqrt(tr$omegap), tr$mun, sqrt(tr$omegan))
    )

  area_sum <- sum(all_areas)
  out <- area_sum / (thetar - thetal)
  if (out > 1 | out < 0) {
    warning("Estimated prevalence outside 0-1.")
  }

  if (conf_alpha) {
    vari <- find_variance(barriers, tr, te, out, "advanced")
    out <- list(
      prevalence = out,
      variance = vari,
      conf_interval = c(
        out + qnorm((1 - conf_alpha) / 2) * sqrt(vari),
        out + qnorm(1 - (1 - conf_alpha) / 2) * sqrt(vari)
      ),
      conf_alpha = conf_alpha
    )
    if (bar_char) {
      out$pdelta <- opt_var$pdelta
      if (print_output) {
        cat(paste("\U03B1-CS = ", round(out$prevalence, 3), "; ",
          round(100 * out$conf_alpha, 3),
          "% CI:[", round(out$conf_interval[1], 3), "-",
          round(out$conf_interval[2], 3), "]", ". \n\U03c9 = ",
          signif(out$variance, 3), "; p\U0394 = ", signif(out$pdelta, 3),
          ". \n",
          sep = ""
        ))
      }
    } else {
      if (print_output) {
        cat(paste("\U03B1-CS = ", round(out$prevalence, 3), "; ",
          round(100 * out$conf_alpha, 3),
          "% CI:[", round(out$conf_interval[1], 3), "-",
          round(out$conf_interval[2], 3), "]", ". \n\U03c9 = ",
          signif(out$variance, 3), ". \n",
          sep = ""
        ))
      }
    }
  } else {
    if (print_output) {
      cat(paste("\U03B1-CS = ", round(out, 3), ".", sep = "."))
    }
  }

  if (print_output) {
    invisible(out)
  } else {
    return(out)
  }
}
