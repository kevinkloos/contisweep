#' Estimate the variance of Continuous Sweep
#'
#' @param barriers The left and right decision boundaries
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
#' @param prev Initial estimate of the prevalence, preferably from adjusted count.
#' @param var_type The type of variance estimate:
#' \itemize{
#'   \item basic - assumes training parameters Fp and Fn as known.
#'   \item advanced - assumes training parameters Fp and Fn as estimated.
#' } 
#' 
#' @return A scalar containing the estimated variance
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
#' # Obtain decision boundaries
#' barriers <- c(-3, 3) # manual
#' barriers <- find_barriers(pdelta = 0.25, tr = tr) # using a pdelta
#' # Apply the function
#' find_variance(barriers, tr, te, prev = 0.6, var_type = "advanced")

find_variance <- function(barriers, tr, te, prev = ac_pmax, var_type = "advanced") {

  stopifnot(var_type == "basic" || var_type ==  "advanced")

  # return infinity when there are no barriers as input
  if (length(barriers) != 2) {
    return(Inf)
  }
  ntest <- length(te)

  Fp <- function(x) {
    stats::pnorm(x, tr$mup, sqrt(tr$omegap), lower.tail = FALSE)
  }
  Fn <- function(x) {
    stats::pnorm(x, tr$mun, sqrt(tr$omegan), lower.tail = FALSE)
  }
  V1 <- function(u, v, trn = tr, tes = te, ntr = tr$ntrain, a0 = tr$prevtrain) {
    Dh <- .delta_h(tes, u, v, trn$mup, trn$omegap, trn$mun, trn$omegan)
    Vmat <- .V_matrix(ntr, u, v, trn$mup, trn$omegap, trn$mun, trn$omegan, a0)
    out <- as.numeric(t(Dh) %*% Vmat %*% Dh)
    return(out)
  }

  V2 <- function(u, v, tr = tr, a0 = prev) {
    n1 <- Fp(pmax(u, v)) - Fp(u) * Fp(v)
    n2 <- Fn(pmax(u, v)) - Fn(u) * Fn(v)
    num <- a0 * n1 + (1 - a0) * n2
    den <- (Fp(u) - Fn(u)) * (Fp(v) - Fn(v))
    return(num / den)
  }

  intgl2 <- pracma::integral2(
    fun = function(x, y) {
      V2(x, y)
    }, barriers[1], barriers[2], barriers[1], barriers[2]
  )$Q

  const2 <- 1 / (ntest * (barriers[2] - barriers[1])^2)
  p2 <- const2 * intgl2
  if(var_type == "basic") {
    return(p2)
  }

  p1 <- V1(barriers[1], barriers[2])

  if(var_type == "advanced") {
    out <- p1 + p2
    return(out)
  }
}


.delta_h <- function(test_scores, thetal, thetar, mup, omegap, mun, omegan) {
  # compute integrals between two thresholds where cac is constant
  .between_two_thresholds <- function(cac, thetal, thetar, mup, omegap, mun, omegan) {
    sdp <- sqrt(omegap)
    sdn <- sqrt(omegan)

    if (dplyr::near(thetal, thetar)) {
      return(rep(0, 7))
    }
    Fp <- function(q) {
      stats::pnorm(q, mup, sdp, lower.tail = FALSE)
    }
    Fn <- function(q) {
      stats::pnorm(q, mun, sdn, lower.tail = FALSE)
    }
    fp <- function(x) {
      stats::dnorm(x, mup, sdp)
    }
    fn <- function(x) {
      stats::dnorm(x, mun, sdn)
    }

    dmup <- stats::integrate(f = function(x) {
      (cac - Fn(x)) / (Fp(x) - Fn(x))^2 * fp(x)
    }, lower = thetal, upper = thetar)$value

    domegap <- stats::integrate(f = function(x) {
      (cac - Fn(x)) / (Fp(x) - Fn(x)) * fp(x) * (x - mup) / (2 * sdp^2)
    }, lower = thetal, upper = thetar)$value

    dmun1 <- stats::integrate(f = function(x) {
      fn(x) / (Fp(x) - Fn(x))
    }, lower = thetal, upper = thetar)$value

    dmun2 <- stats::integrate(f = function(x) {
      fn(x) * (cac - Fn(x)) / (Fp(x) - Fn(x))^2
    }, lower = thetal, upper = thetar)$value

    domegan1 <- stats::integrate(f = function(x) {
      ((x - mun) * fn(x)) / (2 * sdn^2 * (Fp(x) - Fn(x)))
    }, lower = thetal, upper = thetar)$value

    domegan2 <- stats::integrate(f = function(x) {
      ((cac - Fn(x)) * (x - mun) * fn(x)) / (2 * sdn^2 * (Fp(x) - Fn(x))^2)
    }, lower = thetal, upper = thetar)$value

    dcs <- stats::integrate(f = function(x) {
      (cac - Fn(x)) / (Fp(x) - Fn(x))
    }, lower = thetal, upper = thetar)$value

    # summed[1] till summed[7]
    out <- c(dmup, domegap, dmun1, dmun2, domegan1, domegan2, dcs)

    return(out)
  }

  # compute the number of test observations
  n_obs <- length(test_scores)

  # make 3 vectors containing the left bound, right bound, and cac-value
  all_theta_l <- test_scores |>
    sort() |>
    utils::head(-1)
  all_theta_r <- test_scores |>
    sort() |>
    utils::tail(-1)
  all_cac <- seq(1 - 1 / n_obs, 1 / n_obs, -1 / n_obs)

  # check which values are between thetal and thetar
  bool_crit <- all_theta_r >= thetal & all_theta_l <= thetar

  # compute the "number of intervals"
  n_trunc <- sum(bool_crit)

  # replace the values slightly below and sightly above barriers
  # with correct values
  all_theta_l <-
    subset(all_theta_l, bool_crit) |>
    replace(1, thetal)
  all_theta_r <-
    subset(all_theta_r, bool_crit) |>
    replace(n_trunc, thetar)
  all_cac <- subset(all_cac, bool_crit)

  # compute all intervals between thetal and thetar where cac is constant
  all_areas <-
    purrr::pmap(
      .l = list(all_cac, all_theta_l, all_theta_r),
      .f = ~ .between_two_thresholds(
        ..1, ..2, ..3,
        mup, omegap, mun, omegan
      )
    )

  # and sum every part
  summed <- do.call("rbind", all_areas) |> colSums()

  sdp <- sqrt(omegap)
  sdn <- sqrt(omegan)

  Fp <- function(q) {
    stats::pnorm(q, mup, sdp, lower.tail = FALSE)
  }
  Fn <- function(q) {
    stats::pnorm(q, mun, sdn, lower.tail = FALSE)
  }

  # compute the adjusted count at thetal and thetar
  ac_thetal <- (all_cac[1] - Fn(thetal)) / (Fp(thetal) - Fn(thetal))
  ac_thetar <- (all_cac[n_trunc] - Fn(thetar)) / (Fp(thetar) - Fn(thetar))

  # compute every partial derivative using the integrals and other constants
  out <- c(
    -1 / (thetar - thetal) * summed[1],
    -1 / (thetar - thetal) * summed[2],
    -1 / (thetar - thetal) * summed[3] + 1 / (thetar - thetal) * summed[4],
    -1 / (thetar - thetal) * summed[5] + 1 / (thetar - thetal) * summed[6],
    summed[7] / (thetar - thetal)^2 - 1 / (thetar - thetal) * ac_thetal,
    -summed[7] / (thetar - thetal)^2 + 1 / (thetar - thetal) * ac_thetar
  )
  return(out)
}

.V_matrix <- function(ntrain, thetal, thetar, mup, omegap, mun, omegan, prev = prev_train) {

  sdp <- sqrt(omegap)
  sdn <- sqrt(omegan)
  
  fp <- function(x) {
    stats::dnorm(x, mup, sdp)
  }
  fn <- function(x) {
    stats::dnorm(x, mun, sdn)
  }

  Vfp <- function(th) {
    fp(th)^2 * (omegap + 0.5 * (th - mup)^2) / (ntrain * prev)
  }
  Vfn <- function(th) {
    fn(th)^2 * (omegan + 0.5 * (th - mun)^2) / (ntrain * (1 - prev))
  }


  Vaa <-
    matrix(
      c(omegap / (ntrain * prev), 0, 
        0, 2 * omegap^2 / (ntrain * prev)),
      ncol = 2, nrow = 2
    )

  Vbb <-
    matrix(
      c(omegan / (ntrain * (1 - prev)), 0, 
        0, 2 * omegan^2 / (ntrain * (1 - prev))),
      ncol = 2, nrow = 2
    )

  Vab <- matrix(rep(0, 4), ncol = 2, nrow = 2)
  Vba <- matrix(rep(0, 4), ncol = 2, nrow = 2)

  Eas <- matrix(c(fn(thetal) - fp(thetal), 0, 
                  0, fn(thetar) - fp(thetar)), 
                ncol = 2, nrow = 2) # nooit 0?
  
  cov_ca <- matrix(c(Vaa[1,1] * fp(thetal), 
                     Vaa[2,2] * fp(thetal) * (thetal - mup) / (2 * omegap),
                     Vaa[1,1] * fp(thetar),
                     Vaa[2,2] * fp(thetar) * (thetar - mup) / (2 * omegap)),
                   ncol = 2, nrow = 2)
  Vca <- -cov_ca %*% solve(t(Eas))
  Vac <- t(Vca)

  cov_cb <- matrix(c(Vbb[1,1] * fn(thetal), 
                     Vbb[2,2] * fn(thetal) * (thetal - mun) / (2 * omegan),
                     Vbb[1,1] * fn(thetar),
                     Vbb[2,2] * fn(thetar) * (thetar - mun) / (2 * omegan)),
                   ncol = 2, nrow = 2)
  Vcb <- cov_cb %*% solve(t(Eas))
  Vbc <- t(Vcb)

  Vcc <-
    matrix(
      c(
        (Vfp(thetal) + Vfn(thetal)) / (fp(thetal) - fn(thetal))^2, 0,
        0, (Vfp(thetar) + Vfn(thetar)) / (fp(thetar) - fn(thetar))^2
      ),
      ncol = 2, nrow = 2
    )

  Va <- rbind(Vaa, Vab, Vac)
  Vb <- rbind(Vba, Vbb, Vbc)
  Vc <- rbind(Vca, Vcb, Vcc)
  V <- cbind(Va, Vb, Vc)

  return(V)
}
