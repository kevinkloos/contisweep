% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize_variance.R
\name{optimize_variance}
\alias{optimize_variance}
\title{Optimize the variance of Continuous Sweep}
\usage{
optimize_variance(tr, te, var_type = "advanced")
}
\arguments{
\item{tr}{List containing training parameters:
\itemize{
  \item mup - Mean of positive class distribution
  \item omegap - Variance of positive class distribution
  \item mun - Mean of negative class distribution  
  \item omegan - Variance of negative class distribution
  \item ntrain - Number of training instances
  \item prevtrain - Prevalence of positive class in training data
}}

\item{te}{Numeric vector of test scores (target population)}

\item{var_type}{The type of variance estimate:
\itemize{
  \item basic - assumes training parameters Fp and Fn as known.
  \item advanced - assumes training parameters Fp and Fn as estimated.
}}
}
\value{
A list containing:
\itemize{
  \item variance - The estimated variance given optimal pdelta.
  \item pdelta - Optimal pdelta estimate.
  \item barriers - Two decision boundaries.
  \item ac_max - An initial adjusted count estimate using AC-MAX.
  \item var_type - The type of variance estimate (basic or advanced).
}
}
\description{
Optimize the variance of Continuous Sweep
}
\examples{
# Example training parameters
tr <- list(mup = 1, omegap = 1, mun = 0, omegan = 0.5, ntrain = 1000, prevtrain = 0.7)
# Generate test scores
set.seed(123)
te <- c(stats::rnorm(500, mean = tr$mup, sd = sqrt(tr$omegap)), 
        stats::rnorm(500, mean = tr$mun, sd = sqrt(tr$omegan)))
# Apply the function
optimize_variance(tr, te, var_type = "advanced")

}
