context("CS Prevalence Estimation")

test_that("basic functionality works", {
  tr <- list(mup = 2, omegap = 1, mun = -2, omegan = 1.5)
  te <- c(rnorm(500, mean = tr$mup, sd = sqrt(tr$omegap)), 
          rnorm(500, mean = tr$mun, sd = sqrt(tr$omegan)))
  result <- CS(tr, te, barriers = c(-3, 3))
  expect_true(is.numeric(result))
  expect_true(result >= 0 & result <= 1)
})

test_that("handles invalid barrier inputs", {
  tr <- list(mup = 0, omegap = 1, mun = 0, omegan = 1)
  te <- rnorm(100)
  # Test with non-numeric barriers
  expect_error(CS(tr, te, barriers = c("a", "b")))
  # Test with single barrier value
  expect_error(CS(tr, te, barriers = 1))
})

test_that("input validation works", {
  tr <- list(mup = "a", omegap = 1, mun = -2, omegan = 1.5)
  te <- rnorm(100)
  expect_error(CS(tr, te, barriers = c(-3, 3)))
})

test_that("warning triggered for out-of-bounds prevalence", {
  tr <- list(mup = 1, omegap = 0.1, mun = -1, omegan = 0.1)
  te <- rnorm(1000, 1, 1) # Strongly positive skewed sample
  expect_warning(CS(tr, te, barriers = c(-5, 5)), "Estimated prevalence outside 0-1")
})

