# Exact tests for the Pi -> beta mapping (MLE path; the bayesian = TRUE
# branch is covered in test-known-broken.R, fixed in the A2 pass).

test_that("Pi_to_beta_inner recovers beta exactly for a diagonal Pi", {
  Pi = c(diag(3) / 3)
  expect_equal(Pi_to_beta_inner(Pi, 1:3, 1:3, 1), 1)
  expect_equal(Pi_to_beta_inner(Pi, 1:3, 2 * (1:3), 1), 2)
  expect_equal(Pi_to_beta_inner(Pi, 1:3, 5 - (1:3), 1), -1)
})

test_that("Pi_to_beta_inner clamps negative entries and renormalizes", {
  raw = c(0.6, -0.2, 0.1, 0.3)
  clamped = c(0.6, 0, 0.1, 0.3)
  expect_equal(
    Pi_to_beta_inner(raw, 1:2, 1:2, 1),
    Pi_to_beta_inner(clamped / sum(clamped), 1:2, 1:2, 1)
  )
})

test_that("Pi_to_beta_inner aggregates identical control cells to the pooled beta", {
  Pi = c(diag(3) / 3)
  b = Pi_to_beta_inner(list(Pi, Pi), list(1:3, 1:3), list(1:3, 1:3), c(0.5, 0.5))
  expect_equal(b, 1)
})

test_that("se_beta_deltamethod is zero under zero covariance and positive otherwise", {
  Pi = c(diag(3) / 3)
  expect_equal(se_beta_deltamethod(Pi, matrix(0, 9, 9), 1:3, 1:3, 1), 0)
  se = se_beta_deltamethod(Pi, diag(1e-4, 9), 1:3, 1:3, 1)
  expect_true(is.finite(se))
  expect_gt(se, 0)
})

test_that("Pi_to_beta mle path returns beta, SE, and per-cell estimates", {
  Pi = c(diag(3) / 3)
  zero_cov = matrix(0, 9, 9)

  # Single population
  out = Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = TRUE, bayesian = FALSE,
                   Pi_mle = Pi, cov_Pi = zero_cov)
  expect_named(out, c("beta_hat_mle", "se_beta_mle", "betas_hat_mle",
                      "se_betas_mle", "posterior_beta", "posterior_beta_med",
                      "posterior_beta_sd", "posterior_betas",
                      "posterior_betas_med", "posterior_betas_sd",
                      "HPD_draws", "HPDCI"))
  expect_equal(unname(out$beta_hat_mle), 1)
  expect_equal(unname(out$se_beta_mle), 0)
  expect_true(is.na(out$betas_hat_mle))
  expect_true(is.na(out$posterior_beta))

  # Two control cells
  out2 = Pi_to_beta(X_vals = list(1:3, 1:3), Y_vals = list(1:3, 1:3),
                    W_weights = c(600, 400), mle = TRUE, bayesian = FALSE,
                    Pi_mle = list(Pi, Pi), cov_Pi = list(zero_cov, zero_cov))
  expect_equal(unname(out2$beta_hat_mle), 1)
  expect_equal(unname(out2$se_beta_mle), 0)
  expect_equal(unname(out2$betas_hat_mle), c(1, 1))
  expect_equal(unname(out2$se_betas_mle), c(0, 0))
})

test_that("Pi_to_beta input validation catches malformed calls", {
  Pi = c(diag(3) / 3)
  zero_cov = matrix(0, 9, 9)
  # MLE requested without estimates
  expect_error(Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = TRUE),
               "should be provided")
  # List of Pis without cell weights
  expect_error(
    Pi_to_beta(X_vals = list(1:3), Y_vals = list(1:3), mle = TRUE,
               Pi_mle = list(Pi), cov_Pi = list(zero_cov)),
    "W_weights"
  )
})
