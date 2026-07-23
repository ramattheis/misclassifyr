# Regression tests pinning current MLE behavior of misclassifyr() on small
# synthetic DGPs. Tolerances are deliberately loose: these tests are meant to
# catch structural breakage during refactoring, not to certify efficiency.

# The full set of elements returned by misclassifyr() (v0.2.3): 20 estimation
# outputs + 4 plot slots.
misclassifyr_output_names = c(
  "Pi_hat_mle", "Delta_hat_mle", "cov_Pi_mle", "eta_hat_mle",
  "log_likelihood_mle", "W_weights", "optim_counts", "model_to_Pi_jacobian",
  "eta_hessian_mle", "fisher_info_err", "inconsistency_mle",
  "posterior_Pi", "posterior_Delta", "posterior_eta", "ll_history",
  "accepted_proposals", "trace_plots_eta", "trace_plots_Pi",
  "trace_plots_Delta", "misclassification_inputs",
  "Pi_hat_mle_plot", "Delta_hat_mle_plot", "Pi_hat_posterior_plot",
  "Delta_hat_posterior_plot"
)

test_that("MLE with model_to_Delta_NP_ind recovers Pi and beta on synthetic data", {
  set.seed(42)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1e5,
                       dgp_delta = "Nonparametric, independent, strong diagonal",
                       dgp_pi = "Exponential")
  tab = syn$tab[[1]]
  Pi_true = syn$Pi[[1]]

  out = suppressWarnings(suppressMessages(misclassifyr(
    tab = tab, J = 3, K = 3,
    X_names = as.character(1:3),
    Y1_names = as.character(1:3),
    Y2_names = as.character(1:3),
    model_to_Delta = model_to_Delta_NP_ind,
    X_vals = 1:3, Y_vals = 1:3,
    mle = TRUE, bayesian = FALSE, makeplots = FALSE,
    optim_tol = 1e-6
  )))

  # Output object structure (snapshot of names)
  expect_named(out, misclassifyr_output_names)

  # Convergence / bookkeeping
  expect_true(is.finite(out$log_likelihood_mle))
  expect_equal(out$fisher_info_err, "Fisher information matrix is invertible.")
  expect_equal(out$W_weights, sum(tab$n))
  expect_true(is.numeric(out$optim_counts))
  expect_true(is.numeric(out$inconsistency_mle))
  expect_length(out$eta_hat_mle, 3 * 3 - 1 + 2 * 3 * (3 - 1))
  expect_equal(dim(out$cov_Pi_mle), c(9, 9))
  expect_true(all(diag(out$cov_Pi_mle) >= 0))
  # Bayesian and plot slots are NA when not requested
  expect_true(is.na(out$posterior_Pi))
  expect_true(is.na(out$accepted_proposals))
  expect_true(is.na(out$Pi_hat_mle_plot[1]))

  # Pi_hat: on the simplex and close to truth (loose tolerance)
  Pi_hat = matrix(out$Pi_hat_mle, nrow = 3)
  expect_equal(sum(Pi_hat), 1, tolerance = 1e-8)
  expect_true(all(Pi_hat >= 0))
  expect_lt(max(abs(Pi_hat - Pi_true)), 0.05)

  # Delta_hat: each row (given Y*) is a distribution over (Y1, Y2)
  Delta_hat = matrix(out$Delta_hat_mle, nrow = 3)
  expect_equal(unname(rowSums(Delta_hat)), rep(1, 3), tolerance = 1e-8)
  expect_true(all(Delta_hat >= 0))

  # Recovered beta close to the true Pi-implied beta
  beta_true = Pi_to_beta_inner(c(Pi_true), 1:3, 1:3, 1)
  bo = Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = TRUE, bayesian = FALSE,
                  Pi_mle = out$Pi_hat_mle, cov_Pi = out$cov_Pi_mle)
  expect_lt(abs(unname(bo$beta_hat_mle) - beta_true), 0.1)
  expect_true(is.finite(bo$se_beta_mle))
  expect_gt(bo$se_beta_mle, 0)
})

test_that("MLE with model_to_Delta_RL_ind recovers Pi, Delta, and beta on RL data", {
  set.seed(7)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1e5,
                       dgp_delta = "Record Linkage, independent, 10 - 30%",
                       dgp_pi = "Exponential")
  tab = syn$tab[[1]]
  Pi_true = syn$Pi[[1]]
  Delta_true = syn$Delta[[1]]

  out = suppressWarnings(suppressMessages(misclassifyr(
    tab = tab, J = 3, K = 3,
    X_names = as.character(1:3),
    Y1_names = as.character(1:3),
    Y2_names = as.character(1:3),
    model_to_Delta = model_to_Delta_RL_ind,
    X_vals = 1:3, Y_vals = 1:3,
    mle = TRUE, bayesian = FALSE, makeplots = FALSE,
    optim_tol = 1e-6
  )))

  expect_named(out, misclassifyr_output_names)
  expect_true(is.finite(out$log_likelihood_mle))
  expect_length(out$eta_hat_mle, 3 * 3 - 1 + 4 * 3 - 2)

  Pi_hat = matrix(out$Pi_hat_mle, nrow = 3)
  expect_equal(sum(Pi_hat), 1, tolerance = 1e-8)
  expect_lt(max(abs(Pi_hat - Pi_true)), 0.05)

  Delta_hat = matrix(out$Delta_hat_mle, nrow = 3)
  expect_lt(max(abs(Delta_hat - Delta_true)), 0.05)

  beta_true = Pi_to_beta_inner(c(Pi_true), 1:3, 1:3, 1)
  bo = Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = TRUE, bayesian = FALSE,
                  Pi_mle = out$Pi_hat_mle, cov_Pi = out$cov_Pi_mle)
  expect_lt(abs(unname(bo$beta_hat_mle) - beta_true), 0.05)
  expect_true(is.finite(bo$se_beta_mle))
  expect_gt(bo$se_beta_mle, 0)
})

test_that("misclassifyr input validation catches malformed calls", {
  set.seed(11)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1e4)
  tab = syn$tab[[1]]
  nm = as.character(1:3)

  # Neither estimator requested
  expect_error(
    misclassifyr(tab, 3, 3, nm, nm, nm, mle = FALSE, bayesian = FALSE),
    "Either `bayesian` or `mle`"
  )

  # Missing count column
  bad = tab
  bad$n = NULL
  expect_error(
    suppressMessages(misclassifyr(bad, 3, 3, nm, nm, nm, makeplots = FALSE)),
    "four columns"
  )

  # Negative counts
  bad = tab
  bad$n[1] = -5
  expect_error(
    suppressMessages(misclassifyr(bad, 3, 3, nm, nm, nm, makeplots = FALSE)),
    "non-negative"
  )

  # NAs in the table
  bad = tab
  bad$Y1[1] = NA
  expect_error(
    suppressMessages(misclassifyr(bad, 3, 3, nm, nm, nm, makeplots = FALSE)),
    "NA values"
  )

  # Unbalanced table
  expect_error(
    suppressMessages(misclassifyr(tab[-1, ], 3, 3, nm, nm, nm, makeplots = FALSE)),
    "balanced"
  )
})
