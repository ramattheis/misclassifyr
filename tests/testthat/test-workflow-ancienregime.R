# End-to-end smoke test adapted from misc/workflow.R (continuous-outcome
# branch): prep_misclassification_data -> misclassifyr (MLE only, coarse
# tolerance) -> Pi_to_beta on the bundled ancienregime data.
#
# Note: misc/workflow.R passes `estimate_beta = TRUE` to misclassifyr(), which
# is not an argument of the function and errors; the working call omits it.

test_that("ancienregime workflow produces finite estimates and SEs", {
  skip_on_cran()

  env = new.env()
  data("ancienregime", package = "misclassifyr", envir = env)
  ancienregime = env$ancienregime

  set.seed(1)
  inputs = prep_misclassification_data(
    data = ancienregime,
    outcome_1 = "son_income_1780",
    outcome_2 = "son_income_1770",
    regressor = "father_income_1750",
    outcome_1_bin = "son_occupation_1780",
    outcome_2_bin = "son_occupation_1770",
    regressor_bin = "father_occupation_1750",
    weights = "linked_weight",
    record_vals = TRUE,
    round_vals = 0
  )

  # Six occupation bins for both generations, balanced J^2 * K table
  expect_equal(inputs$J, 6)
  expect_equal(inputs$K, 6)
  expect_equal(nrow(inputs$tab), 6^2 * 6)
  expect_equal(sum(inputs$tab$n), sum(ancienregime$linked_weight))
  expect_length(inputs$X_vals, 6)
  expect_length(inputs$Y_vals, 6)
  expect_false(is.unsorted(inputs$X_vals))

  # The J = 6 nonparametric model has 95 free parameters; the stability check
  # can warn about inconsistent optima, which we tolerate in a smoke test.
  out = suppressWarnings(suppressMessages(misclassifyr(
    tab = inputs$tab,
    J = inputs$J,
    K = inputs$K,
    model_to_Delta = model_to_Delta_NP_ind,
    X_names = inputs$X_names,
    Y1_names = inputs$Y1_names,
    Y2_names = inputs$Y2_names,
    X_vals = inputs$X_vals,
    Y_vals = inputs$Y_vals,
    W_names = inputs$W_names,
    mle = TRUE, bayesian = FALSE, makeplots = FALSE,
    optim_tol = 1e-5
  )))

  expect_true(is.finite(out$log_likelihood_mle))
  Pi_hat = matrix(out$Pi_hat_mle, nrow = 6)
  expect_equal(sum(Pi_hat), 1, tolerance = 1e-6)
  expect_true(all(Pi_hat >= 0))

  beta_out = Pi_to_beta(
    X_vals = inputs$X_vals,
    Y_vals = inputs$Y_vals,
    mle = TRUE, bayesian = FALSE,
    Pi_mle = out$Pi_hat_mle,
    cov_Pi = out$cov_Pi_mle
  )

  expect_true(is.finite(unname(beta_out$beta_hat_mle)))
  expect_true(is.finite(unname(beta_out$se_beta_mle)))
  expect_gt(beta_out$se_beta_mle, 0)
  # Loose sanity band around the pinned estimate (~0.72 with this seed);
  # the naive OLS on the same data is ~0.60.
  expect_gt(unname(beta_out$beta_hat_mle), 0.3)
  expect_lt(unname(beta_out$beta_hat_mle), 1.2)
})
