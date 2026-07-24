# A3.1: sparse tabulations must be exactly equivalent to dense ones, and
# lambda_dd is a real argument.

test_that("sparse prep drops only zero cells and indexes them correctly", {
  set.seed(11)
  n = 2000
  md = data.frame(
    f_occ  = sample(1:4, n, replace = TRUE),
    s_occ1 = sample(1:4, n, replace = TRUE),
    s_occ2 = sample(1:4, n, replace = TRUE)
  )
  nm = as.character(1:4)
  dense = prep_misclassification_data(md, "s_occ1", "s_occ2", "f_occ",
                                      X_names = nm, Y1_names = nm, Y2_names = nm)
  sp    = prep_misclassification_data(md, "s_occ1", "s_occ2", "f_occ",
                                      X_names = nm, Y1_names = nm, Y2_names = nm,
                                      sparse = TRUE)
  expect_true(all(sp$tab$n > 0))
  expect_equal(sum(sp$tab$n), sum(dense$tab$n))
  # cell_idx points at the matching row of the dense balanced layout
  expect_equal(dense$tab$n[sp$tab$cell_idx], sp$tab$n)
  expect_equal(dense$tab$X[sp$tab$cell_idx], sp$tab$X)
  expect_equal(dense$tab$Y1[sp$tab$cell_idx], sp$tab$Y1)
  expect_equal(dense$tab$Y2[sp$tab$cell_idx], sp$tab$Y2)
})

test_that("loglikelihood is identical for dense and sparse tabs", {
  set.seed(12)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 300)
  tab = syn$tab[[1]]
  tab$n[c(2, 7, 19)] = 0  # force sparsity deterministically
  sp = tab; sp$cell_idx = seq_len(nrow(sp)); sp = sp[sp$n > 0, ]
  theta = c(model_to_Pi_NP(rnorm(8), 3), model_to_Delta_NP_ind(rnorm(12)))
  expect_equal(loglikelihood(theta, tab, 3, 3, lambda_dd = 100),
               loglikelihood(theta, sp,  3, 3, lambda_dd = 100))
})

test_that("misclassifyr returns identical MLEs on dense and sparse tabs", {
  set.seed(13)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 300)
  tab = syn$tab[[1]]
  sp = tab; sp$cell_idx = seq_len(nrow(sp)); sp = sp[sp$n > 0, ]
  args = list(J = 3, K = 3, X_names = as.character(1:3),
              Y1_names = as.character(1:3), Y2_names = as.character(1:3),
              X_vals = 1:3, Y_vals = 1:3, mle = TRUE, bayesian = FALSE,
              makeplots = FALSE)
  out_d = suppressWarnings(suppressMessages(
    do.call(misclassifyr, c(list(tab = tab), args))))
  out_s = suppressWarnings(suppressMessages(
    do.call(misclassifyr, c(list(tab = sp), args))))
  expect_equal(out_d$Pi_hat_mle, out_s$Pi_hat_mle, tolerance = 1e-8)
  expect_equal(out_d$Delta1_hat_mle, out_s$Delta1_hat_mle, tolerance = 1e-8)
})

test_that("lambda_dd argument is respected", {
  set.seed(14)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1000)
  tab = syn$tab[[1]]
  args = list(tab = tab, J = 3, K = 3, X_names = as.character(1:3),
              Y1_names = as.character(1:3), Y2_names = as.character(1:3),
              X_vals = 1:3, Y_vals = 1:3, mle = TRUE, bayesian = FALSE,
              makeplots = FALSE)
  out_default  = suppressWarnings(suppressMessages(do.call(misclassifyr, args)))
  out_explicit = suppressWarnings(suppressMessages(
    do.call(misclassifyr, c(args, list(lambda_dd = sum(tab$n)^2)))))
  out_zero     = suppressWarnings(suppressMessages(
    do.call(misclassifyr, c(args, list(lambda_dd = 0)))))
  # explicit default value reproduces the default exactly
  expect_equal(out_default$Pi_hat_mle, out_explicit$Pi_hat_mle)
  # lambda_dd = 0 still converges to finite estimates on a clean DGP
  expect_true(all(is.finite(out_zero$Pi_hat_mle)))
})
