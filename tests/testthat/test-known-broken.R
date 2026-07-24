# Formerly known-broken areas, fixed in the A2 bug-fix pass (2026-07-23);
# these tests now pin the corrected behavior.

test_that("Pi_to_beta bayesian = TRUE summarizes the posterior of beta", {
  # A tiny fake posterior with two draws of a
  # diagonal Pi should give posterior_beta == c(1, 1), median 1, sd 0.
  fake_posterior = do.call(rbind, lapply(1:2, function(d) {
    data.frame(
      Pi_hat = c(diag(3) / 3),
      X_name = as.character(rep(1:3, each = 3)),
      Y_name = as.character(rep(1:3, times = 3)),
      X_val = rep(1:3, each = 3),
      Y_val = rep(1:3, times = 3),
      draw = d
    )
  }))
  out = Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = FALSE, bayesian = TRUE,
                   posterior_Pi = fake_posterior)
  expect_equal(unname(out$posterior_beta), c(1, 1))
  expect_equal(unname(out$posterior_beta_med), 1)
  expect_equal(unname(out$posterior_beta_sd), 0)
})

test_that("make_empirical_Delta_RL_common_alpha single-tab branch matches list branch", {
  set.seed(3)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1e4)
  tab = syn$tab[[1]]

  single = make_empirical_Delta_RL_common_alpha(tab)
  aslist = make_empirical_Delta_RL_common_alpha(list(tab))
  psi = c(-1, -0.5)

  expect_equal(single$model_to_Delta(psi), aslist$model_to_Delta[[1]](psi))
  expect_equal(single$log_prior_Delta(psi), sum(dlogis(psi, log = TRUE)))
  # Rows of Delta (given Y*) must be distributions over (Y1, Y2)
  D = matrix(single$model_to_Delta(psi), nrow = 3)
  expect_equal(unname(rowSums(D)), rep(1, 3))
})

test_that("log_prior_Delta_NP_ind responds to the Delta2 block of psi", {
  J = 3
  psi_a = rep(0, 2 * J * (J - 1))
  psi_b = c(rep(0, J * (J - 1)), seq(-1, 1, length.out = J * (J - 1)))
  expect_false(
    isTRUE(all.equal(log_prior_Delta_NP_ind(psi_a), log_prior_Delta_NP_ind(psi_b)))
  )
})

test_that("misclassifyr MH acceptance uses the symmetric-proposal rule", {
  # Short-chain check that the (corrected, symmetric-rule) sampler targets
  # the posterior: posterior mean of Pi near the MLE on a well-identified DGP.
  set.seed(5)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1e5,
                       dgp_delta = "Record Linkage, independent, 10 - 30%")
  tab = syn$tab[[1]]
  out = suppressWarnings(suppressMessages(misclassifyr(
    tab = tab, J = 3, K = 3,
    X_names = as.character(1:3),
    Y1_names = as.character(1:3),
    Y2_names = as.character(1:3),
    model_to_Delta = model_to_Delta_RL_ind,
    log_prior_Delta = log_prior_Delta_RL_ind,
    X_vals = 1:3, Y_vals = 1:3,
    mle = TRUE, bayesian = TRUE, makeplots = FALSE,
    n_mcmc_draws = 12000, n_burnin = 6000
  )))
  # Group with Y_name varying fastest to match Pi_hat_mle's element order
  # (eta_hat_to_Pi_hat emits Y* fastest within X blocks).
  posterior_Pi_mean = aggregate(Pi_hat ~ Y_name + X_name,
                                data = out$posterior_Pi, FUN = mean)
  expect_lt(max(abs(posterior_Pi_mean$Pi_hat - out$Pi_hat_mle)), 0.05)
})

test_that("Pi_to_beta bayesian list branch pools cells and computes HPD outputs", {
  # Two identical control cells with a diagonal Pi over three draws: pooled
  # and per-cell posterior betas are exactly 1, and with an ll_history the
  # CCT outputs are populated (3 draws < 20, so the clamped cutoff keeps all).
  make_cell = function() do.call(rbind, lapply(1:3, function(d)
    data.frame(
      Pi_hat = c(diag(3) / 3),
      X_name = as.character(rep(1:3, each = 3)),
      Y_name = as.character(rep(1:3, times = 3)),
      X_val = rep(1:3, each = 3),
      Y_val = rep(1:3, times = 3),
      draw = d
    )))
  posterior = list(make_cell(), make_cell())
  llh = list(data.frame(ll = c(-10, -1, -2), draw = 1:3),
             data.frame(ll = c(-10, -1, -2), draw = 1:3))

  out = Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = FALSE, bayesian = TRUE,
                   W_weights = c(0.5, 0.5),
                   posterior_Pi = posterior, ll_history = llh)

  expect_equal(unname(out$posterior_beta), c(1, 1, 1))
  expect_length(out$posterior_betas, 2)
  expect_true(all(abs(unlist(out$posterior_betas) - 1) < 1e-8))
  expect_equal(out$HPD_draws, 1:3)
  expect_equal(unname(out$HPDCI), c(1, 1))
})

test_that("Pi_to_beta validates ll_history structure against posterior_Pi", {
  fake_posterior = data.frame(
    Pi_hat = c(diag(3) / 3),
    X_name = as.character(rep(1:3, each = 3)),
    Y_name = as.character(rep(1:3, times = 3)),
    X_val = rep(1:3, each = 3),
    Y_val = rep(1:3, times = 3),
    draw = 1
  )
  expect_error(
    Pi_to_beta(X_vals = 1:3, Y_vals = 1:3, mle = FALSE, bayesian = TRUE,
               posterior_Pi = fake_posterior,
               ll_history = list(data.frame(ll = -1, draw = 1))),
    "same structure"
  )
})

test_that("misclassifyr rejects a burn-in that thinning does not divide", {
  set.seed(7)
  syn = synthetic_data(J = 3, K = 3, I = 1, sample_size = 1e3)
  expect_error(
    suppressWarnings(misclassifyr(
      tab = syn$tab[[1]], J = 3, K = 3,
      X_names = as.character(1:3),
      Y1_names = as.character(1:3),
      Y2_names = as.character(1:3),
      X_vals = 1:3, Y_vals = 1:3,
      mle = FALSE, bayesian = TRUE, makeplots = FALSE,
      n_mcmc_draws = 12000, n_burnin = 6001, thinning_rate = 2
    )),
    "multiple of `thinning_rate`"
  )
})
