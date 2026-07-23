# Known-broken areas, documented as skipped tests. Each test states the
# intended post-fix behavior; un-skip (and adjust if the fixed signature
# differs) once the corresponding bug is fixed. Do NOT write green tests
# around current behavior here: these code paths error or are wrong by
# construction, and are scheduled for surgery.

test_that("Pi_to_beta bayesian = TRUE summarizes the posterior of beta", {
  skip(paste(
    "Known broken: Pi_to_beta's bayesian = TRUE branch references undefined",
    "objects (`tab`, `misclassification_output`) instead of its own arguments,",
    "so it errors on any call. Post-fix it should map posterior draws of Pi",
    "(and an ll_history for the CCT HPD interval) to a posterior for beta."
  ))

  # Intended post-fix behavior: a tiny fake posterior with two draws of a
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
  skip(paste(
    "Known broken: the non-list branch of make_empirical_Delta_RL_common_alpha",
    "builds Delta with outer(FY, rep(1 - alpha)) where the list branch",
    "(correctly) uses rep(alpha), and its log_prior_Delta contains a 'retun('",
    "typo. Post-fix, a single tab and a length-one list must agree."
  ))

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
  skip(paste(
    "Known broken: log_prior_Delta_NP_ind computes log_prior_Delta2 from",
    "Delta1 (copy-paste bug), so the second half of psi never affects the",
    "prior. Post-fix, changing only the Delta2 block must change the value."
  ))

  J = 3
  psi_a = rep(0, 2 * J * (J - 1))
  psi_b = c(rep(0, J * (J - 1)), seq(-1, 1, length.out = J * (J - 1)))
  expect_false(
    isTRUE(all.equal(log_prior_Delta_NP_ind(psi_a), log_prior_Delta_NP_ind(psi_b)))
  )
})

test_that("misclassifyr MH acceptance uses the symmetric-proposal rule", {
  skip(paste(
    "Known broken: the Metropolis-Hastings acceptance in misclassifyr()'s",
    "Gibbs sampler subtracts `gibbs_jump` from the log acceptance ratio.",
    "The proposal is a symmetric normal increment, so no correction term",
    "belongs there; the current chain is biased toward negative jumps.",
    "Post-fix, add a short-chain test that the sampler targets the posterior",
    "(e.g. posterior mean of Pi near the MLE on a well-identified DGP)."
  ))

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
  posterior_Pi_mean = aggregate(Pi_hat ~ X_name + Y_name,
                                data = out$posterior_Pi, FUN = mean)
  expect_lt(max(abs(posterior_Pi_mean$Pi_hat - out$Pi_hat_mle)), 0.05)
})
