# A3.2: EM estimator for the record-linkage model at scale.

# Simulate from the exact model: latent (X, Y*) with a diagonal-heavy sparse
# Pi; each measure equals Y* w.p. 1-alpha_m, else an independent draw from
# the latent marginal rho.
rl_dgp = function(J, N, alpha1, alpha2, stay = 0.8, seed = 1){
  set.seed(seed)
  xi = sample.int(J, N, replace = TRUE)
  moves = sample.int(3, N, replace = TRUE)
  j = ifelse(runif(N) < stay, xi, ((xi - 1 + moves) %% J) + 1)
  rho = tabulate(j, J) / N
  y1 = ifelse(runif(N) < alpha1, sample.int(J, N, replace = TRUE, prob = rho), j)
  y2 = ifelse(runif(N) < alpha2, sample.int(J, N, replace = TRUE, prob = rho), j)
  tab = dplyr::count(data.frame(X = xi, Y1 = y1, Y2 = y2), X, Y1, Y2, name = "n")
  Pi_true = unclass(table(factor(j, 1:J), factor(xi, 1:J))) / N
  list(tab = tab, Pi_true = Pi_true, rho = rho)
}

test_that("EM recovers alpha, rho, and Pi on a small-J DGP", {
  d = rl_dgp(J = 4, N = 2e5, alpha1 = 0.15, alpha2 = 0.30, seed = 21)
  out = misclassifyr_rl_em(d$tab, J = 4, K = 4)
  expect_true(out$converged)
  # loglik is monotone (EM property)
  expect_true(all(diff(out$loglik_trace) > -1e-6))
  expect_lt(abs(out$alpha["alpha1"] - 0.15), 0.02)
  expect_lt(abs(out$alpha["alpha2"] - 0.30), 0.02)
  expect_lt(max(abs(out$rho1 - d$rho)), 0.03)
  # Pi on its support vs truth
  Pi_hat = matrix(0, 4, 4)
  Pi_hat[cbind(out$Pi$j, out$Pi$i)] = out$Pi$p
  expect_lt(max(abs(Pi_hat - d$Pi_true)), 0.02)
})

test_that("EM with fixed rho recovers alpha", {
  d = rl_dgp(J = 4, N = 1e5, alpha1 = 0.25, alpha2 = 0.25, seed = 22)
  out = misclassifyr_rl_em(d$tab, J = 4, K = 4, rho1 = d$rho, rho2 = d$rho)
  expect_true(out$converged)
  expect_lt(abs(out$alpha["alpha1"] - 0.25), 0.02)
  expect_lt(abs(out$alpha["alpha2"] - 0.25), 0.02)
})

test_that("EM handles county-scale J without forming dense objects", {
  d = rl_dgp(J = 300, N = 4e5, alpha1 = 0.20, alpha2 = 0.35, seed = 23)
  t0 = Sys.time()
  out = misclassifyr_rl_em(d$tab, J = 300, K = 300, maxit = 1000)
  runtime = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_true(out$converged)
  expect_lt(abs(out$alpha["alpha1"] - 0.20), 0.03)
  expect_lt(abs(out$alpha["alpha2"] - 0.35), 0.03)
  # Pi accurate on its support
  p_true = d$Pi_true[cbind(out$Pi$j, out$Pi$i)]
  expect_gt(cor(out$Pi$p, p_true), 0.98)
  # and the whole thing is fast enough for MC use
  expect_lt(runtime, 120)
})

test_that("EM validates malformed inputs", {
  d = rl_dgp(J = 4, N = 1e3, alpha1 = 0.2, alpha2 = 0.2, seed = 24)
  expect_error(misclassifyr_rl_em(d$tab[, c("X", "Y1", "n")], J = 4, K = 4),
               "columns")
  bad = d$tab; bad$Y1[1] = 9
  expect_error(misclassifyr_rl_em(bad, J = 4, K = 4), "integer codes")
  expect_error(misclassifyr_rl_em(d$tab, J = 4, K = 4, rho1 = rep(0.5, 4)),
               "probability vectors")
})
