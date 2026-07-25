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

test_that("EM with alpha_fixed holds alpha and still recovers Pi", {
  d = rl_dgp(J = 4, N = 1e5, alpha1 = 0.2, alpha2 = 0.3, seed = 25)
  out = misclassifyr_rl_em(d$tab, J = 4, K = 4, alpha_fixed = c(0.2, 0.3))
  expect_true(out$converged)
  expect_equal(unname(out$alpha), c(0.2, 0.3))
  Pi_hat = matrix(0, 4, 4)
  Pi_hat[cbind(out$Pi$j, out$Pi$i)] = out$Pi$p
  expect_lt(max(abs(Pi_hat - d$Pi_true)), 0.02)
  expect_error(misclassifyr_rl_em(d$tab, J = 4, K = 4, alpha_fixed = c(0.2)),
               "length-2")
})

test_that("EM with T2 separates true transitions from linkage error", {
  set.seed(31)
  J = 50; N = 2e5; a1 = 0.10; a2 = 0.15
  xi = sample.int(J, N, replace = TRUE)
  jj = ifelse(runif(N) < 0.8, xi, ((xi + sample.int(4, N, replace = TRUE) - 1) %% J) + 1)
  rho = tabulate(jj, J) / N
  mv2 = sample.int(3, N, replace = TRUE)
  l_true = ifelse(runif(N) < 0.85, jj, ((jj - 1 + mv2) %% J) + 1)
  y1 = ifelse(runif(N) < a1, sample.int(J, N, replace = TRUE, prob = rho), jj)
  y2 = ifelse(runif(N) < a2, sample.int(J, N, replace = TRUE, prob = rho), l_true)
  tab = dplyr::count(data.frame(X = xi, Y1 = y1, Y2 = y2), X, Y1, Y2, name = "n")
  Tdf = do.call(rbind, lapply(1:J, function(j){
    data.frame(j = j, l = c(j, ((j - 1 + 1:3) %% J) + 1),
               t = c(0.85, rep(0.05, 3)))
  }))
  out_T = misclassifyr_rl_em(tab, J, J, T2 = Tdf)
  expect_true(out_T$converged)
  expect_lt(abs(out_T$alpha["alpha1"] - a1), 0.02)
  expect_lt(abs(out_T$alpha["alpha2"] - a2), 0.02)
  # Ignoring T (identity model) misreads true moves as link failures,
  # inflating alpha2 by roughly the true move rate - the design point.
  out_I = misclassifyr_rl_em(tab, J, J)
  expect_gt(unname(out_I$alpha["alpha2"]), a2 + 0.08)
  # T2 validation
  bad = Tdf; bad$t[1] = 0.5
  expect_error(misclassifyr_rl_em(tab, J, J, T2 = bad), "sum to one")
})
