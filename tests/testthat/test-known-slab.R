# Two-component mixture with a KNOWN rival distribution: the 1935-40
# migration-question design, where the population joint is observed in the
# full count rather than inferred from a repeated measure.

sim_slab = function(J, N, alpha, stay = 0.7, seed = 1){
  set.seed(seed)
  x <- sample.int(J, N, replace = TRUE)
  rho <- rep(1/J, J)
  y_correct <- ifelse(runif(N) < stay, x, (x %% J) + 1)
  y <- ifelse(runif(N) < alpha, sample.int(J, N, replace = TRUE), y_correct)
  tab <- aggregate(list(n = rep(1, N)), by = list(x = x, y = y), FUN = sum)
  T_cond <- do.call(rbind, lapply(1:J, function(i)
    data.frame(x = i, y = c(i, (i %% J) + 1), t = c(stay, 1 - stay))))
  list(tab = tab, T_cond = T_cond, rho = rho)
}

test_that("known-slab estimator recovers alpha and its standard error", {
  d = sim_slab(J = 12, N = 5e4, alpha = 0.20, seed = 11)
  fit = misclassifyr_known_slab(d$tab, d$T_cond, d$rho)
  expect_lt(abs(fit$alpha - 0.20), 0.02)
  # Wald interval covers the truth, and the SE is of a sensible magnitude
  expect_lt(fit$ci[1], 0.20); expect_gt(fit$ci[2], 0.20)
  expect_lt(fit$se, 0.01)
})

test_that("recovery holds across the plausible range of alpha", {
  for(a in c(0.05, 0.35)){
    d = sim_slab(J = 12, N = 5e4, alpha = a, seed = round(100 * a))
    fit = misclassifyr_known_slab(d$tab, d$T_cond, d$rho)
    expect_lt(abs(fit$alpha - a), 0.025)
  }
})

test_that("impossible cells give a model-free lower bound on alpha", {
  # With a tight T (no moves at all), every off-diagonal pair is impossible
  # under a correct link, so the impossible share estimates alpha directly.
  d = sim_slab(J = 12, N = 3e4, alpha = 0.25, stay = 1, seed = 13)
  fit = misclassifyr_known_slab(d$tab, d$T_cond, d$rho)
  expect_lt(abs(fit$share_impossible - 0.25 * (11/12)), 0.02)
  expect_lt(fit$share_impossible, fit$alpha + 1e-6)
})

test_that("known-slab validates its inputs", {
  d = sim_slab(J = 5, N = 1e3, alpha = 0.2, seed = 14)
  expect_error(misclassifyr_known_slab(d$tab[, c("x", "n")], d$T_cond, d$rho),
               "columns")
  expect_error(misclassifyr_known_slab(d$tab, d$T_cond, rep(0.5, 5)),
               "sum to one")
  bad = d$T_cond; bad$t[1] = 0.9
  expect_error(misclassifyr_known_slab(d$tab, bad, d$rho), "sum to one")
})
