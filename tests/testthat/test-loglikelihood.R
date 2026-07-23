# Exact, hand-computed unit tests for the cell-level likelihood machinery:
# softlog(), logit_link_volume(), loglikelihood(), and the diagonal-dominance
# penalty.

test_that("softlog matches log above the clamp and clamps below 1e-20", {
  expect_equal(softlog(1), 0)
  expect_equal(softlog(exp(2)), 2)
  expect_equal(softlog(0), log(1e-20))
  expect_equal(softlog(1e-30), log(1e-20))
  expect_equal(softlog(-5), log(1e-20))
  expect_equal(softlog(c(1, 0, 0.5)), c(0, log(1e-20), log(0.5)))
})

test_that("logit_link_volume matches hand-computed values", {
  # Scalar 0: p = 1/2, |Jacobian| = p(1-p) = 1/4, constant = log(2!)
  # => log(1/4) + log(2) = log(1/2)
  expect_equal(logit_link_volume(0), log(0.5))

  # c(0, 0): p = (1/3, 1/3), Jacobian = [[2/9, -1/9], [-1/9, 2/9]],
  # det = 1/27, constant = log(3!) => log(6/27) = log(2/9)
  expect_equal(logit_link_volume(c(0, 0)), log(2 / 9))
})

test_that("log_prior_Pi_NP and log_prior_Delta_RL_ind reduce to logit_link_volume", {
  expect_equal(log_prior_Pi_NP(c(0, 0)), logit_link_volume(c(0, 0)))
  # J = 2, psi = rep(0, 6): two scalar row-scale volumes at 0 (log(1/2) each),
  # column-scale terms are log(0.5) - log(0.5) = 0
  expect_equal(log_prior_Delta_RL_ind(rep(0, 6)), 2 * log(0.5))
})

test_that("loglikelihood reproduces a hand-computed small case", {
  J = 2
  K = 2
  Pi = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = J)      # rows Y*, cols X
  delta1 = matrix(c(0.8, 0.2, 0.3, 0.7), nrow = J)  # Pr(Y1 = r | Y* = c)
  delta2 = matrix(c(0.9, 0.1, 0.25, 0.75), nrow = J)
  Delta = do.call(cbind, lapply(1:J, function(j) diag(delta2[j, ]) %*% t(delta1)))

  tab = expand.grid(X = 1:K, Y1 = 1:J, Y2 = 1:J, KEEP.OUT.ATTRS = FALSE)
  tab$n = 1:8

  # Independent derivation of the cell probabilities:
  # Pr(X = x, Y1 = y1, Y2 = y2) = sum_ys Pi[ys, x] Pr(y1 | ys) Pr(y2 | ys)
  p = apply(tab[, c("X", "Y1", "Y2")], 1, function(row) {
    sum(sapply(1:J, function(ys)
      Pi[ys, row["X"]] * delta1[row["Y1"], ys] * delta2[row["Y2"], ys]))
  })

  expect_equal(
    loglikelihood(c(c(Pi), c(Delta)), tab, J, K, lambda_dd = 0),
    sum(tab$n * log(p))
  )
})

test_that("loglikelihood clamps zero-probability cells with positive counts", {
  J = 2
  K = 1
  Pi = matrix(c(1, 0), nrow = J)
  delta_id = diag(J)
  Delta = do.call(cbind, lapply(1:J, function(j) diag(delta_id[j, ]) %*% t(delta_id)))

  tab = expand.grid(X = 1, Y1 = 1:2, Y2 = 1:2, KEEP.OUT.ATTRS = FALSE)
  tab$n = c(5, 0, 0, 3)

  # Cell (Y1 = 2, Y2 = 2) has probability 0 but count 3: contributes
  # 3 * log(1e-20); cell (1, 1) has probability 1 and contributes 0.
  expect_equal(
    loglikelihood(c(c(Pi), c(Delta)), tab, J, K, lambda_dd = 0),
    3 * log(1e-20)
  )
})

test_that("diagonal-dominance penalty activates when off-diagonal dominates", {
  J = 2
  K = 1
  Pi = matrix(c(0.5, 0.5), nrow = J)
  tab = expand.grid(X = 1, Y1 = 1:2, Y2 = 1:2, KEEP.OUT.ATTRS = FALSE)
  tab$n = rep(1, 4)

  # In row Y* = 1, block Y2 = 1, the entry for Y1 = 2 (0.3) exceeds the
  # dominant entry for Y1 = 1 (0.2) by 0.1; all other comparisons are slack.
  Delta_bad = rbind(c(0.20, 0.30, 0.40, 0.10),
                    c(0.25, 0.25, 0.25, 0.25))
  ll0 = loglikelihood(c(c(Pi), c(Delta_bad)), tab, J, K, lambda_dd = 0)
  ll_pen = loglikelihood(c(c(Pi), c(Delta_bad)), tab, J, K, lambda_dd = 100)
  expect_equal(ll_pen, ll0 - 100 * 0.1^2)

  # A diagonally dominant Delta incurs no penalty even at huge lambda
  delta_good = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = J)
  Delta_good = do.call(cbind, lapply(1:J, function(j) diag(delta_good[j, ]) %*% t(delta_good)))
  llg0 = loglikelihood(c(c(Pi), c(Delta_good)), tab, J, K, lambda_dd = 0)
  llg = loglikelihood(c(c(Pi), c(Delta_good)), tab, J, K, lambda_dd = 1e6)
  expect_equal(llg, llg0)
})
