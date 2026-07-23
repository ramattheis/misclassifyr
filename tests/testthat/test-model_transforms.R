# Exact unit tests for the parameter -> (Pi, Delta) transforms.

test_that("model_to_Pi_NP applies the logit link exactly", {
  # phi = log(2, 3, 4): probabilities (2, 3, 4, 1) / 10
  expect_equal(model_to_Pi_NP(log(c(2, 3, 4)), J = 2), c(0.2, 0.3, 0.4, 0.1))
  # Flat phi maps to the uniform distribution
  expect_equal(model_to_Pi_NP(rep(0, 5), J = 2), rep(1 / 6, 6))
  # Always lands on the simplex
  set.seed(1)
  p = model_to_Pi_NP(rnorm(8), J = 3)
  expect_equal(sum(p), 1)
  expect_true(all(p > 0))
  expect_length(p, 9)
})

test_that("model_to_Delta_NP_ind builds a valid conditional distribution", {
  J = 3
  psi = rep(0, 2 * J * (J - 1))
  D = matrix(model_to_Delta_NP_ind(psi), nrow = J)

  # J x J^2 matrix; each row (given Y*) is a distribution over (Y1, Y2)
  expect_equal(dim(D), c(J, J^2))
  expect_equal(unname(rowSums(D)), rep(1, J))
  expect_true(all(D >= 0))
  # Flat psi gives the uniform joint distribution
  expect_equal(c(D), rep(1 / J^2, J^3))
})

test_that("model_to_Delta_NP_ind factorizes as Pr(Y1|Y*) * Pr(Y2|Y*)", {
  # Hand-computed J = 2 case: column blocks indexed by Y2, within-block
  # columns by Y1, rows by Y*.
  psi = c(1, 0.5, -0.3, 0.2)
  D = matrix(model_to_Delta_NP_ind(psi), nrow = 2)
  d1 = rbind(plogis(psi[1:2]), 1 - plogis(psi[1:2]))  # Pr(Y1 = r | Y* = c)
  d2 = rbind(plogis(psi[3:4]), 1 - plogis(psi[3:4]))  # Pr(Y2 = r | Y* = c)
  for (ys in 1:2) for (y1 in 1:2) for (y2 in 1:2) {
    expect_equal(D[ys, (y2 - 1) * 2 + y1], d1[y1, ys] * d2[y2, ys])
  }
})

test_that("model_to_Delta_RL_ind builds the record-linkage structure", {
  J = 3
  a = 0.2
  # Uniform row margins, common column scale alpha = 0.2 for both measures
  psi = c(rep(0, J - 1), rep(0, J - 1), rep(qlogis(a), J), rep(qlogis(a), J))
  D = matrix(model_to_Delta_RL_ind(psi), nrow = J)

  expect_equal(dim(D), c(J, J^2))
  expect_equal(unname(rowSums(D)), rep(1, J))

  # Single-measure matrix: (1 - alpha) on the diagonal + alpha * uniform row
  d_single = diag(J) * (1 - a) + outer(rep(1 / J, J), rep(a, J))
  for (ys in 1:J) for (y1 in 1:J) for (y2 in 1:J) {
    expect_equal(D[ys, (y2 - 1) * J + y1], d_single[y1, ys] * d_single[y2, ys])
  }
})
