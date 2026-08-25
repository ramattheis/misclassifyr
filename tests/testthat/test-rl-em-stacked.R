# Stacked shared-alpha EM: conditional rival draws (the Mattheis
# critique): if failed links draw from the candidate pool matching the
# linking keys (e.g. same birthplace), an unconditional slab is too
# diffuse - rivals look like correct links and pooled alpha-hat is
# biased down. Cell-conditional estimation with shared alpha fixes it.

test_that("pooled EM understates alpha under conditional rivals; stacked recovers it", {
  set.seed(41)
  C = 6; J = 30; N = 3e5; a1 = 0.25; a2 = 0.25
  cell = sample.int(C, N, replace = TRUE)
  base = (cell - 1) * 5
  xi = base + sample.int(5, N, replace = TRUE)
  j  = ifelse(runif(N) < 0.8, xi, base + sample.int(5, N, replace = TRUE))
  # Rival draws are CELL-CONDITIONAL (same 5 home counties)
  y1 = ifelse(runif(N) < a1, base + sample.int(5, N, replace = TRUE), j)
  y2 = ifelse(runif(N) < a2, base + sample.int(5, N, replace = TRUE), j)
  df = data.frame(cell = cell, X = xi, Y1 = y1, Y2 = y2)

  pooled_tab = dplyr::count(df, X, Y1, Y2, name = "n")
  out_pooled = misclassifyr_rl_em(pooled_tab, J = J, K = J)

  tabs = lapply(split(df, df$cell), function(d)
    dplyr::count(d, X, Y1, Y2, name = "n"))
  out_stacked = misclassifyr_rl_em_stacked(tabs, J = J, K = J)

  expect_true(out_stacked$converged)
  # Stacked, cell-conditional: recovers the true alpha
  expect_lt(abs(out_stacked$alpha["alpha1"] - a1), 0.02)
  expect_lt(abs(out_stacked$alpha["alpha2"] - a2), 0.02)
  # Pooled, unconditional slab: materially biased downward
  expect_lt(unname(out_pooled$alpha["alpha1"]), a1 - 0.05)
  expect_lt(unname(out_pooled$alpha["alpha2"]), a2 - 0.05)
})

test_that("stacked EM agrees with single-cell EM when there is one cell", {
  set.seed(42)
  J = 5; N = 5e4
  xi = sample.int(J, N, replace = TRUE)
  j  = ifelse(runif(N) < 0.8, xi, sample.int(J, N, replace = TRUE))
  rho = tabulate(j, J) / N
  y1 = ifelse(runif(N) < 0.2, sample.int(J, N, replace = TRUE, prob = rho), j)
  y2 = ifelse(runif(N) < 0.2, sample.int(J, N, replace = TRUE, prob = rho), j)
  tab = dplyr::count(data.frame(X = xi, Y1 = y1, Y2 = y2), X, Y1, Y2, name = "n")
  a = misclassifyr_rl_em(tab, J, J)
  b = misclassifyr_rl_em_stacked(list(tab), J, J)
  expect_equal(unname(a$alpha), unname(b$alpha), tolerance = 1e-6)
})
