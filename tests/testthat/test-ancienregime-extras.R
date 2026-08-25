# The companion datasets and the recoveries the tour vignette shows.
# Subsampled where estimation is involved, so the suite stays fast; the
# tolerances are loose accordingly. The full-sample versions are in
# vignette("tour-ancien-regime").

test_that("the companion datasets load with the documented shape", {
  env <- new.env()
  data("ancienregime_parishes", package = "misclassifyr", envir = env)
  data("ancienregime_lineages", package = "misclassifyr", envir = env)
  data("ancienregime_truth", package = "misclassifyr", envir = env)
  ancienregime_parishes <- env$ancienregime_parishes
  ancienregime_lineages <- env$ancienregime_lineages
  ancienregime_truth <- env$ancienregime_truth

  expect_equal(nrow(ancienregime_parishes), 120000L)
  expect_named(ancienregime_parishes,
               c("province_birth", "parish_1750", "parish_1770_linked",
                 "parish_1780_linked"))
  expect_true(all(ancienregime_parishes$parish_1750 %in% 1:240))

  expect_equal(nrow(ancienregime_lineages), 100000L)
  expect_true(all(ancienregime_lineages$son_occupation_1790_linked %in%
                    ancienregime_truth$occupations))

  expect_equal(nrow(ancienregime_truth$parish_register), 240L)
  expect_equal(ancienregime_truth$lineages$alpha,
               ancienregime_truth$lineages$pi_shared +
                 (1 - ancienregime_truth$lineages$pi_shared) *
                 ancienregime_truth$lineages$b)
  # the vingtieme flow table is row-stochastic within origin
  s <- tapply(ancienregime_truth$vingtieme_flows$t,
              ancienregime_truth$vingtieme_flows$j, sum)
  expect_true(all(abs(s - 1) < 1e-8))
})

test_that("local rivals: stacked recovers alpha, unconditional collapses", {
  skip_on_cran()
  env <- new.env()
  data("ancienregime_parishes", package = "misclassifyr", envir = env)
  data("ancienregime_truth", package = "misclassifyr", envir = env)
  d <- env$ancienregime_parishes
  ancienregime_truth <- env$ancienregime_truth
  ving <- ancienregime_truth$vingtieme_flows[, c("j", "l", "t")]

  tab <- aggregate(list(n = rep(1L, nrow(d))),
                   by = list(X = d$parish_1750, Y1 = d$parish_1770_linked,
                             Y2 = d$parish_1780_linked), FUN = sum)
  fit_u <- misclassifyr_rl_em(tab, J = 240, K = 240, T2 = ving,
                              maxit = 300, tol = 1e-8)
  tabs <- lapply(split(d, d$province_birth), function(z)
    aggregate(list(n = rep(1L, nrow(z))),
              by = list(X = z$parish_1750, Y1 = z$parish_1770_linked,
                        Y2 = z$parish_1780_linked), FUN = sum))
  fit_s <- misclassifyr_rl_em_stacked(tabs, J = 240, K = 240, T2 = ving,
                                      maxit = 300, tol = 1e-8)

  truth <- ancienregime_truth$parishes$alpha
  # the unconditional fit understates severely; the stacked fit recovers
  expect_lt(fit_u$alpha[1], truth[1] / 2)
  expect_lt(abs(fit_s$alpha[1] - truth[1]), 0.02)
  expect_lt(abs(fit_s$alpha[2] - truth[2]), 0.03)
})

test_that("the trajectory model separates alpha, s, mu on a subsample", {
  skip_on_cran()
  env <- new.env()
  data("ancienregime_lineages", package = "misclassifyr", envir = env)
  data("ancienregime_truth", package = "misclassifyr", envir = env)
  ancienregime_lineages <- env$ancienregime_lineages
  tr <- env$ancienregime_truth$lineages
  occ <- env$ancienregime_truth$occupations
  set.seed(93)
  lin <- ancienregime_lineages[sample.int(100000L, 40000L), ]
  code <- function(v) match(v, occ)
  tab <- aggregate(list(n = rep(1L, nrow(lin))),
    by = list(Y1 = code(lin$son_occupation_1770_linked),
              Y2 = code(lin$son_occupation_1780_linked),
              Y3 = code(lin$son_occupation_1790_linked),
              X  = code(lin$father_occupation_1750),
              Xf = code(lin$father_occupation_1745_linked)), FUN = sum)
  fit <- misclassifyr_traj_em(tab, J = 6, rho = tr$occupation_slab,
                              rho_f = tr$father_margin, kernel = tr$kernel,
                              dstep = c(1, 1), dstep_f = 1, maxit = 600)
  expect_true(fit$converged)
  expect_lt(max(abs(fit$alpha - tr$alpha)), 0.03)
  expect_lt(abs(fit$mu[1] - tr$mu), 0.02)
  expect_lt(abs(fit$alpha_f - tr$alpha_f), 0.03)
  expect_lt(abs(fit$s - tr$s), 0.12)
  # and the estimand: the fitted joint sits on the generating joint
  Pi_true <- tr$Pi * tr$father_margin
  expect_lt(max(abs(fit$Pi / sum(fit$Pi) - Pi_true)), 0.01)
})
