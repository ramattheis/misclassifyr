test_that("shared-alpha EM recovers alpha across heterogeneous designs", {
  set.seed(93)
  N <- 12000; alpha <- c(0.2, 0.3)
  xf <- sample.int(12, N, replace = TRUE)
  xc <- (xf - 1) %/% 3 + 1
  fail1 <- runif(N) < alpha[1]; fail2 <- runif(N) < alpha[2]
  rf1 <- sample.int(12, N, TRUE); rf2 <- sample.int(12, N, TRUE)
  f1 <- ifelse(fail1, rf1, xf); f2 <- ifelse(fail2, rf2, xf)
  c1 <- ifelse(fail1, (rf1 - 1) %/% 3 + 1, xc)
  c2 <- ifelse(fail2, (rf2 - 1) %/% 3 + 1, xc)
  tab_f <- aggregate(list(n = rep(1, N)),
                     by = list(X = xf, Y1 = f1, Y2 = f2), FUN = sum)
  tab_c <- aggregate(list(n = rep(1, N)),
                     by = list(X = xc, Y1 = c1, Y2 = c2), FUN = sum)

  fit <- misclassifyr_rl_em_shared(list(
    fine   = list(tab = tab_f, J = 12, K = 12),
    coarse = list(tab = tab_c, J = 4,  K = 4)))

  expect_true(fit$converged)
  expect_equal(unname(fit$alpha[1]), 0.2, tolerance = 0.03)
  expect_equal(unname(fit$alpha[2]), 0.3, tolerance = 0.03)
  expect_named(fit$designs, c("fine", "coarse"))
  # per-design Pi returned on each design's own alphabet
  expect_lte(max(fit$designs$coarse$Pi$j), 4)
  expect_lte(max(fit$designs$fine$Pi$j), 12)
  # composite log likelihood is monotone
  expect_true(all(diff(fit$loglik_trace) > -1e-6))
})

test_that("alpha_fixed is honored and designs keep their own structure", {
  set.seed(11)
  N <- 6000
  x <- sample.int(6, N, TRUE)
  f1 <- runif(N) < 0.25; f2 <- runif(N) < 0.25
  y1 <- ifelse(f1, sample.int(6, N, TRUE), x)
  y2 <- ifelse(f2, sample.int(6, N, TRUE), x)
  tab <- aggregate(list(n = rep(1, N)),
                   by = list(X = x, Y1 = y1, Y2 = y2), FUN = sum)
  fit <- misclassifyr_rl_em_shared(
    list(a = list(tab = tab, J = 6, K = 6),
         b = list(tab = tab, J = 6, K = 6)),
    alpha_fixed = c(0.25, 0.25))
  expect_equal(unname(fit$alpha), c(0.25, 0.25))
})

test_that("a design-specific T2 operator is accepted and applied", {
  set.seed(29)
  N <- 20000; J <- 5; alpha <- c(0.2, 0.2)
  # true change between the two readings: sticky operator
  Tm <- matrix(0.05, J, J); diag(Tm) <- 0.8
  x  <- sample.int(J, N, TRUE)
  ystar2 <- vapply(x, function(j) sample.int(J, 1, prob = Tm[j, ]), 1L)
  f1 <- runif(N) < alpha[1]; f2 <- runif(N) < alpha[2]
  y1 <- ifelse(f1, sample.int(J, N, TRUE), x)
  y2 <- ifelse(f2, sample.int(J, N, TRUE), ystar2)
  tab <- aggregate(list(n = rep(1, N)),
                   by = list(X = x, Y1 = y1, Y2 = y2), FUN = sum)
  T2 <- data.frame(j = rep(1:J, each = J), l = rep(1:J, J),
                   t = as.vector(t(Tm)))
  fit <- misclassifyr_rl_em_shared(list(
    withT = list(tab = tab, J = J, K = J, T2 = T2)))
  expect_true(fit$converged)
  expect_equal(unname(fit$alpha[1]), 0.2, tolerance = 0.06)
  expect_equal(unname(fit$alpha[2]), 0.2, tolerance = 0.06)
})
