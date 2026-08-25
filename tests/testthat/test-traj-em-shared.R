test_that("shared trajectory EM recovers (alpha, s) across two views at k = 3", {
  set.seed(41)
  N <- 12000; k <- 3; J1 <- 12; J2 <- 4
  pis <- 0.06; b <- c(0.08, 0.10, 0.12)
  x1 <- sample.int(J1, N, TRUE); x2 <- (x1 - 1) %/% 3 + 1
  shared <- runif(N) < pis
  r1 <- sample.int(J1, N, TRUE)
  Y <- matrix(0L, N, k)
  for (j in 1:k) {
    fresh <- !shared & runif(N) < b[j]
    Y[, j] <- ifelse(shared, r1, ifelse(fresh, sample.int(J1, N, TRUE), x1))
  }
  tab1 <- aggregate(list(n = rep(1, N)),
    by = list(X = x1, Y1 = Y[, 1], Y2 = Y[, 2], Y3 = Y[, 3]), FUN = sum)
  cg <- function(v) (v - 1L) %/% 3L + 1L
  tab2 <- aggregate(list(n = rep(1, N)),
    by = list(X = x2, Y1 = cg(Y[, 1]), Y2 = cg(Y[, 2]), Y3 = cg(Y[, 3])),
    FUN = sum)
  rho1 <- lapply(1:k, function(j) matrix(1 / J1, J1, J1))
  rho2 <- lapply(1:k, function(j) matrix(1 / J2, J2, J2))

  fit <- misclassifyr_traj_em_shared(list(
    fine   = list(tab = tab1, J = J1, rho = rho1),
    coarse = list(tab = tab2, J = J2, rho = rho2)), maxit = 400)

  truth_a <- pis + (1 - pis) * b
  expect_true(fit$converged)
  expect_equal(unname(fit$alpha), truth_a, tolerance = 0.12)
  expect_equal(fit$s, pis / truth_a[1], tolerance = 0.15)
  # composite trace is monotone (block GEM)
  expect_true(all(diff(fit$loglik_trace) > -1e-6))
  # per-view Pi on each view's own alphabet
  expect_equal(dim(fit$views$fine$Pi), c(J1, J1))
  expect_equal(dim(fit$views$coarse$Pi), c(J2, J2))
})

test_that("views with different k are refused, and owned args are refused", {
  tab3 <- data.frame(X = 1, Y1 = 1, Y2 = 1, Y3 = 1, n = 5)
  tab2 <- data.frame(X = 1, Y1 = 1, Y2 = 1, n = 5)
  rho3 <- lapply(1:3, function(j) matrix(1, 1, 1))
  rho2 <- lapply(1:2, function(j) matrix(1, 1, 1))
  expect_error(
    misclassifyr_traj_em_shared(list(
      a = list(tab = tab3, J = 1, rho = rho3),
      b = list(tab = tab2, J = 1, rho = rho2))),
    "different numbers of linked measures")
  expect_error(
    misclassifyr_traj_em_shared(list(
      a = list(tab = tab3, J = 1, rho = rho3, s_fixed = 0.2))),
    "driver owns")
})

test_that("traj_em exposes E-step aggregates consistent with its alpha", {
  set.seed(5)
  N <- 4000; J <- 6
  x <- sample.int(J, N, TRUE)
  y1 <- ifelse(runif(N) < 0.2, sample.int(J, N, TRUE), x)
  y2 <- ifelse(runif(N) < 0.2, sample.int(J, N, TRUE), x)
  tab <- aggregate(list(n = rep(1, N)),
                   by = list(X = x, Y1 = y1, Y2 = y2), FUN = sum)
  rho <- lapply(1:2, function(j) matrix(1 / J, J, J))
  fit <- suppressWarnings(
    misclassifyr_traj_em(tab, J = J, rho = rho, squarem = FALSE, maxit = 200))
  expect_named(fit$Ecounts, c("R_C", "R_SSS", "R_f", "N", "ll"),
               ignore.order = TRUE)
  expect_equal(fit$Ecounts$N, N)
  # at a fixed point the count ratios reproduce the returned parameters
  fail1 <- sum(fit$Ecounts$R_C[bitwAnd(0:3, 1L) == 0L])
  b1 <- fail1 / (N - fit$Ecounts$R_SSS)
  a1_implied <- fit$pi_shared + (1 - fit$pi_shared) * b1
  expect_equal(unname(fit$alpha[1]), a1_implied, tolerance = 1e-4)
})
