# EST-1: the trajectory EM. These tests mirror the simulation the estimator
# was validated on (gate G0, code/sim/m2 in the spurious-mobility repo): the
# same DGP, the same truths, and the same three claims the simulation made —
# recovery, the k = 2 versus k = 3 identification contrast, and that a
# missing measure must be handled as a missing emission.

# ---------------------------------------------------------------------------
# The DGP, generic in k. One anchor observed directly and once through a link
# of its own; k linked measures of a traced person whose latent state follows
# T; one latent rival shared across links with probability pi_s and a fresh
# rival otherwise; every observed value passed through (1-mu) I + mu K; each
# linked measure observed with probability p_obs.
# ---------------------------------------------------------------------------

traj_T = function(J, stay = 0.7, decay = 0.5, maxd = 2){
  Tm = matrix(0, J, J)
  for(j in seq_len(J)){
    nb = setdiff(max(1, j - maxd):min(J, j + maxd), j)
    w = decay^abs(nb - j)
    Tm[j, nb] = (1 - stay) * w / sum(w)
    Tm[j, j] = stay
  }
  Tm
}
traj_pow = function(M, h){
  R = diag(nrow(M))
  while(h > 0){ if(h %% 2 == 1) R = R %*% M; M = M %*% M; h = h %/% 2 }
  R
}
traj_draw = function(P, x){
  cp = t(apply(P, 1, cumsum))
  as.integer(rowSums(runif(length(x)) > cp[x, , drop = FALSE]) + 1L)
}
traj_local = function(J, q, pi0){
  L = matrix(0, J, J)
  for(i in seq_len(J)){ nb = max(1, i - 1):min(J, i + 1); L[i, nb] = q / length(nb) }
  L + (1 - q) * matrix(pi0, J, J, byrow = TRUE)
}
traj_kernel = function(J, nb = 2L){
  blk = as.integer(cut(seq_len(J), nb, labels = FALSE)); K = matrix(0, J, J)
  for(j in seq_len(J)){ idx = which(blk == blk[j]); K[j, idx] = 1 / length(idx) }
  K
}

traj_dgp = function(k, J = 6, N = 3e4, alpha = 0.3, s = 0.6, alpha_f = 0.25,
                    mu = 0.15, q_local = 0.7, dstep = rep(1L, k - 1L),
                    dstep_f = 1L, h1 = 4L, p_obs = rep(1, k), seed = 1){
  set.seed(seed)
  Tm = traj_T(J)
  w = rgamma(J, 2, 1); pi0 = w / sum(w)
  Th1 = traj_pow(Tm, h1)
  Kmat = traj_kernel(J)
  if(length(mu) == 1L) mu = rep(mu, k + 2L)
  Ms = lapply(mu, function(m) (1 - m) * diag(J) + m * Kmat)
  pis = s * alpha; b = if(pis < 1) (alpha - pis) / (1 - pis) else 0
  Tst = lapply(dstep, function(d) traj_pow(Tm, d))
  Tf = traj_pow(Tm, dstep_f)

  i = sample.int(J, N, TRUE, prob = pi0)
  xu = traj_draw(Tf, i)
  own = vector("list", k); own[[1]] = traj_draw(Th1, i)
  for(j in seq_len(k - 1L)) own[[j + 1L]] = traj_draw(Tst[[j]], own[[j]])

  L = traj_local(J, q_local, pi0)
  chain = function(x0){
    out = vector("list", k); out[[1]] = traj_draw(Th1, x0)
    for(j in seq_len(k - 1L)) out[[j + 1L]] = traj_draw(Tst[[j]], out[[j]])
    out
  }
  riv = chain(traj_draw(L, i))
  fresh = lapply(seq_len(k), function(j) chain(traj_draw(L, i))[[j]])

  shared = runif(N) < pis
  Z = lapply(seq_len(k), function(j){
    fail = shared | (runif(N) < b)
    ifelse(!fail, own[[j]], ifelse(shared, riv[[j]], fresh[[j]]))
  })
  Zu = ifelse(runif(N) < alpha_f, traj_draw(Tf, traj_draw(L, i)), xu)

  d = data.frame(X = traj_draw(Ms[[k + 1L]], i),
                 Xf = traj_draw(Ms[[k + 2L]], as.integer(Zu)))
  for(j in seq_len(k)){
    y = traj_draw(Ms[[j]], as.integer(Z[[j]]))
    d[[paste0("Y", j)]] = ifelse(runif(N) < p_obs[j], y, NA_integer_)
  }
  tab = aggregate(list(n = rep(1, N)),
                  by = lapply(d, function(v) addNA(factor(v, levels = 1:J))),
                  FUN = sum)
  for(nm in setdiff(names(tab), "n"))
    tab[[nm]] = suppressWarnings(as.integer(as.character(tab[[nm]])))

  G = vector("list", k); G[[1]] = L %*% Th1
  for(j in seq_len(k - 1L)) G[[j + 1L]] = G[[j]] %*% Tst[[j]]

  list(tab = tab, J = J, k = k, T = Tm, rho = G, rho_f = L %*% Tf,
       kernel = Kmat, alpha = alpha, s = s, alpha_f = alpha_f, mu = mu,
       pi_shared = pis, dstep = dstep,
       Pi = unclass(table(factor(i, 1:J), factor(own[[1]], 1:J))) / N)
}

traj_fit = function(d, ...){
  args = list(tab = d$tab, J = d$J, rho = d$rho, rho_f = d$rho_f,
              kernel = d$kernel, dstep = d$dstep, alpha_0 = 0.2,
              init = list(s = 0.3, mu = 0.1, alpha_f = 0.2),
              tol = 1e-10, maxit = 1200)
  extra = list(...)
  args[names(extra)] = extra
  suppressWarnings(do.call(misclassifyr_traj_em, args))
}
frob = function(A, B) sqrt(sum((A - B)^2))

# ---------------------------------------------------------------------------

test_that("the trajectory EM recovers alpha, s, mu, T and Pi", {
  d = traj_dgp(3, J = 8, N = 3e4, seed = 11)
  fit = traj_fit(d)
  expect_true(fit$converged)
  expect_equal(fit$k, 3L)
  # the slab
  expect_lt(max(abs(fit$alpha - d$alpha)), 0.06)
  expect_lt(abs(fit$s - d$s), 0.12)
  expect_lt(abs(fit$alpha_f - d$alpha_f), 0.05)
  # the measurement layer
  expect_lt(abs(fit$mu[1] - d$mu[1]), 0.03)
  # the dynamics and the estimand
  expect_lt(frob(fit$T, d$T), 0.10)
  expect_lt(max(abs(fit$Pi - d$Pi)), 0.01)
  expect_equal(sum(fit$Pi), 1)
  expect_equal(unname(rowSums(fit$T)), rep(1, d$J))
  # (alpha, s) and (pi_shared, b) are the same slab in two parameterisations
  expect_equal(unname(fit$alpha),
               unname(fit$pi_shared + (1 - fit$pi_shared) * fit$b))
  expect_equal(fit$s, unname(fit$pi_shared / fit$alpha[1]))
})

test_that("the (alpha, s) profile is flat at k = 2 and curved at k = 3", {
  # The design point of the three-link bundle: at k = 2 alpha and s trade off
  # along a ridge, so profiling pi_shared out moves the likelihood barely at
  # all; a third link separates a shared rival from a coincidence.
  d = traj_dgp(3, J = 6, N = 3e4, seed = 22)
  tab3 = d$tab
  tab2 = aggregate(list(n = tab3$n),
                   by = list(X = tab3$X, Xf = tab3$Xf,
                             Y1 = tab3$Y1, Y2 = tab3$Y2), FUN = sum)
  grid = c(0.02, 0.08, 0.14, 0.20)
  profile = function(tb, kk) vapply(grid, function(pp)
    suppressWarnings(misclassifyr_traj_em(
      tb, J = d$J, rho = d$rho[seq_len(kk)], rho_f = d$rho_f,
      kernel = d$kernel, alpha_0 = 0.25, pi_shared_fixed = pp,
      init = list(mu = 0.12, alpha_f = 0.2), tol = 1e-9,
      maxit = 500))$loglik, numeric(1))

  ll3 = profile(tab3, 3); ll2 = profile(tab2, 2)
  expect_gt(diff(range(ll3)), 3)      # k = 3: a bowl
  expect_lt(diff(range(ll2)), 2)      # k = 2: a ridge
  expect_gt(diff(range(ll3)) / diff(range(ll2)), 3)
})

test_that("a missing measure is handled as a missing emission", {
  d = traj_dgp(3, J = 6, N = 3e4, seed = 11, p_obs = c(1, 0.9, 0.7))
  expect_gt(mean(is.na(d$tab$Y3)), 0)
  fit = traj_fit(d)
  expect_true(fit$converged)
  # the units with a missing third measure still carry the branches their
  # observed pattern supports, so the estimates stay near the truth
  expect_lt(max(abs(fit$alpha - d$alpha)), 0.08)
  expect_lt(abs(fit$s - d$s), 0.15)
  expect_lt(abs(fit$mu[1] - d$mu[1]), 0.03)
  expect_lt(abs(fit$alpha_f - d$alpha_f), 0.05)
  expect_lt(frob(fit$T, d$T), 0.10)

  # NA and 0 are the same missing marker
  dz = d
  for(j in 1:3){ cc = paste0("Y", j); dz$tab[[cc]][is.na(dz$tab[[cc]])] = 0L }
  f_na = traj_fit(d, maxit = 12)
  f_zero = traj_fit(dz, maxit = 12)
  expect_equal(f_zero$loglik, f_na$loglik)
  expect_equal(f_zero$T, f_na$T)

  # dropping the incomplete units instead is a different (smaller) problem
  keep = !is.na(d$tab$Y3)
  expect_lt(sum(d$tab$n[keep]), sum(d$tab$n))
})

test_that("mu can vary by census year through mu_index", {
  d = traj_dgp(3, J = 6, N = 3e4, seed = 11,
               mu = c(0.10, 0.10, 0.10, 0.25, 0.25))
  fit = traj_fit(d, mu_index = c(1, 1, 1, 2, 2), init = list(
    s = 0.3, mu = 0.15, alpha_f = 0.2))
  expect_true(fit$converged)
  expect_equal(length(fit$mu), 2L)
  expect_lt(abs(fit$mu[1] - 0.10), 0.03)
  expect_lt(abs(fit$mu[2] - 0.25), 0.04)
})

test_that("the freeze options are no-ops when frozen at the free estimates", {
  d = traj_dgp(3, J = 6, N = 2e4, seed = 11)
  free = traj_fit(d)
  expect_true(free$converged)

  fT = traj_fit(d, T_fixed = free$T)
  expect_equal(fT$T, free$T)
  expect_lt(abs(fT$loglik - free$loglik), 0.05)
  expect_lt(abs(fT$mu[1] - free$mu[1]), 0.01)

  fA = traj_fit(d, alpha_fixed = unname(free$alpha), s_fixed = free$s)
  expect_equal(unname(fA$alpha), unname(free$alpha), tolerance = 1e-6)
  expect_equal(fA$s, free$s, tolerance = 1e-6)
  expect_lt(abs(fA$loglik - free$loglik), 0.05)

  fM = traj_fit(d, mu_fixed = free$mu)
  expect_equal(fM$mu, free$mu)
  expect_lt(abs(fM$loglik - free$loglik), 0.05)

  fP = traj_fit(d, Pi_fixed = free$Pi)
  expect_equal(fP$Pi, free$Pi)
  expect_lt(abs(fP$loglik - free$loglik), 0.05)

  # a frozen alpha really is frozen, wherever it is put
  fF = traj_fit(d, alpha_fixed = rep(0.4, 3), s_fixed = 0.5, maxit = 30)
  expect_equal(unname(fF$alpha), rep(0.4, 3), tolerance = 1e-8)
  expect_equal(fF$s, 0.5, tolerance = 1e-8)
})

test_that("plain EM increases the log likelihood monotonically", {
  d = traj_dgp(3, J = 6, N = 2e4, seed = 11)
  fit = traj_fit(d, squarem = FALSE, tol = 1e-8, maxit = 60)
  expect_true(all(diff(fit$loglik_trace) > -1e-8))
})

test_that("the estimator is generic in k and in the horizon spacing", {
  f2 = traj_fit(traj_dgp(2, J = 6, N = 2e4, seed = 11))
  expect_equal(f2$k, 2L)
  expect_equal(length(f2$alpha), 2L)
  expect_true(f2$converged)

  d4 = traj_dgp(4, J = 5, N = 1.5e4, h1 = 3L, seed = 11)
  f4 = traj_fit(d4)
  expect_equal(f4$k, 4L)
  expect_lt(max(abs(f4$alpha - d4$alpha)), 0.10)
  expect_lt(abs(f4$mu[1] - d4$mu[1]), 0.04)

  # two measures at the same date, then a two-step gap
  dg = traj_dgp(3, J = 5, N = 1.5e4, h1 = 3L, seed = 11, dstep = c(0L, 2L))
  fg = traj_fit(dg)
  expect_lt(abs(fg$mu[1] - dg$mu[1]), 0.05)
  expect_equal(unname(rowSums(fg$T)), rep(1, dg$J))
})

test_that("the measurement layer and the anchor link are optional", {
  d = traj_dgp(3, J = 6, N = 2e4, seed = 11, mu = 0)
  f0 = suppressWarnings(misclassifyr_traj_em(
    d$tab, J = d$J, rho = d$rho, rho_f = d$rho_f, kernel = NULL,
    alpha_0 = 0.2, init = list(s = 0.3, alpha_f = 0.2), tol = 1e-10,
    maxit = 1200))
  expect_equal(unname(f0$mu), 0)
  expect_true(f0$converged)
  expect_lt(abs(f0$alpha_f - d$alpha_f), 0.06)

  no_anchor_link = aggregate(list(n = d$tab$n),
                             by = list(X = d$tab$X, Y1 = d$tab$Y1,
                                       Y2 = d$tab$Y2, Y3 = d$tab$Y3), FUN = sum)
  fx = suppressWarnings(misclassifyr_traj_em(
    no_anchor_link, J = d$J, rho = d$rho, kernel = NULL, alpha_0 = 0.2,
    init = list(s = 0.3), tol = 1e-10, maxit = 1200))
  expect_true(is.na(fx$alpha_f))
  expect_lt(max(abs(fx$alpha - f0$alpha)), 0.05)
})

test_that("chunking changes nothing but working memory", {
  d = traj_dgp(3, J = 5, N = 8e3, h1 = 3L, seed = 5, p_obs = c(1, 0.85, 0.7))
  whole = traj_fit(d, maxit = 12, chunk_size = 1e6)
  chunked = traj_fit(d, maxit = 12, chunk_size = 100L)
  expect_equal(chunked$loglik, whole$loglik, tolerance = 1e-8)
  expect_equal(chunked$T, whole$T, tolerance = 1e-8)
  expect_equal(chunked$Pi, whole$Pi, tolerance = 1e-8)
})

test_that("malformed inputs are rejected", {
  d = traj_dgp(3, J = 5, N = 2e3, h1 = 3L, seed = 5)
  tb = d$tab

  expect_error(misclassifyr_traj_em(tb[, c("X", "Y1", "n")], J = 5,
                                    rho = d$rho[[1]]),
               "at least two linked-measure columns")
  expect_error(misclassifyr_traj_em(tb[, setdiff(names(tb), "X")], J = 5,
                                    rho = d$rho[[1]]),
               "column `X`")
  bad = tb; bad$Y1[1] = 99L
  expect_error(misclassifyr_traj_em(bad, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f),
               "integer codes in 1..J")
  expect_error(misclassifyr_traj_em(rbind(tb[1, ], tb), J = 5,
                                    rho = d$rho[[1]], rho_f = d$rho_f),
               "duplicated measure patterns")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = matrix(0.5, 5, 5),
                                    rho_f = d$rho_f),
               "rows summing to one")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]]),
               "`rho_f` is required")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f, kernel = matrix(1, 5, 5)),
               "rows summing to one")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f, dstep = c(1L, 1L, 1L)),
               "length `k - 1`")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f, s_fixed = 0.5),
               "must be supplied together")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f, alpha_fixed = rep(0.3, 3),
                                    s_fixed = 0.5, pi_shared_fixed = 0.1),
               "not both")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f, mu_index = c(1, 3, 1, 1, 1),
                                    kernel = d$kernel),
               "consecutive group labels")
  expect_error(misclassifyr_traj_em(tb, J = 5, rho = d$rho[[1]],
                                    rho_f = d$rho_f,
                                    T_fixed = matrix(1, 5, 5)),
               "rows summing to one")
})
