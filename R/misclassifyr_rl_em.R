#' EM estimation of the record-linkage misclassification model at scale
#'
#' Estimates the common-alpha record-linkage model for high-dimensional
#' discrete outcomes (e.g. counties, \eqn{J} in the thousands), where the
#' general-purpose \code{misclassifyr()} machinery is infeasible: both noisy
#' measures follow \eqn{\Delta^{(m)} = (1-\alpha_m) I + \alpha_m \mathbf{1}
#' \rho_m^{\top}} — a link is correct with probability \eqn{1-\alpha_m},
#' otherwise the observed value is an independent draw from \eqn{\rho_m}.
#' The mixture likelihood factorizes into four terms per observed cell, so
#' each EM iteration is O(observed cells + support of Pi) and neither the
#' misclassification matrices nor the balanced tabulation are ever formed.
#'
#' The E-step computes, for each observed cell \eqn{(X=i, Y_1=k, Y_2=l)},
#' the posterior weights of the four linkage configurations
#' (both links correct — requires \eqn{k=l}; only the first correct; only
#' the second; neither). The M-step has closed forms: \eqn{\alpha_m} is the
#' expected share of failed links, \eqn{\rho_m} the expected distribution of
#' failed-link draws, and \eqn{\Pi} the expected latent-cell counts (cells
#' where neither link is correct spread their mass over \eqn{\Pi_{\cdot,i}}).
#' \eqn{\Pi}'s support is fixed at initialization to the union of observed
#' \eqn{(Y_1, X)} and \eqn{(Y_2, X)} pairs (EM keeps zero cells at zero).
#'
#' @param tab A data.frame with columns `X`, `Y1`, `Y2` (positive integer
#'   codes; `Y1`, `Y2` share the same coding with `J` categories) and `n`
#'   (non-negative counts). Only observed cells need be present.
#' @param J An integer, the number of categories of the latent outcome.
#' @param K An integer, the number of categories of the regressor.
#' @param rho1,rho2 Optional length-`J` probability vectors fixing the
#'   failed-link draw distributions (e.g. known population margins). If
#'   `NULL` (default) they are estimated.
#' @param alpha_fixed Optional length-2 vector fixing (alpha1, alpha2)
#'   instead of estimating them — the second step of the two-step
#'   architecture for high-dimensional outcomes, where alpha is first
#'   estimated at a coarser aggregation (dense cells, standard asymptotics)
#'   and the fine-geography Pi is then estimated with alpha held fixed to
#'   avoid the incidental-parameters bias of a free high-dimensional Pi.
#' @param alpha_0 Numeric starting value for both alpha parameters.
#' @param T2 Optional sparse "true transition" operator for the SECOND
#'   measure, as a data.frame with columns `j` (latent category), `l`
#'   (observed category), `t` (probability), rows summing to one within
#'   `j`. When supplied, the second measure's model becomes
#'   \eqn{\Delta^{(2)} = (1-\alpha_2) T + \alpha_2 \mathbf{1}\rho_2^\top}:
#'   a correct link observes a draw from `T` applied to the latent value
#'   (e.g. true migration between the two observation dates, disciplined
#'   by external data such as the 1940 census migration question), while a
#'   failed link draws from \eqn{\rho_2}. `T2 = NULL` (default) is the
#'   identity — the plain record-linkage model. `T2` is held fixed (a
#'   plug-in), not estimated.
#' @param tol Relative log-likelihood convergence tolerance.
#' @param maxit Maximum EM iterations.
#' @param verbose Print the log likelihood every 25 iterations.
#' @return A list: `alpha` (length 2), `rho1`, `rho2`, `Pi` (data.frame
#'   `j`, `i`, `p` over the support), `loglik` (final), `loglik_trace`,
#'   `n_iter`, `converged`.
#' @seealso [misclassifyr_rl_em_stacked()] to share `alpha` across
#'   conditioning cells, and `vignette("sparse-and-large")`.
#' @examples
#' # Simulate the record-linkage model directly: the latent value moves away
#' # from X with probability 0.2, and each measure is a failed link (an
#' # independent draw from rho) with probability alpha_m.
#' set.seed(1)
#' J <- 40; N <- 50000
#' xi <- sample.int(J, N, replace = TRUE)
#' j  <- ifelse(runif(N) < 0.8, xi, sample.int(J, N, replace = TRUE))
#' rho <- tabulate(j, J) / N
#' y1 <- ifelse(runif(N) < 0.15, sample.int(J, N, replace = TRUE, prob = rho), j)
#' y2 <- ifelse(runif(N) < 0.30, sample.int(J, N, replace = TRUE, prob = rho), j)
#' tab <- aggregate(list(n = rep(1, N)),
#'                  by = list(X = xi, Y1 = y1, Y2 = y2), FUN = sum)
#'
#' fit <- misclassifyr_rl_em(tab, J = J, K = J)
#' fit$alpha        # close to (0.15, 0.30)
#' fit$converged
#' head(fit$Pi)     # only the cells on the support are carried
#'
#' # Two-step estimation: hold alpha at a value estimated elsewhere
#' fit2 <- misclassifyr_rl_em(tab, J = J, K = J, alpha_fixed = c(0.15, 0.30))
#' fit2$alpha
#' @export
misclassifyr_rl_em = function(tab, J, K,
                              rho1 = NULL, rho2 = NULL,
                              alpha_fixed = NULL,
                              alpha_0 = 0.2,
                              T2 = NULL,
                              tol = 1e-8, maxit = 500, verbose = FALSE){

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  if(!all(c("X","Y1","Y2","n") %in% colnames(tab))){
    stop("`tab` should have columns `X`, `Y1`, `Y2`, and `n`.")
  }
  tab = tab[tab$n > 0, , drop = FALSE]
  x = as.integer(tab$X); k = as.integer(tab$Y1); l = as.integer(tab$Y2)
  n = as.numeric(tab$n)
  if(any(is.na(x)) || any(x < 1) || any(x > K)) stop("`X` must be integer codes in 1..K.")
  if(any(is.na(k)) || any(k < 1) || any(k > J) ||
     any(is.na(l)) || any(l < 1) || any(l > J)) stop("`Y1`,`Y2` must be integer codes in 1..J.")
  if(any(duplicated(cbind(x, k, l)))) stop("`tab` has duplicated (X, Y1, Y2) cells.")
  for(rr in list(rho1, rho2)){
    if(!is.null(rr) && (length(rr) != J || any(rr < 0) || abs(sum(rr) - 1) > 1e-8)){
      stop("`rho1`/`rho2` must be length-J probability vectors.")
    }
  }
  N = sum(n)

  #------------------------------------------------------------
  # Support of Pi and precomputed indices (fixed across iterations)
  #------------------------------------------------------------

  # Support = union of observed (Y1, X) and (Y2, X) pairs, encoded j + J*(i-1)
  s_code = unique(c(k + J*(x - 1), l + J*(x - 1)))
  s_code = sort(s_code)
  S_j = ((s_code - 1) %% J) + 1
  S_i = ((s_code - 1) %/% J) + 1

  # Per-cell positions of (Y1, X) and (Y2, X) in the support vector
  a_ki = match(k + J*(x - 1), s_code)
  a_li = match(l + J*(x - 1), s_code)

  # Group indices for the X-margin of the support and of the cells
  ui   = sort(unique(S_i))
  gS_i = match(S_i, ui)     # support row -> position in ui
  gc_i = match(x,   ui)     # cell -> position in ui (every observed X is in ui)

  same = as.numeric(k == l)

  # Optional true-transition operator for measure 2
  use_T = !is.null(T2)
  if(use_T){
    if(!all(c("j","l","t") %in% colnames(T2))) stop("`T2` needs columns `j`, `l`, `t`.")
    Tj = as.integer(T2$j); Tl = as.integer(T2$l); Tt = as.numeric(T2$t)
    if(any(Tj < 1) || any(Tj > J) || any(Tl < 1) || any(Tl > J)) stop("`T2` codes must be in 1..J.")
    rs = rowsum(Tt, Tj)
    if(any(abs(rs - 1) > 1e-6)) stop("`T2` rows must sum to one within `j`.")
    Tmat = Matrix::sparseMatrix(i = Tj, j = Tl, x = Tt, dims = c(J, J))
    # Per-cell lookup T_{l|k} for the both-links-correct term (fixed):
    t_code = Tj + as.numeric(J)*(Tl - 1)
    t_kl   = Tt[match(k + as.numeric(J)*(l - 1), t_code)]
    t_kl[is.na(t_kl)] = 0
    # Cell codes for extracting m_{l,i} = (T' Pi)_{l,i} each iteration
    cell_li_code = l + as.numeric(J)*(x - 1)
  } else {
    t_kl = same          # T = I: T_{l|k} = 1{k = l}
  }

  #------------------------------------------------------------
  # Initial values
  #------------------------------------------------------------

  # Pi_0: empirical mass of observed (Y1, X) and (Y2, X) cells, averaged
  p = numeric(length(s_code))
  p_tmp = rowsum(c(n, n) / (2 * N), c(a_ki, a_li))
  p[as.integer(rownames(p_tmp))] = p_tmp

  est_rho1 = is.null(rho1); est_rho2 = is.null(rho2)
  if(est_rho1){ r1_tmp = rowsum(n, k); rho1 = numeric(J); rho1[as.integer(rownames(r1_tmp))] = r1_tmp / N }
  if(est_rho2){ r2_tmp = rowsum(n, l); rho2 = numeric(J); rho2[as.integer(rownames(r2_tmp))] = r2_tmp / N }
  est_alpha = is.null(alpha_fixed)
  if(!est_alpha){
    if(length(alpha_fixed) != 2 || any(alpha_fixed <= 0) || any(alpha_fixed >= 1)){
      stop("`alpha_fixed` must be a length-2 vector with values in (0, 1).")
    }
    alpha1 = alpha_fixed[1]; alpha2 = alpha_fixed[2]
  } else {
    alpha1 = alpha_0; alpha2 = alpha_0
  }

  #------------------------------------------------------------
  # EM iterations
  #------------------------------------------------------------

  loglik_trace = numeric(0)
  ll_old = -Inf
  converged = FALSE

  for(iter in seq_len(maxit)){

    # E-step ---------------------------------------------------
    pplus  = as.numeric(rowsum(p, gS_i))       # Pi_{+,i} over ui
    pk = p[a_ki]
    if(use_T){
      # m_{l,i} = (T' Pi)_{l,i}: the density of observing Y2 = l via a
      # CORRECT link when the latent is distributed as Pi_{.,i}
      P  = Matrix::sparseMatrix(i = S_j, j = S_i, x = p, dims = c(J, K))
      ms = Matrix::summary(Matrix::crossprod(Tmat, P))
      mv = ms$x[match(cell_li_code, ms$i + as.numeric(J)*(ms$j - 1))]
      mv[is.na(mv)] = 0
    } else {
      mv = p[a_li]
    }
    w11 = (1 - alpha1) * (1 - alpha2) * pk * t_kl
    w10 = (1 - alpha1) * alpha2 * rho2[l] * pk
    w01 = alpha1 * rho1[k] * (1 - alpha2) * mv
    w00 = alpha1 * alpha2 * rho1[k] * rho2[l] * pplus[gc_i]
    pc  = w11 + w10 + w01 + w00
    r11 = w11 / pc; r10 = w10 / pc; r01 = w01 / pc; r00 = w00 / pc

    ll = sum(n * log(pc))
    loglik_trace = c(loglik_trace, ll)
    if(verbose && iter %% 25 == 0) cat("iter", iter, "loglik", ll, "\n")
    if(is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + 1)){
      converged = TRUE
      break
    }
    ll_old = ll

    # M-step ---------------------------------------------------
    if(est_alpha){
      alpha1 = sum(n * (r01 + r00)) / N
      alpha2 = sum(n * (r10 + r00)) / N
    }
    if(est_rho1){
      r1_tmp = rowsum(n * (r01 + r00), k)
      rho1 = numeric(J); rho1[as.integer(rownames(r1_tmp))] = r1_tmp
      rho1 = rho1 / sum(rho1)
    }
    if(est_rho2){
      r2_tmp = rowsum(n * (r10 + r00), l)
      rho2 = numeric(J); rho2[as.integer(rownames(r2_tmp))] = r2_tmp
      rho2 = rho2 / sum(rho2)
    }
    # Expected latent-cell counts on the support
    C = numeric(length(s_code))
    if(use_T){
      # j = k contributions land at (Y1, X); the one-wrong (r01) mass is
      # spread over latent predecessors j of l: contribution to (j,i) is
      # pi_{j,i} * (T B)_{j,i} with B_{l,i} = sum of n*r01/m over cells
      C_tmp = rowsum(n * (r11 + r10), a_ki)
      C[as.integer(rownames(C_tmp))] = C_tmp
      Bmat = Matrix::sparseMatrix(i = l, j = x, x = n * r01 / pmax(mv, 1e-300),
                                  dims = c(J, K))
      tb = Matrix::summary(Tmat %*% Bmat)
      tbv = tb$x[match(s_code, tb$i + as.numeric(J)*(tb$j - 1))]
      tbv[is.na(tbv)] = 0
      C = C + p * tbv
    } else {
      C_tmp = rowsum(c(n * (r11 + r10), n * r01), c(a_ki, a_li))
      C[as.integer(rownames(C_tmp))] = C_tmp
    }
    A_i = as.numeric(rowsum(n * r00, gc_i))    # over positions present in gc_i
    A_full = numeric(length(ui)); A_full[sort(unique(gc_i))] = A_i
    p = (C + A_full[gS_i] * p / pmax(pplus[gS_i], 1e-300)) / N

  }

  if(!converged){
    warning("EM did not converge within `maxit` iterations; increase `maxit` or loosen `tol`.")
  }
  if(alpha1 > 0.5 || alpha2 > 0.5){
    warning("An estimated alpha exceeds 0.5: the diagonal-dominance identification condition (Assumption 3) fails at the estimate. Interpret with caution.")
  }

  return(list(
    alpha = c(alpha1 = alpha1, alpha2 = alpha2),
    rho1 = rho1, rho2 = rho2,
    Pi = data.frame(j = S_j, i = S_i, p = p),
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    n_iter = length(loglik_trace),
    converged = converged
  ))
}
