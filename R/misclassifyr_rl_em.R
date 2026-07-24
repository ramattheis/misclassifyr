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
#' @param alpha_0 Numeric starting value for both alpha parameters.
#' @param tol Relative log-likelihood convergence tolerance.
#' @param maxit Maximum EM iterations.
#' @param verbose Print the log likelihood every 25 iterations.
#' @return A list: `alpha` (length 2), `rho1`, `rho2`, `Pi` (data.frame
#'   `j`, `i`, `p` over the support), `loglik` (final), `loglik_trace`,
#'   `n_iter`, `converged`.
#' @export
misclassifyr_rl_em = function(tab, J, K,
                              rho1 = NULL, rho2 = NULL,
                              alpha_0 = 0.2,
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
  alpha1 = alpha_0; alpha2 = alpha_0

  #------------------------------------------------------------
  # EM iterations
  #------------------------------------------------------------

  loglik_trace = numeric(0)
  ll_old = -Inf
  converged = FALSE

  for(iter in seq_len(maxit)){

    # E-step ---------------------------------------------------
    pplus  = as.numeric(rowsum(p, gS_i))       # Pi_{+,i} over ui
    pk = p[a_ki]; pl = p[a_li]
    w11 = (1 - alpha1) * (1 - alpha2) * pk * same
    w10 = (1 - alpha1) * alpha2 * rho2[l] * pk
    w01 = alpha1 * rho1[k] * (1 - alpha2) * pl
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
    alpha1 = sum(n * (r01 + r00)) / N
    alpha2 = sum(n * (r10 + r00)) / N
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
    C_tmp = rowsum(c(n * (r11 + r10), n * r01), c(a_ki, a_li))
    C[as.integer(rownames(C_tmp))] = C_tmp
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
