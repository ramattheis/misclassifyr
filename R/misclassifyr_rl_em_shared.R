#' Shared-alpha EM across heterogeneous designs
#'
#' Estimates the record-linkage misclassification model JOINTLY on two or
#' more designs that observe the SAME links through different outcomes --
#' for example a geography design (counties, with a migration operator)
#' and an occupation design (meso occupations, with its own transition
#' structure) built on the same linked men. The false-link rates
#' (alpha1, alpha2) are a property of the links, not of the outcome being
#' read off the linked record, so they are FORCED equal across designs;
#' everything else -- the latent distribution `Pi`, the failed-link draw
#' distributions `rho`, and the optional true-transition operator `T2` --
#' is design-specific, including the alphabet sizes `J` and `K`.
#'
#' This differs from [misclassifyr_rl_em_stacked()], which shares alpha
#' across conditioning cells of ONE design (common `J`, `K`, `T2`).
#' Here each design carries its own alphabet and operator. Fitting the
#' designs separately and comparing their alphas is a diagnostic; fitting
#' them jointly is the estimator the diagnostic motivates. When the
#' designs disagree under separate fits, the joint fit allocates the
#' disagreement: alpha is pinned by the design that identifies it best,
#' and what the other design was absorbing into alpha is pushed back into
#' its outcome-specific structure, where it belongs.
#'
#' The objective is a composite (pairwise-margin) likelihood: each
#' design's tabulation contributes its own term and the terms share
#' alpha. A man observed in several designs contributes to each, so the
#' trace is a composite log likelihood, monotone under this EM; model-
#' based standard errors from its curvature are NOT valid -- bootstrap
#' over men (or tabulation cells) instead.
#'
#' @param designs A named list. Each element is a list with components:
#'   `tab` (a data.frame with columns `X`, `Y1`, `Y2`, `n`, coded
#'   1..K / 1..J for THIS design), `J`, `K`, and optionally `T2` (a
#'   data.frame `j`, `l`, `t`; rows summing to one within `j`), as in
#'   [misclassifyr_rl_em()].
#' @param alpha_fixed Optional length-2 vector fixing (alpha1, alpha2).
#' @param alpha_0,tol,maxit,verbose As in [misclassifyr_rl_em()].
#' @return A list: `alpha` (shared), `designs` (per-design lists with
#'   `rho1`, `rho2`, `Pi`, `N`, `loglik`), `loglik` (composite),
#'   `loglik_trace`, `n_iter`, `converged`.
#' @seealso [misclassifyr_rl_em()], [misclassifyr_rl_em_stacked()].
#' @examples
#' # One population of men, two outcomes read through the SAME two links:
#' # a fine outcome (J = 12) and a coarse one (J = 4). Links fail at the
#' # same rate regardless of which outcome you read.
#' set.seed(93)
#' N <- 30000; alpha <- c(0.2, 0.3)
#' xf <- sample.int(12, N, replace = TRUE)        # fine latent outcome
#' xc <- (xf - 1) %/% 3 + 1                       # coarse = grouping of fine
#' fail1 <- runif(N) < alpha[1]; fail2 <- runif(N) < alpha[2]
#' rf1 <- sample.int(12, N, TRUE); rf2 <- sample.int(12, N, TRUE)
#' f1 <- ifelse(fail1, rf1, xf); f2 <- ifelse(fail2, rf2, xf)
#' c1 <- ifelse(fail1, (rf1 - 1) %/% 3 + 1, xc)
#' c2 <- ifelse(fail2, (rf2 - 1) %/% 3 + 1, xc)
#' tab_f <- aggregate(list(n = rep(1, N)),
#'                    by = list(X = xf, Y1 = f1, Y2 = f2), FUN = sum)
#' tab_c <- aggregate(list(n = rep(1, N)),
#'                    by = list(X = xc, Y1 = c1, Y2 = c2), FUN = sum)
#' fit <- misclassifyr_rl_em_shared(list(
#'   fine   = list(tab = tab_f, J = 12, K = 12),
#'   coarse = list(tab = tab_c, J = 4,  K = 4)))
#' fit$alpha   # close to (0.2, 0.3)
#' @export
misclassifyr_rl_em_shared = function(designs,
                                     alpha_fixed = NULL,
                                     alpha_0 = 0.2,
                                     tol = 1e-8, maxit = 500,
                                     verbose = FALSE){

  if(!is.list(designs) || is.data.frame(designs) || length(designs) < 1){
    stop("`designs` should be a non-empty list of design specifications.")
  }
  if(is.null(names(designs)) || any(!nzchar(names(designs)))){
    names(designs) = paste0("design", seq_along(designs))
  }

  #------------------------------------------------------------
  # Per-design precomputation: misclassifyr_rl_em_stacked's cell
  # machinery with J, K, and T2 local to the design
  #------------------------------------------------------------

  prep_design = function(d, nm){
    for(need in c("tab", "J", "K")){
      if(is.null(d[[need]])) stop("design `", nm, "` needs `", need, "`.")
    }
    J = as.integer(d$J); K = as.integer(d$K)
    tab = d$tab[d$tab$n > 0, , drop = FALSE]
    x = as.integer(tab$X); k = as.integer(tab$Y1); l = as.integer(tab$Y2)
    n = as.numeric(tab$n)
    if(any(is.na(x)) || any(x < 1) || any(x > K) ||
       any(is.na(k)) || any(k < 1) || any(k > J) ||
       any(is.na(l)) || any(l < 1) || any(l > J)){
      stop("design `", nm, "`: codes must be integers in 1..K / 1..J.")
    }
    e = list(J = J, K = K, use_T = !is.null(d$T2))
    if(e$use_T){
      T2 = d$T2
      if(!all(c("j","l","t") %in% colnames(T2)))
        stop("design `", nm, "`: `T2` needs columns `j`, `l`, `t`.")
      Tj = as.integer(T2$j); Tl = as.integer(T2$l); Tt = as.numeric(T2$t)
      rs = rowsum(Tt, Tj)
      if(any(abs(rs - 1) > 1e-6))
        stop("design `", nm, "`: `T2` rows must sum to one within `j`.")
      e$Tmat = Matrix::sparseMatrix(i = Tj, j = Tl, x = Tt, dims = c(J, J))
      e$Tt = Tt
      e$t_code = Tj + as.numeric(J)*(Tl - 1)
    }
    s_code = sort(unique(c(k + J*(x - 1), l + J*(x - 1))))
    e$s_code = s_code
    e$S_j = ((s_code - 1) %% J) + 1
    e$S_i = ((s_code - 1) %/% J) + 1
    e$a_ki = match(k + J*(x - 1), s_code)
    e$a_li = match(l + J*(x - 1), s_code)
    e$ui   = sort(unique(e$S_i))
    e$x = x; e$k = k; e$l = l; e$n = n; e$N = sum(n)
    e$gS_i = match(e$S_i, e$ui); e$gc_i = match(x, e$ui)
    e$t_kl = if(e$use_T){
      v = e$Tt[match(k + as.numeric(J)*(l - 1), e$t_code)]; v[is.na(v)] = 0; v
    } else as.numeric(k == l)
    e$cell_li_code = l + as.numeric(J)*(x - 1)
    p = numeric(length(s_code))
    p_tmp = rowsum(c(n, n) / (2 * e$N), c(e$a_ki, e$a_li))
    p[as.integer(rownames(p_tmp))] = p_tmp
    e$p = p
    r1 = numeric(J); r1t = rowsum(n, k); r1[as.integer(rownames(r1t))] = r1t / e$N
    r2 = numeric(J); r2t = rowsum(n, l); r2[as.integer(rownames(r2t))] = r2t / e$N
    e$rho1 = r1; e$rho2 = r2
    e
  }
  cells = mapply(prep_design, designs, names(designs), SIMPLIFY = FALSE)
  N_all = sum(sapply(cells, `[[`, "N"))

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
  # EM: E-step per design with the common alpha; alpha pooled
  #------------------------------------------------------------

  loglik_trace = numeric(0); ll_old = -Inf; converged = FALSE
  ll_by_design = numeric(length(cells))

  for(iter in seq_len(maxit)){

    ll = 0; fail1 = 0; fail2 = 0

    for(ci in seq_along(cells)){
      e = cells[[ci]]; J = e$J; K = e$K
      pplus = as.numeric(rowsum(e$p, e$gS_i))
      pk = e$p[e$a_ki]
      if(e$use_T){
        P  = Matrix::sparseMatrix(i = e$S_j, j = e$S_i, x = e$p, dims = c(J, K))
        ms = Matrix::summary(Matrix::crossprod(e$Tmat, P))
        mv = ms$x[match(e$cell_li_code, ms$i + as.numeric(J)*(ms$j - 1))]
        mv[is.na(mv)] = 0
      } else {
        mv = e$p[e$a_li]
      }
      w11 = (1 - alpha1) * (1 - alpha2) * pk * e$t_kl
      w10 = (1 - alpha1) * alpha2 * e$rho2[e$l] * pk
      w01 = alpha1 * e$rho1[e$k] * (1 - alpha2) * mv
      w00 = alpha1 * alpha2 * e$rho1[e$k] * e$rho2[e$l] * pplus[e$gc_i]
      pc  = w11 + w10 + w01 + w00
      r11 = w11/pc; r10 = w10/pc; r01 = w01/pc; r00 = w00/pc

      ll_by_design[ci] = sum(e$n * log(pc))
      ll = ll + ll_by_design[ci]
      fail1 = fail1 + sum(e$n * (r01 + r00))
      fail2 = fail2 + sum(e$n * (r10 + r00))

      r1t = rowsum(e$n * (r01 + r00), e$k)
      rho1 = numeric(J); rho1[as.integer(rownames(r1t))] = r1t
      e$rho1 = if(sum(rho1) > 0) rho1 / sum(rho1) else e$rho1
      r2t = rowsum(e$n * (r10 + r00), e$l)
      rho2 = numeric(J); rho2[as.integer(rownames(r2t))] = r2t
      e$rho2 = if(sum(rho2) > 0) rho2 / sum(rho2) else e$rho2

      C = numeric(length(e$s_code))
      if(e$use_T){
        C_tmp = rowsum(e$n * (r11 + r10), e$a_ki)
        C[as.integer(rownames(C_tmp))] = C_tmp
        Bmat = Matrix::sparseMatrix(i = e$l, j = e$x,
                                    x = e$n * r01 / pmax(mv, 1e-300),
                                    dims = c(J, K))
        tb = Matrix::summary(e$Tmat %*% Bmat)
        tbv = tb$x[match(e$s_code, tb$i + as.numeric(J)*(tb$j - 1))]
        tbv[is.na(tbv)] = 0
        C = C + e$p * tbv
      } else {
        C_tmp = rowsum(c(e$n * (r11 + r10), e$n * r01), c(e$a_ki, e$a_li))
        C[as.integer(rownames(C_tmp))] = C_tmp
      }
      A_i = as.numeric(rowsum(e$n * r00, e$gc_i))
      A_full = numeric(length(e$ui)); A_full[sort(unique(e$gc_i))] = A_i
      e$p = (C + A_full[e$gS_i] * e$p / pmax(pplus[e$gS_i], 1e-300)) / e$N
      cells[[ci]] = e
    }

    loglik_trace = c(loglik_trace, ll)
    if(verbose && iter %% 25 == 0) cat("iter", iter, "loglik", ll, "\n")
    if(is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + 1)){
      converged = TRUE; break
    }
    ll_old = ll

    if(est_alpha){
      alpha1 = fail1 / N_all
      alpha2 = fail2 / N_all
    }
  }

  if(!converged) warning("EM did not converge within `maxit` iterations.")
  if(alpha1 > 0.5 || alpha2 > 0.5){
    warning("An estimated alpha exceeds 0.5: diagonal dominance fails at the estimate.")
  }

  out_designs = lapply(seq_along(cells), function(ci){
    e = cells[[ci]]
    list(rho1 = e$rho1, rho2 = e$rho2,
         Pi = data.frame(j = e$S_j, i = e$S_i, p = e$p),
         N = e$N, loglik = ll_by_design[ci])
  })
  names(out_designs) = names(designs)

  return(list(
    alpha = c(alpha1 = alpha1, alpha2 = alpha2),
    designs = out_designs,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    n_iter = length(loglik_trace),
    converged = converged
  ))
}
