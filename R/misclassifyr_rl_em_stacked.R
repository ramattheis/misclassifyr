#' Stacked EM for the record-linkage model: shared alpha across cells
#'
#' Estimates the record-linkage misclassification model on a LIST of
#' tabulations split by a conditioning variable W (e.g. birthplace), with
#' the false-link rates (alpha1, alpha2) SHARED across cells and the
#' failed-link draw distributions rho and the joint distribution Pi
#' CELL-SPECIFIC. This implements the conditional-independence version of
#' the model: phantom draws for a failed link come from the candidate pool
#' matching the linking keys (e.g. same birthplace), not the unconditional
#' marginal — an unconditional slab that is too diffuse makes false links
#' look correct and biases alpha downward.
#'
#' The E-step runs per cell with the common alpha; the M-step pools the
#' expected failed-link counts across cells for alpha and updates rho_c,
#' Pi_c within cells. An optional common T2 operator (see
#' \code{\link{misclassifyr_rl_em}}) applies to every cell.
#'
#' @param tabs A list of data.frames, each as in `misclassifyr_rl_em`
#'   (columns `X`, `Y1`, `Y2`, `n`, integer codes on a COMMON 1..J / 1..K
#'   coding across cells).
#' @param J,K Integers, the common category counts.
#' @param T2 Optional common true-transition operator (data.frame `j`,
#'   `l`, `t`), as in `misclassifyr_rl_em`.
#' @param alpha_fixed Optional length-2 vector fixing (alpha1, alpha2).
#' @param alpha_0,tol,maxit,verbose As in `misclassifyr_rl_em`.
#' @return A list: `alpha` (shared), `cells` (per-cell lists with `rho1`,
#'   `rho2`, `Pi`, `N`), `loglik`, `loglik_trace`, `n_iter`, `converged`.
#' @seealso [misclassifyr_rl_em()] for the single-tabulation version, and
#'   `vignette("designing-misclassification-models")` for when to stack.
#' @examples
#' # Four birthplace cells; a failed link draws a phantom from within the
#' # SAME cell, because the linking algorithm matched on birthplace.
#' set.seed(41)
#' C <- 4; J <- 20; N <- 40000; alpha <- 0.25
#' cell <- sample.int(C, N, replace = TRUE)
#' base <- (cell - 1) * 5
#' xi <- base + sample.int(5, N, replace = TRUE)
#' j  <- ifelse(runif(N) < 0.8, xi, base + sample.int(5, N, replace = TRUE))
#' y1 <- ifelse(runif(N) < alpha, base + sample.int(5, N, replace = TRUE), j)
#' y2 <- ifelse(runif(N) < alpha, base + sample.int(5, N, replace = TRUE), j)
#' df <- data.frame(cell = cell, X = xi, Y1 = y1, Y2 = y2)
#'
#' tabs <- lapply(split(df, df$cell), function(d)
#'   aggregate(list(n = rep(1, nrow(d))),
#'             by = list(X = d$X, Y1 = d$Y1, Y2 = d$Y2), FUN = sum))
#'
#' stacked <- misclassifyr_rl_em_stacked(tabs, J = J, K = J)
#' stacked$alpha   # close to (0.25, 0.25)
#'
#' # Pooling the cells and using the unconditional margin as the phantom
#' # distribution biases alpha downward
#' pooled_tab <- aggregate(list(n = rep(1, N)),
#'                         by = list(X = df$X, Y1 = df$Y1, Y2 = df$Y2), FUN = sum)
#' misclassifyr_rl_em(pooled_tab, J = J, K = J)$alpha
#' @export
misclassifyr_rl_em_stacked = function(tabs, J, K,
                                      T2 = NULL,
                                      alpha_fixed = NULL,
                                      alpha_0 = 0.2,
                                      tol = 1e-8, maxit = 500,
                                      verbose = FALSE){

  if(!is.list(tabs) || is.data.frame(tabs) || length(tabs) < 1){
    stop("`tabs` should be a non-empty list of tabulations.")
  }

  #------------------------------------------------------------
  # Per-cell precomputation (mirrors misclassifyr_rl_em)
  #------------------------------------------------------------

  use_T = !is.null(T2)
  if(use_T){
    if(!all(c("j","l","t") %in% colnames(T2))) stop("`T2` needs columns `j`, `l`, `t`.")
    Tj = as.integer(T2$j); Tl = as.integer(T2$l); Tt = as.numeric(T2$t)
    rs = rowsum(Tt, Tj)
    if(any(abs(rs - 1) > 1e-6)) stop("`T2` rows must sum to one within `j`.")
    Tmat = Matrix::sparseMatrix(i = Tj, j = Tl, x = Tt, dims = c(J, J))
    t_code = Tj + as.numeric(J)*(Tl - 1)
  }

  prep_cell = function(tab){
    tab = tab[tab$n > 0, , drop = FALSE]
    x = as.integer(tab$X); k = as.integer(tab$Y1); l = as.integer(tab$Y2)
    n = as.numeric(tab$n)
    if(any(is.na(x)) || any(x < 1) || any(x > K) ||
       any(is.na(k)) || any(k < 1) || any(k > J) ||
       any(is.na(l)) || any(l < 1) || any(l > J)){
      stop("cell codes must be integers in 1..K / 1..J.")
    }
    s_code = sort(unique(c(k + J*(x - 1), l + J*(x - 1))))
    S_j = ((s_code - 1) %% J) + 1
    S_i = ((s_code - 1) %/% J) + 1
    a_ki = match(k + J*(x - 1), s_code)
    a_li = match(l + J*(x - 1), s_code)
    ui   = sort(unique(S_i))
    e = list(x = x, k = k, l = l, n = n, N = sum(n),
             s_code = s_code, S_j = S_j, S_i = S_i,
             a_ki = a_ki, a_li = a_li, ui = ui,
             gS_i = match(S_i, ui), gc_i = match(x, ui))
    e$t_kl = if(use_T){
      v = Tt[match(k + as.numeric(J)*(l - 1), t_code)]; v[is.na(v)] = 0; v
    } else as.numeric(k == l)
    e$cell_li_code = l + as.numeric(J)*(x - 1)
    # Initial values
    p = numeric(length(s_code))
    p_tmp = rowsum(c(n, n) / (2 * e$N), c(a_ki, a_li))
    p[as.integer(rownames(p_tmp))] = p_tmp
    e$p = p
    r1 = numeric(J); r1t = rowsum(n, k); r1[as.integer(rownames(r1t))] = r1t / e$N
    r2 = numeric(J); r2t = rowsum(n, l); r2[as.integer(rownames(r2t))] = r2t / e$N
    e$rho1 = r1; e$rho2 = r2
    e
  }
  cells = lapply(tabs, prep_cell)
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
  # EM: E-step per cell with common alpha; alpha pooled in M-step
  #------------------------------------------------------------

  loglik_trace = numeric(0); ll_old = -Inf; converged = FALSE

  for(iter in seq_len(maxit)){

    ll = 0; fail1 = 0; fail2 = 0

    for(ci in seq_along(cells)){
      e = cells[[ci]]
      pplus = as.numeric(rowsum(e$p, e$gS_i))
      pk = e$p[e$a_ki]
      if(use_T){
        P  = Matrix::sparseMatrix(i = e$S_j, j = e$S_i, x = e$p, dims = c(J, K))
        ms = Matrix::summary(Matrix::crossprod(Tmat, P))
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

      ll = ll + sum(e$n * log(pc))
      fail1 = fail1 + sum(e$n * (r01 + r00))
      fail2 = fail2 + sum(e$n * (r10 + r00))

      # Per-cell rho and Pi updates (stored for after alpha update; EM is
      # valid updating all parameters from the same responsibilities)
      r1t = rowsum(e$n * (r01 + r00), e$k)
      rho1 = numeric(J); rho1[as.integer(rownames(r1t))] = r1t
      e$rho1 = if(sum(rho1) > 0) rho1 / sum(rho1) else e$rho1
      r2t = rowsum(e$n * (r10 + r00), e$l)
      rho2 = numeric(J); rho2[as.integer(rownames(r2t))] = r2t
      e$rho2 = if(sum(rho2) > 0) rho2 / sum(rho2) else e$rho2

      C = numeric(length(e$s_code))
      if(use_T){
        C_tmp = rowsum(e$n * (r11 + r10), e$a_ki)
        C[as.integer(rownames(C_tmp))] = C_tmp
        Bmat = Matrix::sparseMatrix(i = e$l, j = e$x,
                                    x = e$n * r01 / pmax(mv, 1e-300),
                                    dims = c(J, K))
        tb = Matrix::summary(Tmat %*% Bmat)
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

  return(list(
    alpha = c(alpha1 = alpha1, alpha2 = alpha2),
    cells = lapply(cells, function(e)
      list(rho1 = e$rho1, rho2 = e$rho2,
           Pi = data.frame(j = e$S_j, i = e$S_i, p = e$p), N = e$N)),
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    n_iter = length(loglik_trace),
    converged = converged
  ))
}
