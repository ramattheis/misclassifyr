#' Shared linkage error across outcome views of the same links
#'
#' Estimates the trajectory model of [misclassifyr_traj_em()] JOINTLY on
#' two or more \emph{views} -- tabulations of the SAME linked men through
#' different outcome alphabets, for example counties in one view and
#' occupation groups in another. The record-linkage component -- the
#' shared-rival probability \eqn{\pi_s} and the fresh-failure rates
#' \eqn{b_1, \dots, b_k}, hence the marginal false-link rates
#' \eqn{\alpha_j} and the shared share \eqn{s} -- is FORCED equal across
#' views, because a link fails or holds once, whatever outcome is read
#' off the linked record. Everything else is view-specific: `Pi`, `T`,
#' the measurement-error layer (`kernel`, `mu`), the anchor-link rate
#' `alpha_f`, and the alphabet size `J`.
#'
#' Fitting the views separately and comparing their alphas is a
#' diagnostic; this function is the estimator that diagnostic motivates.
#' When one view's outcome carries measurement error of its own (an
#' occupation view, say, with coding error a county view does not have),
#' the separate fit parks that error in alpha, where it does not belong.
#' The joint fit pins the linkage component with every view's
#' information at once and pushes what remains into the view's own
#' measurement layer -- so give the noisy view a `kernel`, or the
#' constraint will simply strain against the misfit.
#'
#' @section Algorithm:
#' Block generalized EM on the composite (sum over views) log
#' likelihood. Each outer iteration runs, per view, one EM step of
#' [misclassifyr_traj_em()] with the shared parameters frozen (updating
#' the view-locals), then pools the views' expected-count aggregates --
#' the sufficient statistics whose ratios are the \eqn{\pi_s} and
#' \eqn{b_j} M-steps -- and updates the shared parameters exactly. Both
#' blocks ascend, so the composite trace is monotone. A man observed in
#' several views contributes to each view's term, so the objective is a
#' composite likelihood: bootstrap for standard errors, never the
#' curvature.
#'
#' @param views A named list. Each element is a list of arguments for
#'   [misclassifyr_traj_em()] -- `tab`, `J`, `rho`, and optionally
#'   `rho_f`, `kernel`, `dstep`, `dstep_f`, `mu_index`, `mu_fixed`,
#'   `T_fixed`, `Pi_fixed`, `chunk_size`. All views must observe the
#'   same number of linked measures `k`. Do not pass `alpha_fixed`,
#'   `s_fixed`, `init`, `maxit`, `squarem`, or `tol` here; the driver
#'   owns them.
#' @param alpha_0 Starting marginal false-link rate for every measure.
#' @param s_0 Starting shared-failure share.
#' @param tol Relative tolerance on the composite log likelihood.
#' @param maxit Maximum outer iterations.
#' @param verbose Print the composite log likelihood every 10 iterations.
#' @return A list: `alpha` (shared marginals), `s`, `pi_shared`, `b`,
#'   `views` (per-view lists with `T`, `Pi`, `mu`, `alpha_f`, `N`,
#'   `loglik`), `loglik` (composite), `loglik_trace`, `n_iter`,
#'   `converged`.
#' @seealso [misclassifyr_traj_em()] for one view;
#'   [misclassifyr_rl_em_shared()] for the two-link analogue.
#' @examples
#' \donttest{
#' # The same men, the same two links, read through a fine outcome
#' # (12 states) and a coarse one (4); a shared rival behind 40% of
#' # first-link failures. One alpha process serves both views.
#' set.seed(7)
#' N <- 8000; k <- 2; J1 <- 12; J2 <- 4
#' pis <- 0.08; b <- c(0.10, 0.15)
#' x1 <- sample.int(J1, N, TRUE); x2 <- (x1 - 1) %/% 3 + 1
#' shared <- runif(N) < pis
#' r1 <- sample.int(J1, N, TRUE)   # the shared rival's fine state
#' f <- sapply(b, function(bb) !shared & runif(N) < bb)
#' y1f <- ifelse(shared, r1, ifelse(f[, 1], sample.int(J1, N, TRUE), x1))
#' y2f <- ifelse(shared, r1, ifelse(f[, 2], sample.int(J1, N, TRUE), x1))
#' tab1 <- aggregate(list(n = rep(1, N)),
#'                   by = list(X = x1, Y1 = y1f, Y2 = y2f), FUN = sum)
#' cg <- function(v) (v - 1) %/% 3 + 1
#' tab2 <- aggregate(list(n = rep(1, N)),
#'                   by = list(X = x2, Y1 = cg(y1f), Y2 = cg(y2f)), FUN = sum)
#' rho1 <- lapply(1:2, function(j) matrix(1 / J1, J1, J1))
#' rho2 <- lapply(1:2, function(j) matrix(1 / J2, J2, J2))
#' fit <- misclassifyr_traj_em_shared(list(
#'   fine   = list(tab = tab1, J = J1, rho = rho1),
#'   coarse = list(tab = tab2, J = J2, rho = rho2)))
#' fit$alpha; fit$s
#' }
#' @export
misclassifyr_traj_em_shared = function(views,
                                       alpha_0 = 0.2, s_0 = 0.3,
                                       tol = 1e-8, maxit = 200,
                                       verbose = FALSE){

  if(!is.list(views) || is.data.frame(views) || length(views) < 1){
    stop("`views` should be a non-empty list of view specifications.")
  }
  if(is.null(names(views)) || any(!nzchar(names(views)))){
    names(views) = paste0("view", seq_along(views))
  }
  owned = c("alpha_fixed", "s_fixed", "init", "maxit", "squarem", "tol")
  for(nm in names(views)){
    bad = intersect(owned, names(views[[nm]]))
    if(length(bad)) stop("view `", nm, "` passes `", bad[1],
                         "`, which the driver owns.")
  }

  V = length(views)
  locals = vector("list", V)

  # k from the tabulations directly, so the shared parameters can be
  # sized and pinned from the very first inner call
  kv = vapply(views, function(vw)
    sum(grepl("^Y[0-9]+$", colnames(vw$tab))), integer(1))
  if(length(unique(kv)) != 1L)
    stop("views observe different numbers of linked measures (",
         paste(kv, collapse = ", "),
         "); shared linkage error requires the same links.")
  k = kv[1]

  pis = min(max(s_0 * alpha_0, 1e-8), 0.95)
  bshared = rep(min(max((alpha_0 - pis) / (1 - pis), 1e-6), 0.99), k)

  loglik_trace = numeric(0); ll_old = -Inf; converged = FALSE
  fits = vector("list", V)

  for(iter in seq_len(maxit)){

    ll = 0; sss = 0; nn = 0
    failnum = NULL

    for(v in seq_len(V)){
      a_now = pis + (1 - pis) * bshared
      s_now = pis / a_now[1]
      args = c(views[[v]],
               list(alpha_fixed = a_now,
                    s_fixed = s_now,
                    init = if(is.null(locals[[v]])) list() else locals[[v]],
                    alpha_0 = alpha_0,
                    maxit = 1L, squarem = FALSE, tol = 0))
      fit = suppressWarnings(do.call(misclassifyr_traj_em, args))
      fits[[v]] = fit

      locals[[v]] = list(T = fit$T, Pi = fit$Pi,
                         mu = if(all(fit$mu == 0)) NULL else fit$mu,
                         alpha_f = if(is.na(fit$alpha_f)) NULL else fit$alpha_f)
      locals[[v]] = locals[[v]][!vapply(locals[[v]], is.null, TRUE)]

      ec = fit$Ecounts
      ll = ll + ec$ll
      sss = sss + ec$R_SSS
      nn = nn + ec$N
      fj = vapply(seq_len(k), function(j)
        sum(ec$R_C[bitwAnd(seq_len(2L^k) - 1L,
                           bitwShiftL(1L, j - 1L)) == 0L]), numeric(1))
      failnum = if(is.null(failnum)) fj else failnum + fj
    }

    # Pooled M-step for the shared component: exact count ratios
    pis = min(max(sss / nn, 1e-10), 0.95)
    ind = nn - sss
    bshared = pmin(pmax(failnum / max(ind, 1e-12), 1e-10), 0.999)

    loglik_trace = c(loglik_trace, ll)
    if(verbose && iter %% 10 == 0) cat("outer", iter, "composite loglik", ll, "\n")
    if(is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + 1)){
      converged = TRUE; break
    }
    ll_old = ll
  }

  if(!converged)
    warning("shared EM did not converge within `maxit` outer iterations.")

  alpha = pis + (1 - pis) * bshared
  names(alpha) = paste0("alpha", seq_len(k))
  if(any(alpha > 0.5))
    warning("An estimated false-link rate exceeds 0.5: diagonal dominance fails at the estimate.")
  bvec = bshared; names(bvec) = paste0("b", seq_len(k))

  out_views = lapply(seq_len(V), function(v){
    f = fits[[v]]
    list(T = f$T, Pi = f$Pi, mu = f$mu, alpha_f = f$alpha_f,
         N = f$Ecounts$N, loglik = f$Ecounts$ll)
  })
  names(out_views) = names(views)

  return(list(
    alpha = alpha,
    s = if(alpha[1] > 0) unname(pis / alpha[1]) else 0,
    pi_shared = pis,
    b = bvec,
    views = out_views,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    n_iter = length(loglik_trace),
    converged = converged
  ))
}
