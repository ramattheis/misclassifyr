#' Internal: integer matrix power by repeated squaring
#'
#' @noRd
mat_power = function(M, d){
  R = diag(nrow(M))
  while(d > 0L){
    if(d %% 2L == 1L) R = R %*% M
    M = M %*% M
    d = d %/% 2L
  }
  R
}

#' Internal: T applied repeatedly to a backward matrix, keeping intermediates
#'
#' Returns a list of length `d + 1` whose element `t + 1` is
#' \eqn{T^{d-t} B}, so element `d + 1` is `B` itself and element `1` is
#' \eqn{T^d B}. `NULL` stands for the all-ones matrix, which `T` maps to
#' itself (rows sum to one), so the `NULL` marker propagates for free.
#'
#' @noRd
traj_tchain = function(Tm, d, B){
  out = vector("list", d + 1L)
  out[d + 1L] = list(B)                 # single bracket: NULL must not delete
  if(d > 0L){
    for(t in seq.int(d, 1L)){
      out[t] = list(if(is.null(out[[t + 1L]])) NULL else Tm %*% out[[t + 1L]])
    }
  }
  out
}

#' Internal: the own-chain backward trellis over correct/failed patterns
#'
#' `bw[[j]][[q + 1]]` is the `J x C` matrix of own-chain backward values at
#' measure `j`, where `q` encodes which of measures `j..k` are correctly
#' linked (bit 0 = measure `j`). `NULL` marks the all-ones matrix, i.e. a
#' suffix carrying no own emission. `bo[[j]][[q + 1]]` is the matching
#' per-cell logical "this suffix carries at least one emission", which is
#' what lets the T M-step skip transitions no observation informs.
#'
#' @noRd
traj_backward = function(Tm, dstep, mown, obsm, k, keep_mid){
  bw = vector("list", k); bo = vector("list", k)
  mid = if(keep_mid) vector("list", k) else NULL
  bw[[k]] = list(NULL, mown[[k]])
  bo[[k]] = list(rep(FALSE, length(obsm[[k]])), obsm[[k]])
  if(k > 1L){
    for(j in seq.int(k - 1L, 1L)){
      nq = 2L^(k - j)
      out = vector("list", 2L * nq); ob = vector("list", 2L * nq)
      md = if(keep_mid) vector("list", nq) else NULL
      for(q in seq_len(nq) - 1L){
        chain = traj_tchain(Tm, dstep[j], bw[[j + 1L]][[q + 1L]])
        st = chain[[1L]]
        if(keep_mid) md[[q + 1L]] = chain
        out[2L * q + 1L] = list(st)                   # measure j failed
        out[2L * q + 2L] = list(if(is.null(st)) mown[[j]] else mown[[j]] * st)
        ob[[2L * q + 1L]] = bo[[j + 1L]][[q + 1L]]
        ob[[2L * q + 2L]] = bo[[j + 1L]][[q + 1L]] | obsm[[j]]
      }
      bw[[j]] = out; bo[[j]] = ob
      if(keep_mid) mid[[j]] = md
    }
  }
  list(bw = bw, bo = bo, mid = mid)
}

#' Internal: re-emit the trellis after swapping measure `j`'s own emission
#'
#' Levels above `j` are untouched (the observation pattern is unchanged), and
#' at level `j` the "measure j failed" entries carry no emission at `j`, so
#' they are reused as well. Only levels below `j` need fresh `T` products.
#'
#' @noRd
traj_backward_swap = function(Tm, dstep, mown2, k, j, bw){
  for(lev in seq.int(j, 1L)){
    nq = 2L^(k - lev)
    out = vector("list", 2L * nq)
    for(q in seq_len(nq) - 1L){
      st = if(lev == j) bw[[lev]][[2L * q + 1L]] else
        traj_tchain(Tm, dstep[lev], bw[[lev + 1L]][[q + 1L]])[[1L]]
      out[2L * q + 1L] = list(st)
      out[2L * q + 2L] = list(if(is.null(st)) mown2[[lev]] else mown2[[lev]] * st)
    }
    bw[[lev]] = out
  }
  bw
}

#' Internal: the son-side factor P(Y_1..Y_k | anchor state i)
#'
#' Returns the `J x C` factor and, when `keep = TRUE`, the per-branch own-chain
#' factors `A` and fresh-rival factors `Phi` the M-step reuses.
#'
#' @noRd
traj_son = function(Pi, pi0, mown, efresh, bw, wC, pis, rho1, k, J, Cc,
                    keep = FALSE){
  Son = matrix(0, J, Cc)
  nm = 2L^k
  Alist = if(keep) vector("list", nm) else NULL
  Plist = if(keep) vector("list", nm) else NULL
  for(mask in seq_len(nm) - 1L){
    b1 = bw[[1L]][[mask + 1L]]
    A = if(is.null(b1)) matrix(pi0, J, Cc) else Pi %*% b1
    Phi = NULL
    for(j in seq_len(k)){
      if(bitwAnd(mask, bitwShiftL(1L, j - 1L)) == 0L){
        Phi = if(is.null(Phi)) efresh[[j]] else Phi * efresh[[j]]
      }
    }
    Son = Son + wC[mask + 1L] * (if(is.null(Phi)) A else A * Phi)
    if(keep){ Alist[[mask + 1L]] = A; Plist[mask + 1L] = list(Phi) }
  }
  bfull = bw[[1L]][[nm]]
  Rv = if(is.null(bfull)) matrix(1, J, Cc) else rho1 %*% bfull
  Son = Son + pis * (pi0 * Rv)
  if(!keep) return(list(Son = Son))
  list(Son = Son, A = Alist, Phi = Plist, Rv = Rv)
}

#' EM for multi-link bundles: structured slab, anchor branch, and a
#' measurement-error layer
#'
#' Estimates the trajectory model behind the three-link design: one unit
#' contributes an anchor measure observed directly (in the application, the
#' father's occupation in the base census), an optional second anchor
#' measure reached by a record link of its own, and `k` linked measures of a
#' second person whose latent state follows a Markov chain (the son, linked
#' into `k` later censuses). Three layers are estimated jointly — the
#' intergenerational joint distribution `Pi`, the transition matrix `T`
#' governing latent drift between censuses, and the linkage and
#' measurement-error rates — and every M-step is closed form.
#'
#' @section The structured slab:
#' A failed link attaches the unit to somebody else's record. The rank-one
#' model behind [misclassifyr_rl_em()] treats each failure as an independent
#' draw, which makes two failed links disagree as often as chance allows.
#' That is false: a linking algorithm confused by one rival tends to stay
#' confused by the *same* rival. This estimator gives each unit **one latent
#' rival**, whose own state follows the same chain `T`. With probability
#' \eqn{\pi_s} all `k` links fail to that shared rival; otherwise each link
#' fails independently with probability \eqn{b_j} to a fresh rival drawn from
#' the rival law `rho`. The marginal false-link rate is
#' \eqn{\alpha_j = \pi_s + (1-\pi_s) b_j} and the share of failures that are
#' shared is \eqn{s = \pi_s / \alpha_1}. At `k = 2` this is a five-component
#' mixture (both correct; either one correct; both failed to fresh rivals;
#' both failed to the shared rival); the general-`k` analogue has
#' \eqn{2^k + 1} components and is what the code enumerates.
#'
#' At `k = 2` the likelihood in \eqn{(\alpha, s)} is a ridge — the two trade
#' off almost exactly — and `s` is identified only at `k >= 3`, where the
#' run-length pattern of agreements across three links separates a shared
#' rival from a coincidence. That contrast is the design motivation for the
#' three-link bundle and is exercised in the package's tests.
#'
#' @section The measurement-error layer:
#' Every *observed* category, the anchor's included, passes through
#' \eqn{M = (1-\mu) I + \mu K} for a supplied confusion kernel `K`: with
#' probability \eqn{\mu} the recorded value is not the true one but a draw
#' from `K`'s row. `mu` may vary by census year through `mu_index`; `K` is
#' data, not a parameter. Set `kernel = NULL` to switch the layer off
#' (\eqn{\mu = 0} for every slot), which is the pure linkage model.
#'
#' The kernel matters. Simulation (gate G0) shows that fitting a
#' within-block kernel to rank-distance confusion halves the implied error
#' rate and dumps the residue on `T`, whose mean diagonal fell to 0.655
#' against a truth of 0.750 — a quarter of true persistence misread as
#' mobility, while the estimated `Pi` was barely touched. Never report `T`
#' without a kernel-robustness pass.
#'
#' @section Missing measures:
#' `NA` (or `0`) in any measure column means unobserved, and is handled as a
#' missing emission: the unit contributes exactly the branches its observed
#' pattern supports, with no imputation and no dropping. This is the
#' deliberate design. Modelling "not observed" as an extra latent state
#' instead degenerates — in simulation `s` pinned at its upper boundary and
#' the false-link rates scattered over (0.07, 0.60, 0.87) against a truth of
#' 0.30 — because a high-mass "no occupation" category is agreed on by two
#' links for reasons that have nothing to do with linkage, and the model
#' reads that agreement as correct linking. Do not model an unobserved
#' measure as a state.
#'
#' @section Closed forms:
#' Every M-step is a count ratio, and no Newton step is used anywhere:
#' \eqn{\pi_s} and \eqn{\alpha_f} are branch-responsibility shares,
#' \eqn{b_j} the share of independent-branch mass failing at measure `j`,
#' `Pi` the expected latent-cell counts, `T` the row-normalised pooled
#' expected transition counts (son-chain steps, shared-rival-chain steps and
#' the anchor's own link step), and `mu` the expected share of emissions
#' drawn from `K`, obtained exactly from the emission-swap identity
#' \eqn{E[\#K \text{ at slot } t] = \sum_c (n_c/P_c) P^{(t,\mu K)}_c}. This
#' holds **only because the rival laws `rho` are data**. Rebuilding them
#' from the current `T` each iteration breaks the factorization and, in
#' simulation, breaks recovery; `rho` is never updated.
#'
#' @param tab A data.frame of cell counts with columns `Y1`, ..., `Yk` (the
#'   linked measures of the traced person, integer codes in `1..J`), `X`
#'   (the anchor's directly observed measure), optionally `Xf` (the anchor's
#'   own linked measure), and `n` (non-negative counts). `NA` or `0` marks a
#'   missing measure. `k` is inferred from the `Y` columns. Only observed
#'   cells need be present; unit-level rows with `n = 1` are fine, but
#'   identical patterns must be aggregated first.
#' @param J An integer, the number of latent states. The anchor and the
#'   traced person share the state space, because they share `T`.
#' @param rho The rival (rival) laws, treated as DATA and never updated: a
#'   `J x J` row-stochastic matrix whose row `i` is the distribution of a
#'   fresh rival's latent state given anchor state `i`, or a list of `k` such
#'   matrices (one per linked measure), or a length-`J` probability vector
#'   for an anchor-independent slab. Derive these from observed margins,
#'   localised the way the linking keys localise the candidate pool. When
#'   `kernel` is supplied these are laws over the LATENT state: the model
#'   applies the measurement layer to them itself, so an observed margin
#'   must be deconvolved first (or `kernel` left `NULL`).
#' @param rho_f The rival law for the anchor's own link, as a `J x J`
#'   matrix or a length-`J` vector. Required when `tab` has an `Xf` column.
#' @param kernel The `J x J` row-stochastic confusion kernel `K`. `NULL`
#'   (default) switches the measurement-error layer off.
#' @param dstep Integer vector of length `k - 1`: the number of unit `T`
#'   steps between consecutive linked measures. Defaults to all ones (one
#'   step per measure). `0` means two measures at the same date.
#' @param dstep_f Integer, the number of unit `T` steps between the anchor's
#'   base measure and its linked measure. Defaults to `1`.
#' @param mu_index Integer vector of length `k + 1` (or `k + 2` when `Xf` is
#'   present) mapping each measurement slot — in the order `Y1`, ..., `Yk`,
#'   `X`, `Xf` — to a `mu` parameter. Defaults to one common `mu`. Supply
#'   e.g. `c(1, 2, 3, 1, 1)` for a per-census-year `mu`.
#' @param alpha_0 Numeric starting value for the marginal false-link rates.
#' @param init Optional list of starting values with components `s`, `mu`,
#'   `alpha_f`, `T` and `Pi`. Anything absent gets a default: `Pi` from the
#'   observed `(X, Yj)` joints, `T` a diagonally dominant matrix, and scalars
#'   at mild values. Multi-start is not optional in practice — the
#'   \eqn{(\alpha, s)} block plateaus badly.
#' @param alpha_fixed Optional length-`k` vector of marginal false-link
#'   rates to hold fixed, `NA` for entries left free. Must be supplied
#'   together with `s_fixed`, and `alpha_fixed[1]` must be non-`NA`: the slab
#'   is parameterised by \eqn{(\pi_s, b)}, and only the pair
#'   \eqn{(\alpha_1, s)} pins \eqn{\pi_s = s\,\alpha_1} in closed form.
#'   Fixing both is how the \eqn{(\alpha, s)} profile is traced.
#' @param s_fixed Optional shared-failure share to hold fixed; see
#'   `alpha_fixed`.
#' @param pi_shared_fixed Optional value of \eqn{\pi_s} (the probability
#'   that every link fails to the same rival) to hold fixed while `b` is
#'   re-maximised. This is the profiling handle for `s`: sweep it over a grid
#'   and read the implied `s` off each fit. Mutually exclusive with the
#'   (`alpha_fixed`, `s_fixed`) pair.
#' @param alpha_f_fixed Optional false-link rate for the anchor's own link.
#' @param mu_fixed Optional vector of `mu` values to hold fixed, one per
#'   `mu_index` group (recycled if length one), `NA` for free entries.
#' @param T_fixed Optional `J x J` row-stochastic transition matrix held
#'   fixed as a plug-in instead of estimated.
#' @param Pi_fixed Optional `J x J` joint distribution (summing to one) held
#'   fixed.
#' @param chunk_size Integer, the number of cells processed at a time.
#'   Working memory is `O(2^k J chunk_size)`; lower this on wide state
#'   spaces. `NULL` processes every cell at once. Time scales the same way:
#'   the branch enumeration has `2^k + 1` components, so `k` much beyond 5
#'   is impractical.
#' @param tol Relative log-likelihood convergence tolerance. The default is
#'   tighter than the package's other EM estimators use on purpose: the
#'   \eqn{(\alpha, s)} block crawls across a near-plateau, and a loose
#'   tolerance stops on it. In simulation, `tol = 1e-8` halted at
#'   \eqn{\alpha = 0.18} where the MLE sat at 0.46, with only 0.9 log units
#'   between them.
#' @param maxit Maximum EM iterations. Hundreds are normal; convergence in
#'   a few dozen usually means the tolerance is too loose.
#' @param squarem Logical: accelerate with SQUAREM. Strongly recommended —
#'   plain EM plateaus on the \eqn{(\alpha, s)} block.
#' @param verbose Print the log likelihood every 25 iterations.
#' @return A list: `alpha` (length `k` marginal false-link rates), `s`,
#'   `pi_shared`, `b`, `alpha_f`, `mu`, `T`, `Pi` (a `J x J` matrix, rows the
#'   anchor's latent state, columns the traced person's latent state at the
#'   first linked measure), `k`, `loglik`, `loglik_trace`, `n_iter`,
#'   `n_feval` and `converged`.
#' @seealso [misclassifyr_rl_em()] for the rank-one two-link model this
#'   generalises, and [misclassifyr_known_slab()] for the single-link case
#'   with a known rival law.
#' @examples
#' # Three links of one traced person, plus an anchor observed twice.
#' # Latent states drift one step per census; a failed link attaches the
#' # unit to a rival, and with probability s to the SAME rival every time.
#' set.seed(4)
#' J <- 4; N <- 6000
#' Tm <- 0.7 * diag(J) + 0.3 / J; Tm <- Tm / rowSums(Tm)
#' rho <- matrix(1 / J, J, J)
#' draw <- function(P, x) {
#'   cp <- t(apply(P, 1, cumsum))
#'   as.integer(rowSums(runif(length(x)) > cp[x, , drop = FALSE]) + 1L)
#' }
#' alpha <- 0.3; s <- 0.6; pis <- s * alpha; b <- (alpha - pis) / (1 - pis)
#'
#' i <- sample.int(J, N, TRUE)                 # anchor latent state
#' s1 <- draw(Tm %*% Tm, i)                    # traced person, measure 1
#' s2 <- draw(Tm, s1); s3 <- draw(Tm, s2)
#' r1 <- sample.int(J, N, TRUE)                # the one latent rival
#' r2 <- draw(Tm, r1); r3 <- draw(Tm, r2)
#' shared <- runif(N) < pis
#' fail <- function() shared | (runif(N) < b)
#' pick <- function(f, own, riv) ifelse(!f, own,
#'   ifelse(shared, riv, sample.int(J, N, TRUE)))
#' f1 <- fail(); f2 <- fail(); f3 <- fail()
#' dat <- data.frame(X = i, Xf = ifelse(runif(N) < 0.25,
#'                                      sample.int(J, N, TRUE), draw(Tm, i)),
#'                   Y1 = pick(f1, s1, r1), Y2 = pick(f2, s2, r2),
#'                   Y3 = pick(f3, s3, r3))
#' tab <- aggregate(list(n = rep(1, N)), by = as.list(dat), FUN = sum)
#'
#' fit <- misclassifyr_traj_em(tab, J = J, rho = rho, rho_f = rho,
#'                             maxit = 300)
#' fit$alpha        # near 0.3
#' fit$s            # near 0.6
#' fit$alpha_f      # near 0.25
#' round(fit$Pi, 3) # the estimand: anchor state x traced state
#' @export
misclassifyr_traj_em = function(tab, J, rho, rho_f = NULL, kernel = NULL,
                                dstep = NULL, dstep_f = 1L, mu_index = NULL,
                                alpha_0 = 0.2, init = list(),
                                alpha_fixed = NULL, s_fixed = NULL,
                                pi_shared_fixed = NULL,
                                alpha_f_fixed = NULL, mu_fixed = NULL,
                                T_fixed = NULL, Pi_fixed = NULL,
                                chunk_size = 20000L,
                                tol = 1e-10, maxit = 1000, squarem = TRUE,
                                verbose = FALSE){

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  if(!is.data.frame(tab)) stop("`tab` should be a data.frame.")
  nms = colnames(tab)
  if(!("n" %in% nms)) stop("`tab` should have a column `n` of counts.")
  if(!("X" %in% nms))
    stop("`tab` should have a column `X`: the anchor's directly observed measure.")
  ycol = grep("^Y[0-9]+$", nms, value = TRUE)
  if(length(ycol) < 2L)
    stop("`tab` should have at least two linked-measure columns `Y1`, `Y2`, ...")
  ycol = ycol[order(as.integer(sub("^Y", "", ycol)))]
  if(!identical(as.integer(sub("^Y", "", ycol)), seq_along(ycol)))
    stop("the linked-measure columns must be named `Y1` ... `Yk` with no gaps.")
  k = length(ycol)
  has_f = "Xf" %in% nms

  J = as.integer(J)
  if(is.na(J) || J < 2L) stop("`J` must be an integer of at least 2.")

  tab = tab[!is.na(tab$n) & tab$n > 0, , drop = FALSE]
  if(nrow(tab) == 0L) stop("`tab` has no rows with a positive `n`.")
  n = as.numeric(tab$n); N = sum(n); nc = length(n)

  code_col = function(v, nm){
    if(is.factor(v)) v = as.character(v)
    v = suppressWarnings(as.integer(v))
    v[!is.na(v) & v == 0L] = NA_integer_
    if(any(!is.na(v) & (v < 1L | v > J)))
      stop(sprintf("`%s` must hold integer codes in 1..J (0 or NA marks a missing measure).", nm))
    v
  }
  Y  = lapply(seq_len(k), function(j) code_col(tab[[ycol[j]]], ycol[j]))
  Xv = code_col(tab$X, "X")
  Xf = if(has_f) code_col(tab$Xf, "Xf") else NULL

  key_parts = c(list(Xv), if(has_f) list(Xf), Y)
  key = do.call(paste, c(lapply(key_parts, function(v)
    ifelse(is.na(v), "NA", as.character(v))), list(sep = "\r")))
  if(any(duplicated(key)))
    stop("`tab` has duplicated measure patterns; aggregate `n` over identical cells first.")

  # ---- transition horizons ----
  if(is.null(dstep)) dstep = rep(1L, k - 1L)
  dstep = as.integer(dstep)
  if(length(dstep) == 1L) dstep = rep(dstep, k - 1L)
  if(length(dstep) != k - 1L)
    stop("`dstep` must have length `k - 1` (one gap between consecutive measures).")
  if(any(is.na(dstep)) || any(dstep < 0L))
    stop("`dstep` must be non-negative integers.")
  dstep_f = as.integer(dstep_f)
  if(length(dstep_f) != 1L || is.na(dstep_f) || dstep_f < 0L)
    stop("`dstep_f` must be a single non-negative integer.")

  # ---- rival laws (DATA) ----
  as_rho = function(r, nm){
    if(is.null(r)) stop(sprintf("`%s` is required.", nm))
    if(is.null(dim(r))){
      if(length(r) != J) stop(sprintf("`%s` as a vector must have length `J`.", nm))
      r = matrix(as.numeric(r), J, J, byrow = TRUE)
    }
    r = as.matrix(r)
    if(!identical(dim(r), c(J, J))) stop(sprintf("`%s` must be a J x J matrix.", nm))
    if(any(r < 0) || any(abs(rowSums(r) - 1) > 1e-6))
      stop(sprintf("`%s` must have non-negative rows summing to one.", nm))
    r
  }
  if(is.list(rho) && !is.data.frame(rho)){
    if(length(rho) != k) stop("`rho` as a list must have one entry per linked measure.")
    rho = lapply(seq_len(k), function(j) as_rho(rho[[j]], sprintf("rho[[%d]]", j)))
  } else {
    rho = rep(list(as_rho(rho, "rho")), k)
  }
  if(has_f) rho_f = as_rho(rho_f, "rho_f")

  # ---- confusion kernel and mu groups ----
  use_mu = !is.null(kernel)
  if(use_mu){
    kernel = as.matrix(kernel)
    if(!identical(dim(kernel), c(J, J))) stop("`kernel` must be a J x J matrix.")
    if(any(kernel < 0) || any(abs(rowSums(kernel) - 1) > 1e-6))
      stop("`kernel` must have non-negative rows summing to one.")
  }
  nslot = k + 1L + as.integer(has_f)
  if(is.null(mu_index)) mu_index = rep(1L, nslot)
  mu_index = as.integer(mu_index)
  if(length(mu_index) != nslot)
    stop(sprintf("`mu_index` must have length %d (Y1..Yk, X%s).", nslot,
                 if(has_f) ", Xf" else ""))
  if(any(is.na(mu_index)) || any(mu_index < 1L) ||
     !identical(sort(unique(mu_index)), seq_len(max(mu_index))))
    stop("`mu_index` must use consecutive group labels starting at 1.")
  Gmu = max(mu_index)

  #------------------------------------------------------------
  # Freezes
  #------------------------------------------------------------

  fix_pis = NULL; fix_b = rep(NA_real_, k)
  if(!is.null(pi_shared_fixed)){
    if(!is.null(alpha_fixed) || !is.null(s_fixed))
      stop("supply either `pi_shared_fixed` or the (`alpha_fixed`, `s_fixed`) pair, not both.")
    if(length(pi_shared_fixed) != 1L || is.na(pi_shared_fixed) ||
       pi_shared_fixed < 0 || pi_shared_fixed >= 1)
      stop("`pi_shared_fixed` must be a single value in [0, 1).")
    fix_pis = as.numeric(pi_shared_fixed)
  }
  if(!is.null(alpha_fixed) || !is.null(s_fixed)){
    if(is.null(alpha_fixed) || is.null(s_fixed))
      stop("`alpha_fixed` and `s_fixed` must be supplied together: the slab is parameterised by (pi_shared, b), and only the pair (alpha_1, s) pins it in closed form.")
    if(length(alpha_fixed) == 1L) alpha_fixed = rep(alpha_fixed, k)
    if(length(alpha_fixed) != k)
      stop("`alpha_fixed` must have length `k` (NA for entries left free).")
    if(is.na(alpha_fixed[1]))
      stop("`alpha_fixed[1]` must be supplied: pi_shared = s_fixed * alpha_fixed[1].")
    if(length(s_fixed) != 1L || is.na(s_fixed) || s_fixed < 0 || s_fixed > 1)
      stop("`s_fixed` must be a single value in [0, 1].")
    fix_pis = s_fixed * alpha_fixed[1]
    ok = !is.na(alpha_fixed)
    if(any(alpha_fixed[ok] < fix_pis - 1e-12))
      stop("each fixed alpha must be at least pi_shared = s_fixed * alpha_fixed[1].")
    fix_b[ok] = (alpha_fixed[ok] - fix_pis) / (1 - fix_pis)
  }
  if(!is.null(alpha_f_fixed)){
    if(length(alpha_f_fixed) != 1L || is.na(alpha_f_fixed) ||
       alpha_f_fixed < 0 || alpha_f_fixed >= 1)
      stop("`alpha_f_fixed` must be a single value in [0, 1).")
  }
  fix_mu = rep(NA_real_, Gmu)
  if(!is.null(mu_fixed)){
    if(!use_mu && any(!is.na(mu_fixed) & mu_fixed != 0))
      stop("`mu_fixed` is non-zero but no `kernel` was supplied.")
    if(length(mu_fixed) == 1L) mu_fixed = rep(mu_fixed, Gmu)
    if(length(mu_fixed) != Gmu)
      stop("`mu_fixed` must have one entry per `mu_index` group.")
    if(any(!is.na(mu_fixed) & (mu_fixed < 0 | mu_fixed >= 1)))
      stop("`mu_fixed` entries must lie in [0, 1).")
    fix_mu = as.numeric(mu_fixed)
  }
  if(!use_mu) fix_mu = rep(0, Gmu)
  if(!is.null(T_fixed)){
    T_fixed = as.matrix(T_fixed)
    if(!identical(dim(T_fixed), c(J, J))) stop("`T_fixed` must be a J x J matrix.")
    if(any(T_fixed < 0) || any(abs(rowSums(T_fixed) - 1) > 1e-6))
      stop("`T_fixed` must have non-negative rows summing to one.")
  }
  if(!is.null(Pi_fixed)){
    Pi_fixed = as.matrix(Pi_fixed)
    if(!identical(dim(Pi_fixed), c(J, J))) stop("`Pi_fixed` must be a J x J matrix.")
    if(any(Pi_fixed < 0) || abs(sum(Pi_fixed) - 1) > 1e-6)
      stop("`Pi_fixed` must be non-negative and sum to one.")
  }
  est_T = is.null(T_fixed); est_Pi = is.null(Pi_fixed)
  est_pis = is.null(fix_pis); est_af = has_f && is.null(alpha_f_fixed)

  #------------------------------------------------------------
  # Cells, chunks, and the per-slot observation flags
  #------------------------------------------------------------

  obsY = lapply(Y, function(v) !is.na(v))
  safeY = lapply(Y, function(v){ v[is.na(v)] = 1L; v })
  obsX = !is.na(Xv); safeX = Xv; safeX[!obsX] = 1L
  if(has_f){ obsXf = !is.na(Xf); safeXf = Xf; safeXf[!obsXf] = 1L }

  nobs_slot = numeric(nslot)
  for(j in seq_len(k)) nobs_slot[j] = sum(n[obsY[[j]]])
  nobs_slot[k + 1L] = sum(n[obsX])
  if(has_f) nobs_slot[k + 2L] = sum(n[obsXf])

  if(is.null(chunk_size) || chunk_size >= nc){
    chunks = list(seq_len(nc))
  } else {
    chunk_size = max(1L, as.integer(chunk_size))
    chunks = split(seq_len(nc), ceiling(seq_len(nc) / chunk_size))
  }

  #------------------------------------------------------------
  # Initial values
  #------------------------------------------------------------

  Ident = diag(J)
  Pi0 = matrix(0, J, J)
  for(j in seq_len(k)){
    keep = obsX & obsY[[j]]
    if(!any(keep)) next
    tmp = rowsum(n[keep], safeX[keep] + J * (safeY[[j]][keep] - 1L))
    Pi0[as.integer(rownames(tmp))] = Pi0[as.integer(rownames(tmp))] + tmp
  }
  Pi0 = Pi0 + 1e-8 * max(Pi0, 1)
  Pi0 = Pi0 / sum(Pi0)

  T0 = 0.7 * Ident + 0.3 / J; T0 = T0 / rowSums(T0)

  getinit = function(nm, default) if(!is.null(init[[nm]])) init[[nm]] else default
  p = list(
    mu  = { m0 = getinit("mu", 0.10); if(length(m0) == 1L) rep(m0, Gmu) else as.numeric(m0) },
    pis = if(!est_pis) fix_pis else getinit("s", 0.30) * alpha_0,
    b   = rep(alpha_0, k),
    af  = if(!est_af) (if(has_f) alpha_f_fixed else 0) else getinit("alpha_f", alpha_0),
    T   = if(est_T) getinit("T", T0) else T_fixed,
    Pi  = if(est_Pi) getinit("Pi", Pi0) else Pi_fixed)
  if(length(p$mu) != Gmu) stop("`init$mu` must have one entry per `mu_index` group.")
  p$mu[!is.na(fix_mu)] = fix_mu[!is.na(fix_mu)]
  # start the free b's where the implied marginal false-link rate is alpha_0
  ok = !is.na(fix_b)
  p$b[!ok] = min(max((alpha_0 - p$pis) / (1 - p$pis), 1e-6), 0.99)
  p$b[ok] = fix_b[ok]
  p$T = p$T / rowSums(p$T)
  p$Pi = p$Pi / sum(p$Pi)

  #------------------------------------------------------------
  # One EM step: accumulate sufficient statistics chunk by chunk
  #------------------------------------------------------------

  branch_w = function(pp){
    w = numeric(2L^k)
    for(mask in seq_len(2L^k) - 1L){
      v = 1 - pp$pis
      for(j in seq_len(k))
        v = v * (if(bitwAnd(mask, bitwShiftL(1L, j - 1L)) != 0L) 1 - pp$b[j] else pp$b[j])
      w[mask + 1L] = v
    }
    w
  }

  emis_cols = function(Mt, idx, obs, zero_missing = FALSE){
    E = Mt[, idx, drop = FALSE]
    if(any(!obs)) E[, !obs] = if(zero_missing) 0 else 1
    E
  }

  chunk_stats = function(pp, ci, need_counts){
    Cc = length(ci)
    nn = n[ci]
    Mg = lapply(seq_len(Gmu), function(g)
      if(use_mu) (1 - pp$mu[g]) * Ident + pp$mu[g] * kernel else Ident)
    MgK = if(use_mu) lapply(seq_len(Gmu), function(g) pp$mu[g] * kernel) else NULL
    Ms = lapply(seq_len(nslot), function(t) Mg[[mu_index[t]]])
    MsK = if(use_mu) lapply(seq_len(nslot), function(t) MgK[[mu_index[t]]]) else NULL

    obsm = lapply(seq_len(k), function(j) obsY[[j]][ci])
    mown = lapply(seq_len(k), function(j)
      emis_cols(Ms[[j]], safeY[[j]][ci], obsm[[j]]))
    RM = lapply(seq_len(k), function(j) rho[[j]] %*% Ms[[j]])
    efresh = lapply(seq_len(k), function(j)
      emis_cols(RM[[j]], safeY[[j]][ci], obsm[[j]]))

    obsxc = obsX[ci]
    fb = emis_cols(Ms[[k + 1L]], safeX[ci], obsxc)
    if(has_f){
      obsfc = obsXf[ci]
      Tf = if(dstep_f == 0L) Ident else mat_power(pp$T, dstep_f)
      TMf = Tf %*% Ms[[k + 2L]]
      GMf = rho_f %*% Ms[[k + 2L]]
      Lnk = (1 - pp$af) * TMf + pp$af * GMf
      fl = emis_cols(Lnk, safeXf[ci], obsfc)
      flf = emis_cols(pp$af * GMf, safeXf[ci], obsfc)
      if(any(!obsfc)) flf[, !obsfc] = pp$af
    } else {
      obsfc = rep(FALSE, Cc); fl = 1; flf = NULL
    }
    Fm = if(has_f) fb * fl else fb

    keep_mid = need_counts && est_T
    tr = traj_backward(pp$T, dstep, mown, obsm, k, keep_mid)
    bw = tr$bw; bo = tr$bo; mid = tr$mid
    pi0 = rowSums(pp$Pi)
    wC = branch_w(pp)
    sf = traj_son(pp$Pi, pi0, mown, efresh, bw, wC, pp$pis, rho[[1L]], k, J, Cc,
                  keep = need_counts)
    Son = sf$Son

    P = colSums(Fm * Son)
    P = pmax(P, 1e-300)
    ll = sum(nn * log(P))
    if(!need_counts) return(list(ll = ll))

    d = nn / P
    u = Fm * rep(d, each = J)               # u[i, c] = d_c * Father[i, c]

    R_C = numeric(2L^k)
    Pi_count = matrix(0, J, J)
    T_count = matrix(0, J, J)

    for(mask in seq_len(2L^k) - 1L){
      b1 = bw[[1L]][[mask + 1L]]
      A = sf$A[[mask + 1L]]
      Phi = sf$Phi[[mask + 1L]]
      V = if(is.null(Phi)) u else u * Phi   # V[i, c] = d_c F[i,c] Phi_C[i,c]
      R_C[mask + 1L] = wC[mask + 1L] * sum(A * V)
      if(est_Pi){
        Pi_count = Pi_count + wC[mask + 1L] *
          (if(is.null(b1)) pp$Pi * rowSums(V) else pp$Pi * tcrossprod(V, b1))
      }
      if(est_T && k > 1L && wC[mask + 1L] > 0){
        fw = wC[mask + 1L] * crossprod(pp$Pi, V)          # forward at measure 1
        for(j in seq_len(k - 1L)){
          if(bitwAnd(mask, bitwShiftL(1L, j - 1L)) != 0L) fw = fw * mown[[j]]
          q = bitwShiftR(mask, j)                          # pattern for j+1..k
          info = bo[[j + 1L]][[q + 1L]]
          chain = mid[[j]][[q + 1L]]
          cur = fw
          if(dstep[j] > 0L){
            for(t in seq_len(dstep[j])){
              Bnext = chain[[t + 1L]]
              if(!is.null(Bnext) && any(info)){
                cw = if(all(info)) cur else cur * rep(info, each = J)
                T_count = T_count + pp$T * tcrossprod(cw, Bnext)
              }
              cur = crossprod(pp$T, cur)
            }
          }
          fw = cur
        }
      }
    }

    # shared-rival branch
    mass_s = u * (pp$pis * (pi0 * sf$Rv))
    R_SSS = sum(mass_s)
    if(est_Pi) Pi_count = Pi_count + pp$Pi * (rowSums(mass_s) / pmax(pi0, 1e-300))
    if(est_T && k > 1L && pp$pis > 0){
      fw = pp$pis * crossprod(rho[[1L]], u * pi0)
      for(j in seq_len(k - 1L)){
        fw = fw * mown[[j]]
        q = 2L^(k - j) - 1L
        info = bo[[j + 1L]][[q + 1L]]
        chain = mid[[j]][[q + 1L]]
        cur = fw
        if(dstep[j] > 0L){
          for(t in seq_len(dstep[j])){
            Bnext = chain[[t + 1L]]
            if(!is.null(Bnext) && any(info)){
              cw = if(all(info)) cur else cur * rep(info, each = J)
              T_count = T_count + pp$T * tcrossprod(cw, Bnext)
            }
            cur = crossprod(pp$T, cur)
          }
        }
        fw = cur
      }
    }

    # anchor branch
    R_f = 0
    if(has_f){
      R_f = sum(fb * flf * Son * rep(d, each = J))
      if(est_T && dstep_f > 0L){
        bwf = emis_cols(Ms[[k + 2L]], safeXf[ci], obsfc)
        fchain = traj_tchain(pp$T, dstep_f, bwf)
        cur = (1 - pp$af) * fb * Son * rep(d, each = J)
        if(any(!obsfc)) cur[, !obsfc] = 0
        for(t in seq_len(dstep_f)){
          T_count = T_count + pp$T * tcrossprod(cur, fchain[[t + 1L]])
          cur = crossprod(pp$T, cur)
        }
      }
    }

    # mu numerators via the emission-swap identity
    mu_num = numeric(Gmu)
    if(use_mu && any(is.na(fix_mu))){
      for(j in seq_len(k)){
        g = mu_index[j]
        if(!is.na(fix_mu[g])) next
        mown2 = mown; efresh2 = efresh
        mown2[[j]] = emis_cols(MsK[[j]], safeY[[j]][ci], obsm[[j]], TRUE)
        efresh2[[j]] = emis_cols(rho[[j]] %*% MsK[[j]], safeY[[j]][ci], obsm[[j]], TRUE)
        bw2 = traj_backward_swap(pp$T, dstep, mown2, k, j, bw)
        Son2 = traj_son(pp$Pi, pi0, mown2, efresh2, bw2, wC, pp$pis,
                        rho[[1L]], k, J, Cc)$Son
        mu_num[g] = mu_num[g] + sum(u * Son2)
      }
      gx = mu_index[k + 1L]
      if(is.na(fix_mu[gx])){
        fbK = emis_cols(MsK[[k + 1L]], safeX[ci], obsxc, TRUE)
        FmK = if(has_f) fbK * fl else fbK
        mu_num[gx] = mu_num[gx] + sum(FmK * Son * rep(d, each = J))
      }
      if(has_f){
        gf = mu_index[k + 2L]
        if(is.na(fix_mu[gf])){
          LnkK = (1 - pp$af) * (if(dstep_f == 0L) MsK[[k + 2L]] else
            mat_power(pp$T, dstep_f) %*% MsK[[k + 2L]]) +
            pp$af * (rho_f %*% MsK[[k + 2L]])
          flK = emis_cols(LnkK, safeXf[ci], obsfc, TRUE)
          mu_num[gf] = mu_num[gf] + sum(fb * flK * Son * rep(d, each = J))
        }
      }
    }

    list(ll = ll, R_C = R_C, R_SSS = R_SSS, R_f = R_f, mu_num = mu_num,
         Pi_count = Pi_count, T_count = T_count)
  }

  em_step = function(pp, need_counts = TRUE){
    ll = 0
    R_C = numeric(2L^k); R_SSS = 0; R_f = 0
    mu_num = numeric(Gmu)
    Pi_count = matrix(0, J, J); T_count = matrix(0, J, J)
    for(ci in chunks){
      st = chunk_stats(pp, ci, need_counts)
      ll = ll + st$ll
      if(!need_counts) next
      R_C = R_C + st$R_C; R_SSS = R_SSS + st$R_SSS; R_f = R_f + st$R_f
      mu_num = mu_num + st$mu_num
      Pi_count = Pi_count + st$Pi_count; T_count = T_count + st$T_count
    }
    if(!need_counts) return(list(ll = ll, par = pp))

    q = pp
    clip = function(x, lo = 1e-10, hi = 1 - 1e-6) pmin(pmax(x, lo), hi)

    if(use_mu){
      for(g in seq_len(Gmu)){
        if(!is.na(fix_mu[g])) next
        den = sum(nobs_slot[mu_index == g])
        q$mu[g] = if(den > 0) clip(mu_num[g] / den, 1e-9, 0.95) else pp$mu[g]
      }
    }
    if(est_pis) q$pis = clip(R_SSS / N, 1e-10, 0.95)
    ind = N - R_SSS
    for(j in seq_len(k)){
      if(!is.na(fix_b[j])) next
      fail_j = sum(R_C[bitwAnd(seq_len(2L^k) - 1L, bitwShiftL(1L, j - 1L)) == 0L])
      q$b[j] = if(ind > 0) clip(fail_j / ind, 1e-10, 0.999) else pp$b[j]
    }
    if(est_af) q$af = clip(R_f / N, 1e-10, 0.95)
    if(est_Pi) q$Pi = Pi_count / sum(Pi_count)
    if(est_T){
      rs = rowSums(T_count); keep = rs > 1e-9
      Tn = pp$T; Tn[keep, ] = T_count[keep, ] / rs[keep]; q$T = Tn
    }
    list(ll = ll, par = q,
         counts = list(R_C = R_C, R_SSS = R_SSS, R_f = R_f, N = N))
  }

  #------------------------------------------------------------
  # SQUAREM helpers
  #------------------------------------------------------------

  free_b = which(is.na(fix_b))
  free_mu = if(use_mu) which(is.na(fix_mu)) else integer(0)
  pack = function(pp){
    c(pp$mu[free_mu], if(est_pis) pp$pis, pp$b[free_b], if(est_af) pp$af,
      if(est_T) as.numeric(pp$T), if(est_Pi) as.numeric(pp$Pi))
  }
  unpack = function(v, pp){
    i = 1L
    if(length(free_mu)){ pp$mu[free_mu] = v[i:(i + length(free_mu) - 1L)]
                         i = i + length(free_mu) }
    if(est_pis){ pp$pis = v[i]; i = i + 1L }
    if(length(free_b)){ pp$b[free_b] = v[i:(i + length(free_b) - 1L)]
                        i = i + length(free_b) }
    if(est_af){ pp$af = v[i]; i = i + 1L }
    if(est_T){ pp$T = matrix(v[i:(i + J * J - 1L)], J, J); i = i + J * J }
    if(est_Pi) pp$Pi = matrix(v[i:(i + J * J - 1L)], J, J)
    pp
  }
  project = function(pp){
    cl = function(x, lo = 1e-10, hi = 0.95) pmin(pmax(x, lo), hi)
    if(length(free_mu)) pp$mu[free_mu] = cl(pp$mu[free_mu], 1e-9, 0.95)
    if(est_pis) pp$pis = cl(pp$pis)
    if(length(free_b)) pp$b[free_b] = cl(pp$b[free_b], 1e-10, 0.999)
    if(est_af) pp$af = cl(pp$af)
    if(est_T){ pp$T = pmax(pp$T, 1e-12); pp$T = pp$T / rowSums(pp$T) }
    if(est_Pi){ pp$Pi = pmax(pp$Pi, 1e-14); pp$Pi = pp$Pi / sum(pp$Pi) }
    pp
  }

  #------------------------------------------------------------
  # EM iterations
  #------------------------------------------------------------

  loglik_trace = numeric(0); ll_old = -Inf; converged = FALSE; nfe = 0L

  for(iter in seq_len(maxit)){
    s1 = em_step(p); nfe = nfe + 1L
    ll = s1$ll
    loglik_trace = c(loglik_trace, ll)
    if(verbose && iter %% 25 == 0) cat("iter", iter, "loglik", ll, "\n")
    if(is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + 1)){
      converged = TRUE
      break
    }
    ll_old = ll
    p1 = s1$par
    if(!squarem){ p = p1; next }
    s2 = em_step(p1); nfe = nfe + 1L
    p2 = s2$par
    v0 = pack(p); v1 = pack(p1); v2 = pack(p2)
    r = v1 - v0; uu = v2 - v1 - r
    nu = sqrt(sum(uu * uu))
    if(!is.finite(nu) || nu < 1e-14){ p = p2; next }
    alp = -sqrt(sum(r * r)) / nu
    ok = FALSE
    for(bt in seq_len(6)){
      pn = project(unpack(v0 - 2 * alp * r + alp * alp * uu, p))
      sn = em_step(pn); nfe = nfe + 1L
      if(is.finite(sn$ll) && sn$ll >= s2$ll){ p = sn$par; ok = TRUE; break }
      alp = (alp - 1) / 2
      if(alp > -1.0001) break
    }
    if(!ok) p = p2
  }

  if(!converged)
    warning("EM did not converge within `maxit` iterations; increase `maxit` or loosen `tol`.")

  alpha = p$pis + (1 - p$pis) * p$b
  if(any(alpha > 0.5))
    warning("An estimated false-link rate exceeds 0.5: the diagonal-dominance identification condition fails at the estimate. Interpret with caution.")

  names(alpha) = paste0("alpha", seq_len(k))
  bvec = p$b; names(bvec) = paste0("b", seq_len(k))

  # One extra E-pass at the final parameters, exposing the expected-count
  # aggregates that the shared-parameter M-steps are ratios of. This is
  # what lets misclassifyr_traj_em_shared() pool (pi_shared, b) across
  # views exactly rather than approximately.
  sfin = em_step(p, need_counts = TRUE)

  return(list(
    alpha = alpha,
    s = if(alpha[1] > 0) unname(p$pis / alpha[1]) else 0,
    pi_shared = p$pis,
    b = bvec,
    alpha_f = if(has_f) p$af else NA_real_,
    mu = if(use_mu) p$mu else rep(0, Gmu),
    T = p$T,
    Pi = p$Pi,
    k = k,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    n_iter = length(loglik_trace),
    n_feval = nfe,
    converged = converged,
    Ecounts = c(sfin$counts, list(ll = sfin$ll))
  ))
}
