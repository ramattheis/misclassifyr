#' Estimate a false-link rate against a known phantom distribution
#'
#' A two-component mixture estimator for the cleanest identification
#' available in record linkage: a single link whose correctness can be
#' judged against an outcome whose population distribution is KNOWN.
#'
#' The setting. For each linked pair we observe a category `y` reported by
#' the linked record and a category `x` observed independently of the link
#' (so `x` is correct whatever the link did). If the link is correct,
#' \eqn{y} follows a transition distribution \eqn{T(y \mid x)} — for
#' example, true migration between the two dates. If the link is
#' incorrect, the linked record belongs to somebody else, so \eqn{y} is a
#' draw from the population joint of the linked file, \eqn{\rho(y \mid
#' \cdot)}, independent of \eqn{x}. Both \eqn{T} and \eqn{\rho} are
#' supplied by the user from external data; only the false-link rate
#' \eqn{\alpha} is estimated:
#' \deqn{P(y \mid x) = (1-\alpha)\, T(y \mid x) + \alpha\, \rho(y).}
#'
#' The canonical application is the 1940 census migration question: the
#' 1930 record fixes a man's county in 1930 independently of any link,
#' while the linked 1940 record reports where its person lived in 1935.
#' Coherent pairs reflect true 1930-35 migration; incoherent ones are
#' draws from the population, which the full count reveals exactly. Unlike
#' the repeated-measures estimators in this package, nothing here needs a
#' second noisy measure — one link suffices, because the phantom
#' distribution is observed rather than inferred.
#'
#' Estimation is by direct maximization of a one-parameter likelihood, so
#' the standard error is exact rather than asymptotic in a nuisance
#' dimension.
#'
#' @param tab A data.frame with columns `x`, `y` (integer category codes)
#'   and `n` (counts of linked pairs).
#' @param T_cond A data.frame with columns `x`, `y`, `t`: the conditional
#'   distribution of the reported category given the link-independent
#'   category when the link is CORRECT. Must sum to one within `x` (rows
#'   for unobserved `(x, y)` pairs may be omitted and are treated as
#'   zero).
#' @param rho A numeric vector, or a data.frame with columns `y` and `p`,
#'   giving the population distribution of the reported category among
#'   records the linker could have hit. Must sum to one.
#' @param alpha_bounds Length-2 numeric, the search interval for
#'   \eqn{\alpha}.
#' @return A list with `alpha`, `se` (from the observed information),
#'   `ci` (95% Wald interval, clipped to \[0, 1\]), `loglik`, `n`, and
#'   `share_impossible`: the share of linked pairs whose reported category
#'   has zero probability under `T`, which are necessarily false links and
#'   provide a model-free lower bound on \eqn{\alpha}.
#' @examples
#' # Truth: 20% of links are false; correct links move "one step" w.p. 0.3
#' set.seed(7)
#' J <- 10; N <- 20000; alpha <- 0.2
#' x <- sample.int(J, N, replace = TRUE)
#' rho_true <- rep(1/J, J)
#' y_correct <- ifelse(runif(N) < 0.7, x, (x %% J) + 1)
#' y <- ifelse(runif(N) < alpha, sample.int(J, N, replace = TRUE), y_correct)
#' tab <- aggregate(list(n = rep(1, N)), by = list(x = x, y = y), FUN = sum)
#'
#' T_cond <- do.call(rbind, lapply(1:J, function(i)
#'   data.frame(x = i, y = c(i, (i %% J) + 1), t = c(0.7, 0.3))))
#'
#' fit <- misclassifyr_known_slab(tab, T_cond, rho_true)
#' fit$alpha
#' fit$ci
#' fit$share_impossible   # a model-free lower bound on alpha
#' @export
misclassifyr_known_slab = function(tab, T_cond, rho,
                                   alpha_bounds = c(1e-6, 1 - 1e-6)){

  #------------------------------------------------------------
  # Input checks
  #------------------------------------------------------------

  if(!all(c("x", "y", "n") %in% colnames(tab)))
    stop("`tab` should have columns `x`, `y`, and `n`.")
  if(!all(c("x", "y", "t") %in% colnames(T_cond)))
    stop("`T_cond` should have columns `x`, `y`, and `t`.")

  tab = tab[tab$n > 0, , drop = FALSE]
  x = as.integer(tab$x); y = as.integer(tab$y); n = as.numeric(tab$n)
  if(any(is.na(x)) || any(is.na(y)) || any(x < 1) || any(y < 1))
    stop("`x` and `y` must be positive integer codes.")

  if(is.data.frame(rho)){
    if(!all(c("y", "p") %in% colnames(rho)))
      stop("If `rho` is a data.frame it should have columns `y` and `p`.")
    rv = numeric(max(c(y, rho$y)))
    rv[as.integer(rho$y)] = as.numeric(rho$p)
    rho = rv
  }
  if(abs(sum(rho) - 1) > 1e-6) stop("`rho` must sum to one.")
  if(any(rho < 0)) stop("`rho` must be non-negative.")
  if(max(y) > length(rho)) stop("`rho` is shorter than the largest `y` code.")

  ts = rowsum(as.numeric(T_cond$t), as.integer(T_cond$x))
  if(any(abs(ts - 1) > 1e-6))
    stop("`T_cond` rows must sum to one within `x`.")

  #------------------------------------------------------------
  # Per-cell T and rho lookups
  #------------------------------------------------------------

  key   = function(a, b) a + (b - 1) * (max(c(x, T_cond$x)) + 1)
  t_map = as.numeric(T_cond$t)
  names(t_map) = key(as.integer(T_cond$x), as.integer(T_cond$y))
  t_cell = t_map[as.character(key(x, y))]
  t_cell[is.na(t_cell)] = 0
  r_cell = rho[y]

  # Cells impossible under a correct link are necessarily false links; the
  # share of them is a model-free lower bound on alpha.
  share_impossible = sum(n[t_cell == 0]) / sum(n)

  #------------------------------------------------------------
  # One-parameter likelihood
  #------------------------------------------------------------

  negll = function(a) -sum(n * log((1 - a) * t_cell + a * r_cell))

  opt = stats::optimize(negll, interval = alpha_bounds, tol = 1e-10)
  alpha = opt$minimum

  # Observed information: -d2/da2 of the log likelihood, in closed form
  d = r_cell - t_cell
  p = (1 - alpha) * t_cell + alpha * r_cell
  info = sum(n * (d / p)^2)
  se = if(info > 0) sqrt(1 / info) else NA_real_

  list(alpha = alpha,
       se = se,
       ci = c(max(0, alpha - 1.96 * se), min(1, alpha + 1.96 * se)),
       loglik = -opt$objective,
       n = sum(n),
       share_impossible = share_impossible)
}
