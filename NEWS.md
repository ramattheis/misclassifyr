# misclassifyr (development version)

## Bug fixes (2026-07-23)

* `Pi_to_beta()`'s `bayesian = TRUE` branch now works: it previously
  referenced objects that were never passed in (`tab`,
  `misclassification_output`). The Chen–Christensen–Tamer HPD interval now
  takes the MCMC likelihood history through the new `ll_history` argument
  (structure-checked against `posterior_Pi`), and with fewer than 20
  posterior draws the 5% cutoff clamps to the minimum instead of returning
  an empty interval.
* The Metropolis–Hastings acceptance rule in `misclassifyr()`'s Gibbs
  sampler no longer subtracts the proposal increment (`gibbs_jump`): the
  proposal is a symmetric random walk in the unconstrained parameter space,
  so no Hastings correction belongs in the ratio. Chains from earlier
  versions were biased toward negative jumps.
* `misclassifyr()` now rejects a burn-in that the thinning rate does not
  divide — the posterior draw keys and the likelihood history would
  otherwise be disjoint and downstream joins silently empty.
* `make_empirical_Delta_RL_common_alpha()`: the single-tabulation branch
  now matches the list branch (it built the record-linkage error term with
  `1 - alpha` in place of `alpha`) and its log prior no longer hits a
  `retun()` typo.
* `log_prior_Delta_NP_ind()` computed the prior for the second
  misclassification matrix from the first one's block of `psi`.
* Internal `class(x) == "list"` comparisons (which error for multi-class
  objects such as data.tables) replaced with a single `is_cell_list()`
  helper; the input-type bookkeeping in `misclassifyr()` uses the same
  predicate, so multi-class inputs no longer break the control-cell
  copying machinery.
* `misclassifyr()` checks `tab` for `NA`s before the count-range check,
  which previously failed with an unhelpful error when `n` contained `NA`.
* Removed stale arguments (`estimate_beta`, `Y_names`) from the `misc/`
  example scripts.

## Testing

* Restored a full testthat suite (200+ tests) covering tabulation
  invariants, exact likelihood/prior/transform values, seeded estimator
  regression tests, the Bayesian sampler, and an end-to-end workflow on
  the packaged `ancienregime` data. CI via GitHub Actions (R CMD check).
