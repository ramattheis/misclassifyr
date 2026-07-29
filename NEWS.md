# misclassifyr 0.3.0

First release prepared for CRAN. Everything below is new since 0.2.3.

## New features

* `misclassifyr_rl_em()`: a fast expectation-maximization estimator for the
  common-alpha record-linkage model, for outcomes with thousands of
  categories, where the general-purpose `misclassifyr()` machinery is
  infeasible. The mixture likelihood factorizes into four linkage
  configurations per observed cell, so each iteration is
  `O(observed cells + support of Pi)` and neither the `J` x `J`^2
  misclassification matrix nor the balanced tabulation is ever formed.
  Optional arguments fix the failed-link draw distributions (`rho1`,
  `rho2`) at known population margins, fix `(alpha1, alpha2)` at values
  estimated elsewhere (`alpha_fixed`, the second step of the two-step
  architecture for high-dimensional outcomes), and supply a true-transition
  operator for the second measure (`T2`), which separates genuine change
  between the two observation dates from linkage error.
* `misclassifyr_rl_em_stacked()`: the same estimator over a list of
  tabulations split by a conditioning variable, with the false-link rates
  shared across cells and `rho` and `Pi` cell-specific. This is the
  conditional-independence version of the model: phantom draws for a failed
  link come from the candidate pool matching the linking keys, not the
  unconditional marginal. An unconditional slab that is too diffuse makes
  false links look correct and biases `alpha` downward.
* `misclassifyr_traj_em()`: the trajectory EM for multi-link bundles, and
  the estimator behind the three-link design. A unit contributes an anchor
  measure observed directly, optionally a second anchor measure reached by a
  link of its own, and `k` linked measures of a second person whose latent
  state follows a Markov chain `T`. It relaxes the two assumptions the
  rank-one model in `misclassifyr_rl_em()` makes about failed links. First,
  failures are *structured*: each unit has one latent rival, and with
  probability `s` every failed link lands on that same rival, whose own
  state follows `T` — which is why two false links agree far more often
  than the rank-one model allows, and why `s` needs `k >= 3` (at `k = 2`
  `alpha` and `s` trade off along a ridge). Second, the recorded category
  need not be the true one: every observed value passes through
  `(1 - mu) I + mu K` for a supplied confusion kernel `K`, with `mu`
  allowed to vary by census year. `Pi`, `T`, `mu`, the false-link rates and
  the shared share are estimated jointly, and every M-step is a closed-form
  count ratio — no Newton step anywhere — provided the phantom laws `rho`
  are supplied as data rather than rebuilt from `T` each iteration.
  Freezes (`alpha_fixed` with `s_fixed`, `pi_shared_fixed`, `mu_fixed`,
  `T_fixed`, `Pi_fixed`, `alpha_f_fixed`) make each staged hand-off the
  same function restricted, and `pi_shared_fixed` is the profiling handle
  for `s`. A missing measure is coded `NA` (or `0`) and handled as a
  missing emission: the unit contributes the branches its observed pattern
  supports. Do not model it as an extra latent state — a high-mass
  "not observed" category is agreed on by two links for reasons unrelated
  to linkage, and the model reads that agreement as correct linking.
* `prep_misclassification_data(sparse = TRUE)` returns a tabulation holding
  only observed cells, plus a `cell_idx` column giving each cell's position
  in the balanced `(Y2, Y1, X)` layout. `misclassifyr()` accepts either
  form and returns identical estimates. The balanced table has `J`^2 * `K`
  rows, which is impractical well before the outcome gets interesting.
* `misclassifyr()` exposes `lambda_dd`, the weight on the diagonal-dominance
  penalty for `Delta`. It defaults to `sum(tab$n)^2` (a near-hard
  constraint) and can be lowered, or set to 0, to weaken or remove the
  restriction. Starting values are now boundary-safe.
* `Pi_to_beta()` gains `ll_history`, the MCMC log-posterior history returned
  by `misclassifyr(bayesian = TRUE)`, which it uses to build the
  Chen-Christensen-Tamer Monte Carlo confidence set. The set stays valid
  under partial identification, where the delta-method interval does not.

## Bug fixes

* `make_empirical_Delta_RL_common_alpha_mixed_NP()`'s single-tabulation
  branch errored with "incorrect number of dimensions" on every call: a
  stray `c()` flattened `Delta2` to a vector before the reference row was
  dropped. The list branch was already correct and is the reference
  implementation.
* `Pi_to_beta()`'s `bayesian = TRUE` branch now works: it previously
  referenced objects that were never passed in (`tab`,
  `misclassification_output`). The Chen-Christensen-Tamer HPD interval now
  takes the MCMC likelihood history through the new `ll_history` argument
  (structure-checked against `posterior_Pi`), and with fewer than 20
  posterior draws the 5% cutoff clamps to the minimum instead of returning
  an empty interval.
* The Metropolis-Hastings acceptance rule in `misclassifyr()`'s Gibbs
  sampler no longer subtracts the proposal increment (`gibbs_jump`): the
  proposal is a symmetric random walk in the unconstrained parameter space,
  so no Hastings correction belongs in the ratio. Chains from earlier
  versions were biased toward negative jumps.
* `misclassifyr()` now rejects a burn-in that the thinning rate does not
  divide -- the posterior draw keys and the likelihood history would
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

## Documentation

* Four vignettes: `getting-started` (the workflow end to end on the
  packaged `ancienregime` data, with the naive/corrected comparison),
  `designing-misclassification-models` (error channels, parameterized
  forms for `Delta`, writing the closure and its matching log prior,
  verification harness, worked prompts for using a large language model
  as a drafting assistant, and stacked shared-alpha estimation),
  `sparse-and-large` (sparse tabulations and `misclassifyr_rl_em()` at
  `J = 300`, `alpha_fixed`, `T2`), and `inference` (delta method,
  singular information, Chen-Christensen-Tamer Monte Carlo sets).
* Runnable `@examples` added to every exported function, plus the
  `ancienregime` dataset, which was previously undocumented (its roxygen
  block was missing the trailing object name, so no Rd was generated).
* `misclassifyr()`'s `@return` block documented only `Pi_hat_MLE`; it now
  documents all 24 returned components. Its non-existent `split_eta`
  argument is no longer documented.
* `softlog()` is exported but was marked `@noRd`, so `R CMD check`
  reported it as undocumented; it now has an Rd page.
* `make_empirical_Delta_RL_common_alpha_mixed_NP()` gained its missing
  `@param J`, and `Pi_to_beta()` its missing `@return`.
* Added `inst/CITATION` and a package-level help page; rewrote `README.md`.

## CRAN readiness

* `DESCRIPTION`: the `License` field was the literal string
  `` `use_mit_license()` ``; it is now `MIT + file LICENSE` with the
  standard template `LICENSE` and full text in `LICENSE.md`. Title no
  longer ends in a period, the description is expanded, `URL` and
  `BugReports` added, the unused `Rcpp` import dropped, and `stats` and
  `utils` declared with matching imports.
* Column names referred to non-standardly (inside `dplyr` verbs, `ggplot2`
  aesthetics, and `subset()` calls) are declared with
  `utils::globalVariables()`, clearing 40-odd "no visible binding for
  global variable" notes.
* `Pi_to_beta()` passed `weight =` to `stats::lm()`, which partially
  matched `weights`. Now spelled in full; behaviour is unchanged.
* `se_beta_deltamethod()` defined `Pi_to_beta_inner_wrapper` twice in the
  same function with different formal arguments; the single-cell branch's
  wrapper is now named distinctly.
* `misclassifyr()`'s `bayesian` and `mle` argument checks use
  `is.logical()` rather than comparing `class()` to a string.
* `data/ancienregime.rda` is re-saved with bzip2 compression, 758Kb to
  454Kb. The stored data frame is unchanged.
* Added `cran-comments.md`.

## Testing

* A full testthat suite (250+ tests) covering tabulation invariants, exact
  likelihood/prior/transform values, seeded estimator regression tests, the
  sparse and dense paths against each other, the EM estimators, the
  Bayesian sampler, and an end-to-end workflow on the packaged
  `ancienregime` data. CI via GitHub Actions (R CMD check).
* `misclassifyr_traj_em()` is tested against the simulation it was
  validated on: recovery of `alpha`, `s`, `mu`, `T` and `Pi` at the truth;
  the `k = 2` versus `k = 3` identification contrast (the profiled
  likelihood in the shared share is flat at `k = 2` and has an interior
  bowl at `k = 3`); missing-emission handling, including that `NA` and `0`
  are the same marker; the freeze options as no-ops when frozen at the free
  estimates; monotonicity of the log likelihood under plain EM; and that
  chunking changes nothing but working memory.
