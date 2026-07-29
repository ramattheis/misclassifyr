# misclassifyr

<!-- badges: start -->
[![R-CMD-check](https://github.com/ramattheis/misclassifyr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ramattheis/misclassifyr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Estimation and inference for misclassification models, as described in
Mattheis (2024).

## What it does

You have linked historical records, and some of the links are wrong. A wrong
link does not put a random number in your outcome column — it puts somebody
else's outcome there. Downstream regression estimates are attenuated by
roughly the link failure rate, and no amount of clustering or robust
standard errors will fix it.

`misclassifyr` models the error and corrects the estimate. Given

- a correctly measured discrete regressor $X$,
- **two** noisy measures $Y_1, Y_2$ of a discrete latent outcome $Y^*$, and
- optional discrete controls $W$,

it estimates, by maximum likelihood:

- $\Pi$, the joint distribution of $X$ and the latent outcome $Y^*$;
- $\Delta$, the misclassification distribution $(Y_1, Y_2) \mid Y^*$; and
- $\beta$, the coefficient from the linear projection of $Y^*$ on $X$, with
  delta-method standard errors and Monte Carlo confidence sets that stay
  valid under partial identification.

The two-measures requirement is the essential one. With one noisy measure
you cannot tell a genuinely low outcome from a mis-linked one. With two,
their pattern of disagreement identifies the error rates, and the error
rates identify everything else.

## Installation

```r
# install.packages("devtools")
devtools::install_github("ramattheis/misclassifyr")
```

## A worked example

The packaged `ancienregime` data holds 100,000 synthetic father–son pairs
with the son's income observed twice, in 1770 and 1780.

```r
library(misclassifyr)
data(ancienregime)

inputs <- prep_misclassification_data(
  data = ancienregime,
  outcome_1 = "son_income_1780",  outcome_1_bin = "son_occupation_1780",
  outcome_2 = "son_income_1770",  outcome_2_bin = "son_occupation_1770",
  regressor = "father_income_1750", regressor_bin = "father_occupation_1750",
  weights = "linked_weight", record_vals = TRUE, round_vals = 0
)

fit <- misclassifyr(
  tab = inputs$tab, J = inputs$J, K = inputs$K,
  X_names = inputs$X_names, Y1_names = inputs$Y1_names, Y2_names = inputs$Y2_names,
  X_vals = inputs$X_vals, Y_vals = inputs$Y_vals,
  model_to_Delta = model_to_Delta_RL_ind, makeplots = FALSE
)

Pi_to_beta(X_vals = inputs$X_vals, Y_vals = inputs$Y_vals,
           Pi_mle = fit$Pi_hat_mle, cov_Pi = fit$cov_Pi_mle)$beta_hat_mle
#> MLE beta
#>    0.724
```

The naive slope on the same data and the same discretization is 0.586: the
uncorrected estimate of intergenerational persistence is about 19% too low.

## Vignettes

```r
vignette("getting-started", package = "misclassifyr")
```

- **[Getting started](vignettes/getting-started.Rmd)** — the workflow end to
  end on the packaged data, what the output means, and how much the
  correction moves the answer.
- **[Designing misclassification models](vignettes/designing-misclassification-models.Rmd)**
  — the main reference. How to enumerate the error channels in your own
  setting (transcription, coding drift, linkage failure, heaping), turn them
  into a parameterized $\Delta$, write it as an R closure with a matching
  log prior, and verify it on simulated data. Includes worked prompts for
  using a large language model as a drafting assistant, and — more
  importantly — how to check what it produces. Also covers stacked
  estimation with a shared error rate across conditioning cells.
- **[High-dimensional outcomes](vignettes/sparse-and-large.Rmd)** — sparse
  tabulations and a dedicated expectation-maximization estimator for the
  record-linkage model when the outcome has hundreds or thousands of
  categories, including two-step estimation and separating true transitions
  from linkage error.
- **[Inference](vignettes/inference.Rmd)** — delta-method standard errors,
  what a singular information matrix means (boundary versus weak
  identification), Chen–Christensen–Tamer Monte Carlo confidence sets, and
  which to report.

## Main functions

| Function | Purpose |
|---|---|
| `prep_misclassification_data()` | Tabulate microdata into the form the estimators expect; `sparse = TRUE` for high-dimensional outcomes |
| `misclassifyr()` | General maximum likelihood estimation of $\Pi$ and $\Delta$, with an optional Gibbs sampler for posterior inference |
| `Pi_to_beta()` | Map $\hat\Pi$ to a regression coefficient with standard errors and confidence sets |
| `model_to_Delta_RL_ind()`, `model_to_Delta_NP_ind()`, `model_to_Delta_NP()` | Built-in shapes for the misclassification distribution |
| `make_empirical_Delta_RL*()` | Factories that plug the empirical margin into a record-linkage design |
| `misclassifyr_rl_em()` | Fast EM estimation of the record-linkage model when the outcome has thousands of categories |
| `misclassifyr_rl_em_stacked()` | The same, with the error rate shared across conditioning cells |
| `misclassifyr_known_slab()` | The false-link rate for a single link judged against a known phantom distribution |
| `misclassifyr_traj_em()` | EM for multi-link bundles: one shared latent rival, a Markov latent state, and a measurement-error layer |
| `synthetic_data()` | Simulate from the model to check a design before using it |

## Citation

```r
citation("misclassifyr")
```

> Mattheis, R. (2024). *Misclassification models for linked historical
> records.* Working paper.

## Status

This package is under active development. No promises are made about
backwards compatibility yet.

## License

MIT © Ross Mattheis
