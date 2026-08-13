# cran-comments

## Submission

This is a new submission. misclassifyr 0.3.1 is the first version submitted
to CRAN. Relative to the internal 0.3.0, it adds two documentation
datasets and a comprehensive vignette; the check below is for 0.3.1.

## Test environments

* Local: macOS 26.5 (aarch64-apple-darwin23), R 4.6.0 --
  `R CMD check --as-cran`
* GitHub Actions (`.github/workflows/R-CMD-check.yaml`):
  * macos-latest, R release
  * ubuntu-latest, R release

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE

## The remaining NOTE

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ross Mattheis <ramattheis@gmail.com>'

New submission
```

This is the expected note for a package that is not yet on CRAN. The
maintainer address is correct and reachable.

## Notes on content

* The `Description` field cites Mattheis (2024), the working paper the
  methods come from. No DOI is given because the paper is not yet published;
  the reference will be updated to the `authors (year) <doi:...>` form as
  soon as one is assigned.

* Examples, tests, and vignettes all run without downloading anything or
  writing outside `tempdir()`. Total example CPU time is under two seconds,
  with no single example above half a second; the longest-running pieces of
  the package (the MCMC sampler and the high-dimensional EM estimator) are
  demonstrated in vignettes at sizes that keep each vignette under 20
  seconds.

* The package has no compiled code.
