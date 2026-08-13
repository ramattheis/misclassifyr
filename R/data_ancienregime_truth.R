#' The generating truth behind the Ancien Regime datasets
#'
#' The parameters used to generate [ancienregime_parishes] and
#' [ancienregime_lineages], together with the parish register and the
#' auxiliary tables an analyst would bring to the linked data. Shipping
#' the truth keeps the worked examples honest: every estimate in
#' `vignette("tour-ancien-regime")` is compared against the value that
#' generated the data, and the package's tests assert the comparison.
#' Real data never comes with this object; that is the point of
#' synthetic data.
#'
#' @format A named list:
#' \describe{
#'   \item{occupations}{The six occupational classes, lowest to highest.}
#'   \item{provinces}{The twelve provinces of the `ancienregime` world.}
#'   \item{parish_register}{A data frame mapping the integer parish codes
#'     of [ancienregime_parishes] to parish names, provinces, and 1750
#'     population shares.}
#'   \item{parishes}{The parish world's parameters: `alpha` (the two
#'     false-link rates), `T_decade` (the 240 x 240 decade-long
#'     parish-to-parish transition operator), and the latent parish
#'     margins `margin_1750`, `margin_1770`, `margin_1780`.}
#'   \item{vingtieme_flows}{A flow table `(j, l, t, n)` built from the
#'     (fictional) vingtieme of 1780, which asked every household where
#'     it resided in 1770: the decade operator as an analyst would
#'     observe it, sampling noise included. Pass `(j, l, t)` as the `T2`
#'     argument of [misclassifyr_rl_em()].}
#'   \item{lineages}{The lineage world's parameters: `alpha` (the three
#'     marginal false-link rates of the son's links), `s` and `pi_shared`
#'     and `b` (the dependence structure: `alpha = pi_shared +
#'     (1 - pi_shared) * b`), `alpha_f` (the father's link), `mu` (the
#'     clerk's miscoding share), `Pi` (the father-to-son transmission
#'     matrix -- the estimand), `T_occ` (the decade operator on
#'     occupations), `occupation_slab` (the latent slab failed links draw
#'     from), `father_margin`, and `kernel` (the clerk's adjacent-rung
#'     confusion kernel).}
#' }
#' @source Generated alongside the datasets by
#'   `data-raw/make_ancienregime_extras.R`.
#' @seealso [ancienregime_parishes], [ancienregime_lineages], and
#'   `vignette("tour-ancien-regime")`.
#' @examples
#' data(ancienregime_truth)
#' ancienregime_truth$lineages$alpha     # the rates the vignette recovers
#' head(ancienregime_truth$parish_register)
"ancienregime_truth"
