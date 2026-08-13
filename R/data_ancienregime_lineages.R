#' Synthetic father-son lineages from Ancien Regime France: three links,
#' a fallible clerk, and the same wrong man
#'
#' A companion to [ancienregime] for the two error processes that a
#' plain record-linkage model cannot see. Father-son pairs of the Third
#' Estate: the father's occupation is recorded in the 1750 capitation
#' and again -- through a link -- in the 1745 dixieme; the son is linked
#' into the rolls of 1770, 1780, and 1790, so his occupation is measured
#' three times. Three linked measures are exactly what
#' [misclassifyr_traj_em()] needs, and this dataset is built to its
#' model:
#'
#' * **Dependent link failures.** With probability `pi_shared`, *every*
#'   one of a son's failed links lands on the same rival -- a specific
#'   other man whose own occupation follows the same decade-to-decade
#'   transition process. Otherwise failures are independent draws from
#'   the occupation slab. Two measures are provably blind to this
#'   dependence; the third identifies it.
#' * **Clerk miscoding.** Independently of linkage, every recorded
#'   occupation -- the father's, the son's, a rival's -- is miscoded to
#'   an adjacent rung of the occupational ladder with probability `mu`.
#'   A two-measure model has nowhere to put this error except the
#'   false-link rate; the trajectory model separates it via the
#'   measurement kernel.
#'
#' Occupations use the six classes of the `ancienregime` world, ordered
#' from lowest to highest standing: `"Vagabond"`, `"Metayer"`,
#' `"Journalier"`, `"Petit Metiers"`, `"Petite Bourgeoisie"`,
#' `"Haute Bourgeoisie"`. The generating parameters -- including the
#' transmission matrix \eqn{\Pi} that a mobility researcher is after --
#' are in [ancienregime_truth].
#'
#' @format A data frame with 100,000 rows and 6 variables:
#' \describe{
#'   \item{province}{The province of the lineage (flavor; the error
#'     processes here do not vary by province).}
#'   \item{father_occupation_1750}{The father's occupation as recorded
#'     in the 1750 capitation: his true class, save for clerk miscoding.}
#'   \item{father_occupation_1745_linked}{The father's occupation on the
#'     1745 dixieme record the clerks linked him to -- a second, fallible
#'     reading of the father.}
#'   \item{son_occupation_1770_linked}{The son's occupation on the
#'     linked 1770 roll.}
#'   \item{son_occupation_1780_linked}{The son's occupation on the
#'     linked 1780 roll.}
#'   \item{son_occupation_1790_linked}{The son's occupation on the
#'     linked 1790 roll.}
#' }
#' @source Synthetic data generated for the package by
#'   `data-raw/make_ancienregime_extras.R`.
#' @seealso [ancienregime_truth] for the generating parameters,
#'   [ancienregime_parishes] for the sparse-outcome setting, and
#'   `vignette("tour-ancien-regime")` for the worked example.
#' @examples
#' data(ancienregime_lineages)
#' str(ancienregime_lineages)
#'
#' # The two earliest measures of the son disagree more often than clerk
#' # error alone could explain -- the surplus is linkage:
#' mean(ancienregime_lineages$son_occupation_1770_linked !=
#'      ancienregime_lineages$son_occupation_1780_linked)
"ancienregime_lineages"
