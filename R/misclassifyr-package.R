#' @keywords internal
#' @importFrom stats dlogis lm median optim rexp rmultinom rnorm runif sd
#' @importFrom utils globalVariables
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Column names referred to non-standardly -- inside dplyr verbs, ggplot2
# aesthetics, and subset() calls -- rather than as ordinary R variables. R's
# code analysis cannot see that these always resolve to columns of a data frame
# built inside the package, so it reports them as "no visible binding for global
# variable". Declaring them here is the documented remedy (Writing R Extensions,
# "Checking and building packages").
#
#   X, Y, Y1, Y2                  tabulation columns
#   X_bin, Y1_bin, Y2_bin, weight columns built by prep_misclassification_data()
#   Pi_hat, Delta_hat, Ys         estimate columns in the plotting/posterior frames
#   X_name, Y_name, Y1_name,
#   Y2_name, Ys_name              label columns in the same frames
#   draw                          the MCMC draw index
#   object, type                  columns of misclassifyr()'s input_types frame
#
# Called unqualified, with a matching @importFrom above: this is a top-level
# call, so it runs at install time and is not retained in the namespace, and a
# `utils::` prefix here would leave R CMD check reporting utils as a declared
# but unused import.
globalVariables(c(
  "Delta_hat", "Pi_hat", "X", "X_bin", "X_name", "Y", "Y1", "Y1_bin",
  "Y1_name", "Y2", "Y2_bin", "Y2_name", "Y_name", "Ys", "Ys_name",
  "draw", "object", "type", "weight"
))
