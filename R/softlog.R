#' Logarithm with a lower bound
#'
#' A numerically safe logarithm used throughout the package's likelihoods:
#' probabilities that an optimizer pushes to (or slightly below) zero return
#' `log(1e-20)` instead of `-Inf` or `NaN`, which keeps `optim()` on a finite
#' objective surface.
#'
#' @param x A numeric vector.
#' @return A numeric vector composed of the elements of \code{log(x)}, or
#'   \code{log(1e-20)} for any element less than \code{1e-20}.
#' @examples
#' softlog(c(1, 0.5, 0))
#' log(c(1, 0.5, 0))  # -Inf in the last position
#' @keywords internal
#' @export
softlog = function(x) {
  return(log(pmax(x,1e-20)) )
}
