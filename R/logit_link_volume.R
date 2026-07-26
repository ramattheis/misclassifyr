#' Evaluates the log of the determinant of the Jacobian of the logit link.
#'
#' @param x A numeric vector in logit-space corresponding to some .
#' @return A numeric value equal to the log of the determinant of the Jacobian
#'   of the logit transform (plus a constant that does not depend on `x`).
#' @examples
#' # The volume factor that turns a flat prior on the simplex into the
#' # corresponding prior in unconstrained logit space
#' logit_link_volume(c(0, 0))
#'
#' # It equals the numerical log |det| of the softmax-with-reference map
#' softmax_free <- function(x) { p <- exp(c(x, 0)); p <- p / sum(p); p[-length(p)] }
#' x <- c(-1.2, 0.4)
#' log(abs(det(numDeriv::jacobian(softmax_free, x))))
#' logit_link_volume(x) - lfactorial(3)
#' @keywords internal
#' @export
logit_link_volume = function(x){

  # Applying the logit link
  p = exp(x)/(1+sum(exp(x)))

  # Building the Jacobian
  # The off-diagonal elements of the Jacobian of the logit transform are -p_i p_j
  # The diagonal elements are p_i(1-p_i)
  J_ = -1*p %*% t(p)
  diag(J_) = p*(1-p)

  # Computing the log of the determinant of the Jacobian of the logit transform
  # the log of the determinant is the sum of the log of the eigenvalues
  log_det_logit = sum(softlog(eigen(J_)$values))

  # Adding the flat prior constant (just for kicks)
  # the density in an N-simplex is 1/N!
  # Because the constant does not depend on the parameter (only its dimension),
  # it's folded into the constant divisor of the posterior
  flat_prior_constant = lfactorial(length(p) + 1)

  # Returning the log difference: log(|J_logit|) - log(1/N!)
  return(log_det_logit + flat_prior_constant)
}
