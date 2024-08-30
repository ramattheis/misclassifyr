#' Computes variance of beta_hat
#'
#' @param Pi_cov A numeric vector or list of numeric vectors containing the elements of the covariance of Pi.
#' @param Pi A numeric vector or list of numeric vectors containing the elements of Pi.
#' @param X_vals A numeric vector or a list of numeric vectors representing the scalar values associated with X.
#' @param Y_vals A numeric vector or a list of numeric vectors representing the scalar values associated with Y.
#' @param W_weights A numeric vector representing the sample size of each control cell.
#' @return A scalar equal to beta.
#' @export
se_beta_deltamethod = function(Pi_cov, Pi, X_vals, Y_vals, W_weights = NA){




  # Returning SE of beta
  return(Cov_XY/Var_X)

}
