#' Computes the standard error of beta as function of the joint distribution of X and Y*, Pi
#'
#' @import Matrix
#'
#' @param Pi A numeric vector or list of numeric vectors containing the elements of Pi.
#' @param cov_Pi A numeric vector or a list of numeric vectors representing the covariance of estimates of the elements of Pi.
#' @param X_vals A numeric vector or a list of numeric vectors representing the scalar values associated with X.
#' @param Y_vals A numeric vector or a list of numeric vectors representing the scalar values associated with Y.
#' @param W_weights A numeric vector representing the sample size of each control cell.
#' @return A scalar equal to the standard error of beta.
se_beta_deltamethod = function(Pi, cov_Pi, X_vals, Y_vals, W_weights){

  # Is Pi a list?
  if(class(Pi) == "list"){

    # Recording the length of the list
    n_control_cells = length(Pi)

    # If Pi is a list, defining a wrapper function in order to take the gradient w.r.t. the full set of Pis
    Pi_to_beta_wrapper = function(Pis){

      # Forming Pis back into a list
      Pi_list = split(Pis, rep(1:n_control_cells, each = length(Pis)/n_control_cells))

      # Estimating Pi
      return(Pi_to_beta(Pi_list,X_vals,Y_vals,W_weights))
    }

    # Forming the list of Pis into one long vector
    Pis = c(unlist(Pi))

    # Computing the gradient of beta with respect to Pi
    grad_Pi_to_beta = numDeriv::grad(Pi_to_beta_wrapper, Pis)

    # Forming a block-diagonal matrix of cov_Pi
    cov_Pi = Matrix::bdiag(cov_Pi)

    # Computing the variance of beta via the delta method
    var_beta = as.numeric(t(grad_Pi_to_beta) %*% cov_Pi %*% grad_Pi_to_beta)

    # Returning the standard error of beta
    se_beta = sqrt(var_beta)

  } else {

    # If Pi is not a list, computing the standard error directly via the delta method

    # Defining a wrapper for the Pi_to_beta function to compute the gradient
    Pi_to_beta_wrapper = function(Pi_){ return(Pi_to_beta(Pi_,X_vals,Y_vals,1)) }

    # Computing the gradient of beta with respect to Pi
    grad_Pi_to_beta = numDeriv::grad(Pi_to_beta_wrapper, Pi)

    # Computing the variance of beta via the delta method
    var_beta = t(grad_Pi_to_beta) %*% cov_Pi %*% grad_Pi_to_beta

    # Returning the standard error of beta
    se_beta = sqrt(var_beta)
  }

  # Making sure se_beta is numeric (as opposed to a matrix)
  se_beta = as.numeric(se_beta)

  # Returning the standard error of beta
  return(se_beta)

}
