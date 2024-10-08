#' Maps the joint distribution, Pi, of X and Y* to a scalar, beta
#'
#' @param Pi A numeric vector or list of numeric vectors containing the elements of Pi.
#' @param X_vals A numeric vector or a list of numeric vectors representing the scalar values associated with X.
#' @param Y_vals A numeric vector or a list of numeric vectors representing the scalar values associated with Y.
#' @param W_weights A numeric vector representing the sample size of each control cell.
#' @return A scalar equal to beta.
#' @keywords internal
#' @export
Pi_to_beta_inner = function(Pi, X_vals, Y_vals, W_weights){

  #------------------------------------------------------------
  # Defining a function to compute moments X, Y, XY, and XX w/n cells
  #------------------------------------------------------------

  make_moments = function(Pi_, X_vals_, Y_vals_){

    # Making Pi_ non-negative
    Pi_[Pi_ < 0] = 0
    Pi_ = Pi_/sum(Pi_)

    # Converting Pi_ to a matrix
    Pi_ = matrix(Pi_, nrow = length(Y_vals_))

    # Computing the expected value of X, Y, XX, and XY
    E_X = sum(Pi_ %*% diag(X_vals_))
    E_Y = sum(diag(Y_vals_) %*% Pi_)
    E_XX = sum(Pi_ %*% diag(X_vals_^2))
    E_XY = sum(diag(Y_vals_) %*% Pi_ %*% diag(X_vals_))

    # Returning moments as a row of a dataframe
    return(as.data.frame(cbind(E_X,E_Y,E_XX,E_XY)))

  }


  #------------------------------------------------------------
  # Computing and returning beta
  #------------------------------------------------------------

  if(class(Pi) == "list"){

    # Computing moments within each cell
    moments = do.call(rbind, lapply(seq_along(Pi), function(j) make_moments(Pi[[j]],X_vals[[j]],Y_vals[[j]])))

    # Computing Var(X)
    Var_X = sum(W_weights * moments$E_XX) - sum((W_weights*moments$E_X))^2

    # Computing Cov(X,Y) by the law of total covariance
    Cov_XY = sum(W_weights* (moments$E_XY - moments$E_X*moments$E_Y)) +
      sum(moments$E_Y*moments$E_X*W_weights) - sum(moments$E_Y*W_weights)*sum(moments$E_X*W_weights)

  } else {

    # Computing moments for the full population
    moments = make_moments(Pi,X_vals,Y_vals)

    # Computing Var(X)
    Var_X = moments$E_XX - moments$E_X^2

    # Computing Cov(X,Y)
    Cov_XY = moments$E_XY - moments$E_X*moments$E_Y

  }

  # Returning beta
  return(Cov_XY/Var_X)

}

