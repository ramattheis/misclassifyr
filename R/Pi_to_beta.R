#' Maps the joint distribution, Pi, of X and Y* to a scalar, beta
#'
#' Longer description of what it does...
#'
#' @param Pi A list of numeric vectors containing the elements of Pi
#' @param Pi_vals A list of numeric vectors representing the scalar values associated with each row/column of Pi. The first J elements correspond to the values of Y* and the next K elements correspond to the values of X.
#' @param W_weights A numeric vector representing the sample size of each control cell.
#' @return A scalar equal to beta.
#' @details some details.
#' @examples
#' \dontrun{
#' some example code # Should return something
#' }
#' @export
Pi_to_beta = function(Pi, Pi_vals, W_weights){

  # Importing J from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  if(!exists("J", envir = misclassifyr_env)){stop("Error: `J` missing from `misclassifyr_env`")}
  J = misclassifyr_env$J

  # Replacing Pi_vals with a list if Pi is a list and Pi_vals is a vector
  if(is.list(Pi) & !is.list(Pi_vals)){Pi_vals = replicate(length(Pi), Pi_vals, simplify=F)}

  #------------------------------------------------------------
  # Defining a function to compute moments X, Y, XY, and XX w/n cells
  #------------------------------------------------------------

  make_moments = function(Pi_, Pi_vals_){

    # Splitting Pi_vals into X and Y components
    Pi_vals_Y = Pi_vals_[1:J]
    Pi_vals_X = Pi_vals_[(J+1):length(Pi_vals_)]

    # Converting Pi_ to a matrix
    Pi_ = matrix(Pi_, nrow = J)

    # Computing the expected value of X, Y, XX, and XY
    E_X = sum(Pi_ %*% diag(Pi_vals_X))
    E_Y = sum(diag(Pi_vals_Y) %*% Pi_)
    E_XX = sum(Pi_ %*% diag(Pi_vals_X^2))
    E_XY = sum(diag(Pi_vals_Y) %*% Pi_ %*% diag(Pi_vals_X))

    # Returning moments as a row of a dataframe
    return(as.data.frame(cbind(E_X,E_Y,E_XX,E_XY)))
  }


  #------------------------------------------------------------
  # Computing and returning beta
  #------------------------------------------------------------

  if(is.list(Pi)){
    # Computing moments within each cell
    moments = do.call(rbind, lapply(seq_along(Pi), function(j) make_moments(Pi[[j]],Pi_vals[[j]])))

    # Computing Var(X)
    Var_X = sum(W_weights * moments$E_XX) - sum((W_weights*moments$E_X))^2

    # Computing Cov(X,Y) by the law of total covariance
    Cov_XY = sum(W_weights* (moments$E_XY - moments$E_X*moments$E_Y)) +
      sum(moments$E_Y*moments$E_X*W_weights) - sum(moments$E_Y*W_weights)*sum(moments$E_X*W_weights)
  } else {
    # Computing moments for the full population
    moments = make_moments(Pi,Pi_vals)

    # Computing Var(X)
    Var_X = moments$E_XX - moments$E_X^2

    # Computing Cov(X,Y)
    Cov_XY = moments$E_XY - moments$E_X*moments$E_Y
  }

  # Returning beta
  return(Cov_XY/Var_X)

}

