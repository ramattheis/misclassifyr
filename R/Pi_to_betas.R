#' Maps the joint distribution, Pi, of X and Y* to a vector, representing beta in each covariate cell
#'
#' Longer description of what it does...
#'
#' @param Pi A list of numeric vectors containing the elements of Pi.
#' @param X_vals A list of numeric vectors representing the scalar values associated with X.
#' @param Y_vals A list of numeric vectors representing the scalar values associated with Y.
#' @return A scalar equal to beta.
#' @details some details.
#' @examples
#' \dontrun{
#' some example code # Should return something
#' }
#' @export
Pi_to_betas = function(Pi, X_vals, Y_vals){

  # Returning an error if Pi is not a list
  if(!is.list(Pi)){stop("Error: Pi should be a list in `Pi_to_betas()`. Use `Pi_to_beta()` if there are no controls.")}

  # Importing J from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  if(!exists("J", envir = misclassifyr_env)){stop("Error: `J` missing from `misclassifyr_env`")}
  J = misclassifyr_env$J
  if(!exists("K", envir = misclassifyr_env)){stop("Error: `K` missing from `misclassifyr_env`")}
  K = misclassifyr_env$K

  #------------------------------------------------------------
  # Defining a function to compute moments X, Y, XY, and XX w/n cells
  #------------------------------------------------------------

  make_beta_from_moments = function(Pi_, X_vals_, Y_vals_){
    # Converting Pi_ to a matrix
    Pi_ = matrix(Pi_, nrow = J)

    # Computing the expected value of X, Y, XX, and XY
    E_X = sum(Pi_ %*% diag(X_vals_))
    E_Y = sum(diag(Y_vals_) %*% Pi_)
    E_XX = sum(Pi_ %*% diag(X_vals_^2))
    E_XY = sum(diag(Y_vals_) %*% Pi_ %*% diag(X_vals_))

    # Compute beta
    beta = (E_XY - E_X*E_Y)/(E_XX-E_X^2)

    # Returning moments as a row of a dataframe
    return(beta)
  }


  #------------------------------------------------------------
  # Computing and returning vector of betas
  #------------------------------------------------------------

  # Computing moments within each cell
  betas = sapply(seq_along(Pi), function(j) make_beta_from_moments(Pi[[j]],X_vals[[j]],Y_vals[[j]]))

  # Returning betas
  return(betas)

}

