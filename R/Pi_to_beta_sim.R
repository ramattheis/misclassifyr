#' Maps the joint distribution, Pi, of X and Y* to a scalar, beta via simulation
#'
#' Longer description of what it does...
#'
#' @param Pi A list of numeric vectors containing the elements of Pi
#' @param X_vals A numeric vector or a list of numeric vectors representing the scalar values associated with X.
#' @param Y_vals A numeric vector or a list of numeric vectors representing the scalar values associated with Y.
#' @param W_weights A numeric vector representing the sample size of each control cell.
#' @return A scalar equal to beta.
#' @details some details.
#' @examples
#' \dontrun{
#' some example code # Should return something
#' }
#' @export
Pi_to_beta_sim = function(Pi, X_vals, Y_vals, W_weights){

  # Importing J from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  if(!exists("J", envir = misclassifyr_env)){stop("Error: `J` missing from `misclassifyr_env`")}
  J = misclassifyr_env$J
  if(!exists("K", envir = misclassifyr_env)){stop("Error: `K` missing from `misclassifyr_env`")}
  K = misclassifyr_env$K

  #------------------------------------------------------------
  # Defining a function to generate synthetic data
  #------------------------------------------------------------

  make_synthetic = function(Pi_, X_vals_, Y_vals_, W_weight = 1){

    # Generate some tabulated data
    X = c(sapply(Pi_vals_X, function(x) rep(x, J)))
    Y = rep(Y_vals_,K)
    W = rmultinom(1,size=1e8,prob=Pi_) * W_weight

    # Returning moments as a row of a dataframe
    out = as.data.frame(cbind(X,Y,W))
    colnames(out) = c("X","Y","W")
    return(out)
  }


  #------------------------------------------------------------
  # Computing and returning beta
  #------------------------------------------------------------

  if(is.list(Pi)){
    # Generating synthetic data
    synthetic = do.call(rbind, lapply(seq_along(Pi), function(j) make_synthetic(Pi[[j]],X_vals[[j]],Y_vals[[j]],W_weights[j])))
  } else {
    # Generating synthetic data
    synthetic = make_synthetic(Pi,X_vals,Y_vals)
  }

  # Returning beta
  return(unname(lm(Y~X,weights = W,synthetic)$coefficients[2]))
}
