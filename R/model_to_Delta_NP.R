#' Maps model parameters, psi, to the conditional distribution Y1, Y2 | Y*, Delta.
#'
#' @param psi A numeric vector containing Pi, Delta^{(1)}, and Delta^{(2)}.
#' @return A numeric vector corresponding to the (JxJ)xJ matrix Delta.
#' @export
model_to_Delta_NP = function(psi){

  # Importing J from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  if(!exists("J", envir = misclassifyr_env)){stop("Error: `J` missing from `misclassifyr_env`")}
  J = misclassifyr_env$J


  # Exponentiating to return to levels (psi in logs for numerical performance)
  psi = exp(psi)

  # Building Delta
  Delta = matrix(psi, nrow = J)
  # The last column ensures rows of Delta sum to 1
  Delta = cbind(Delta, 1-apply(Delta,1, sum))

  return(c(Delta))

}

# Adding a name as an attribute
attr(model_to_Delta_NP, "name") = "model_to_Delta_NP"

