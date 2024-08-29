#' Maps model parameters, phi, to the joint distribution of X and Y*, Pi.
#'
#' @param phi A numeric vector.
#' @param J An integer corresponding to the dimension of Y.
#' @param K An integer corresponding to the dimension of X.
#' @return A numeric vector corresponding to the JxK matrix Pi
#' @export
model_to_Pi_NP = function(phi){

  # Importing J, K from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  if(!exists("J", envir = misclassifyr_env)){stop("Error: `J` missing from `misclassifyr_env`")}
  J = misclassifyr_env$J
  if(!exists("K", envir = misclassifyr_env)){stop("Error: `K` missing from `misclassifyr_env`")}
  K = misclassifyr_env$K

  # Building Pi
  phi = exp(phi)  # Exponentiating to return to levels
  phi = c(phi[1:(K-1)], 1-sum(phi), phi[K:(J*K-1)] ) # Bottom-left corner forces probabilities to sum to 1
  Pi = matrix(phi,nrow = J)   # Converting to matrix, Y* rows, X cols

  return(c(Pi))

}

# Adding a name as an attribute
attr(model_to_Pi_NP, "name") = "model_to_Pi_NP"


