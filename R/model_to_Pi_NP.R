#' Maps model parameters, phi, to the joint distribution of X and Y*, Pi.
#'
#' @param phi A numeric vector.
#' @param J An integer corresponding to the dimension of Y.
#' @param ... Additional, optional arguments.
#' @return A numeric vector corresponding to the JxK matrix Pi
#' @export
model_to_Pi_NP = function(phi,J,...){

  # Building Pi
  phi = c(exp(phi)/(1 + sum(exp(phi))), 1/(1 + sum(exp(phi)))) # Logit link
  Pi = matrix(phi,nrow = J)   # Converting to matrix, Y* rows, X cols

  return(c(Pi))

}

# Adding a name as an attribute
attr(model_to_Pi_NP, "name") = "model_to_Pi_NP"


