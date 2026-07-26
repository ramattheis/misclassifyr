#' Maps model parameters, phi, to the joint distribution of X and Y*, Pi.
#'
#' @param phi A numeric vector.
#' @param J An integer corresponding to the dimension of Y.
#' @param ... Additional, optional arguments.
#' @return A numeric vector of length `J * K` corresponding to the `J` by `K`
#'   matrix `Pi`, stored column-major: rows index \eqn{Y^*}, columns index \eqn{X}.
#' @examples
#' # phi = log(2, 3, 4) maps to probabilities (2, 3, 4, 1) / 10
#' model_to_Pi_NP(log(c(2, 3, 4)), J = 2)
#'
#' # Flat parameters give the uniform joint distribution
#' matrix(model_to_Pi_NP(rep(0, 3 * 3 - 1), J = 3), nrow = 3)
#' @export
model_to_Pi_NP = function(phi,J,...){

  # Building Pi
  phi = c(exp(phi)/(1 + sum(exp(phi))), 1/(1 + sum(exp(phi)))) # Logit link
  Pi = matrix(phi,nrow = J)   # Converting to matrix, Y* rows, X cols

  return(c(Pi))

}

# Adding a name as an attribute
attr(model_to_Pi_NP, "name") = "model_to_Pi_NP"


