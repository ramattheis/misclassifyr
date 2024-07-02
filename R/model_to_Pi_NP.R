#' Maps model parameters, phi, to the joint distribution of X and Y*, Pi.
#'
#' Longer description of what it does...
#'
#'
#' @param phi A numeric vector.
#' @param J An integer corresponding to the dimension of Y.
#' @param K An integer corresponding to the dimension of X.
#' @return something
#' @details some details.
#' @examples
#' \dontrun{
#' some example code # Should return something
#' }
#' @export
model_to_Pi_NP = function(phi){

  # Importing J, K from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  J = misclassifyr_env$J

  # Building Pi
  phi = exp(phi)              # Exponentiating to return to levels
  phi = c(phi, 1-sum(phi))    # The last entry forced probabilities to sum to 1
  Pi = matrix(phi,nrow = J)   # Converting to matrix, Y* rows, X cols

  return(c(Pi))
}

# Adding a name as an attribute
attr(model_to_Pi_NP, "name") = "model_to_Pi_NP"


