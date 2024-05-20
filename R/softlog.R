#' Logarithm with a lower bound
#'
#' If \code{x} is greater than \code{1e-20}, it returns \code{log(x)}.
#' Otherwise, it returns \code{log(1e-20)}.
#'
#' @param x A numeric vector.
#' @return A numeric vector composed of the elements of \code{log(x)} or \code{log(1e-20)} for element is less than \code{1e-20}.
#' @details The arguments to the \code{ll} function are the log value of the probabilities in Delta and Pi of the misclassification model. To prevent convergence issues, \code{ll} enforces a lower bound on these probabilities of \code{1e-20}. This function is useful for mapping probabilities that may be zero to \code{ll}. It will throw an error if any element of \code{x} is negative.
#' @examples
#' \dontrun{
#' softlog(c(0.5, 0.1, 0))  # Should return the log values including log(1e-20) for 0
#' softlog(c(0.5, 0.1, -0.1))   # Should throw an error
#' }
#' @export
softlog <- function(x) {
  if (any(x < 0)) {
    stop("Error: softlog attempt to compute log of a negative value.")
  }
  return(ifelse(x > 1e-20, log(x), log(1e-20)))
}
