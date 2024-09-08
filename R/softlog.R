#' Logarithm with a lower bound
#'
#' @param x A numeric vector.
#' @return A numeric vector composed of the elements of \code{log(x)} or \code{log(1e-20)} for element is less than \code{1e-20}.
#' @keywords internal
#' @export
softlog = function(x) {
  return(log(pmax(x,1e-20)) )
}
