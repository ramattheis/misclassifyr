#' Maps model parameters, psi, to Delta, the distribution of Y1, Y2 | Y* under conditional independence of Y1, Y2 on Y
#'
#' @param psi A numeric vector of length `2 * J * (J - 1)` holding the
#'   parameters of \eqn{\Delta^{(1)}} and \eqn{\Delta^{(2)}}, the two
#'   single-measure misclassification matrices, stacked in that order.
#' @return A numeric vector of length `J`^3 corresponding to the values of the
#'   `J`x`J`^2 matrix `Delta`. Column `(y2 - 1) * J + y1` of row `y*` holds
#'   `Pr(Y1 = y1 | Y* = y*) * Pr(Y2 = y2 | Y* = y*)`.
#' @examples
#' J <- 3
#' # Flat parameters give the uniform (fully uninformative) error distribution
#' D <- matrix(model_to_Delta_NP_ind(rep(0, 2 * J * (J - 1))), nrow = J)
#' dim(D)                     # 3 x 9
#' rowSums(D)                 # all ones
#'
#' # A diagonal-dominant design: recover Pr(Y1 | Y*) by marginalizing over Y2
#' set.seed(1)
#' D <- matrix(model_to_Delta_NP_ind(rnorm(2 * J * (J - 1))), nrow = J)
#' round(D %*% do.call(rbind, replicate(J, diag(J), simplify = FALSE)), 3)
#' @export
model_to_Delta_NP_ind = function(psi){

  # J is a deterministic function of psi for model_to_Delta_NP_ind
  J = as.integer((1 + sqrt(1 + 2 * length(psi))) / 2)

  # Building Delta^{(1)} and Delta^{(2)}
  Delta = matrix(psi, nrow = J-1)
  Delta = rbind(Delta, rep(0, ncol(Delta))) # Adding back the reference value
  Delta = apply(Delta,2,function(d) exp(d)/sum(exp(d)))
  Delta1 = Delta[,1:J]
  Delta2 = Delta[,(J+1):(2*J)]

  # Building the joint distribution Y1, Y2 conditional on Y*
  Delta = lapply(1:nrow(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
  Delta = do.call(cbind, Delta)

  return(c(Delta))

}

# Adding a name as an attribute
attr(model_to_Delta_NP_ind, "name") = "model_to_Delta_NP_ind"


