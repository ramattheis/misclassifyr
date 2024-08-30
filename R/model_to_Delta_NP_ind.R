#' Maps model parameters, psi, to the conditional distribution Y1, Y2 | Y*, Delta.
#'
#' @param psi A numeric vector of length 2xJx(J-1) containing Delta^{(1)} and Delta^{(2)}.
#' @return A numeric vector corresponding to the (JxJ)xJ matrix Delta.
#' @keywords internal
#' @export
model_to_Delta_NP_ind = function(psi){

  # J is a deterministic function of psi for model_to_Delta_NP_ind
  J = as.integer((1 + sqrt(1 + 2 * length(psi))) / 2)

  # Exponentiating to return to levels (psi in logs for numerical performance)
  psi = exp(psi)

  # Building Delta^{(1)} and Delta^{(2)}
  Delta = matrix(psi, nrow = J-1)
  Delta1 = Delta[,1:J]
  Delta2 = Delta[,(J+1):(2*J)]
  # The last row ensures columns of Delta^{(1)} and Delta^{(2)} sum to 1
  Delta1 = rbind(Delta1, 1-apply(Delta1,2, sum))
  Delta2 = rbind(Delta2, 1-apply(Delta2,2, sum))

  # Building the joint distribution Y1, Y2 conditional on Y*
  Delta = lapply(1:nrow(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
  Delta = do.call(cbind, Delta)

  return(c(Delta))

}

# Adding a name as an attribute
attr(model_to_Delta_NP_ind, "name") = "model_to_Delta_NP_ind"


