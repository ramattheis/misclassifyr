#' Maps model parameters, psi, to the joint distribution of the data, theta.
#'
#' @param psi A numeric vector of length 2(J-1)+2J corresponding to the column and row scales of the record linkage.
#' @return something
#' @keywords internal
#' @export
model_to_Delta_RL_ind = function(psi){

  # J is a deterministic function of psi for model_to_Delta_RL_ind
  J = as.integer((length(psi) + 2)/4)

  # Exponentiating to return to levels (psi in logs for numerical performance)
  psi = exp(psi)

  # Extracting the row scales for Delta^{(1)} and Delta^{(2)}
  row1 = psi[1:(J-1)]
  row2 = psi[((J-1)+1):(2*(J-1))]
  row1 = c(row1, 1-sum(row1)) # The last element forces the sum to 1
  row2 = c(row2, 1-sum(row2)) # The last element forces the sum to 1

  # Extracting the column scales for Delta^{(1)} and Delta^{(2)}
  psi = psi[-(1:(2*(J-1)))]
  col1 = psi[1:J]
  col2 = psi[(J+1):(2*J)]

  # Building Delta^{(1)} and Delta^{(2)}
  Delta1 = diag(J)*(1-col1) + outer(row1,col1,"*")
  Delta2 = diag(J)*(1-col2) + outer(row2,col2,"*")

  # Building the joint distribution Y1, Y2 conditional on Y*
  Delta = lapply(1:ncol(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
  Delta = do.call(cbind, Delta)

  return(c(Delta))
}

# Adding a name as an attribute
attr(model_to_Delta_RL_ind, "name") = "model_to_Delta_RL_ind"


