#' Maps model parameters, psi, to Delta, the distribution of Y1, Y2 | Y* under record linkage error structure.
#'
#' @param psi A numeric vector of length 2(`J`-1)+2`J` corresponding to the column and row scales of the record linkage.
#' @return A numeric vector of length `J`^3 corresponding to the values of the
#'   `J`x`J`^2 matrix `Delta`, built from single-measure matrices of the
#'   record-linkage "slab" form `(1 - alpha) I + rho alpha'`.
#' @examples
#' J <- 3
#' alpha <- 0.2
#' # Uniform draw distribution (first 2 * (J - 1) entries) and a common
#' # false-link rate alpha for both measures (last 2 * J entries)
#' psi <- c(rep(0, J - 1), rep(0, J - 1),
#'          rep(qlogis(alpha), J), rep(qlogis(alpha), J))
#' D <- matrix(model_to_Delta_RL_ind(psi), nrow = J)
#' rowSums(D)  # all ones
#'
#' # The implied single-measure matrix: (1 - alpha) on the diagonal plus a
#' # uniform slab of total mass alpha in every column
#' round(diag(J) * (1 - alpha) + outer(rep(1 / J, J), rep(alpha, J)), 3)
#' @export
model_to_Delta_RL_ind = function(psi){

  # J is a deterministic function of psi for model_to_Delta_RL_ind
  J = as.integer((length(psi) + 2)/4)

  # Extracting the row scales for Delta^{(1)} and Delta^{(2)}
  psi1 = c(psi[1:(J-1)],0) # Adding the reference value
  psi2 = c(psi[((J-1)+1):(2*(J-1))],0) # Adding the reference value
  row1 = exp(psi1)/sum(exp(psi1)) # logit link
  row2 = exp(psi2)/sum(exp(psi2)) # logit link

  # Extracting the column scales for Delta^{(1)} and Delta^{(2)}
  psi = psi[-(1:(2*(J-1)))]
  psi1 = psi[1:J]
  psi2 = psi[(J+1):(2*J)]
  col1 = exp(psi1)/(1 + exp(psi1)) # Column margins are probabilities but need not sum to one
  col2 = exp(psi2)/(1 + exp(psi2)) # Column margins are probabilities but need not sum to one

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


