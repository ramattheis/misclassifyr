#' Maps model parameters, psi, to Delta, the fully non-parametric  distribution of Y1, Y2 | Y*
#'
#' @param psi A numeric vector of length`J`^2*(`J`-1) containing all but the last row of `Delta`.
#' @return A numeric vector corresponding to the `J`x`J`^2 matrix `Delta`.
#' @export
model_to_Delta_NP = function(psi){

  # J is a deterministic function of psi for model_to_Delta_NP
  # Using Cardano's formula to solve for J...
  cardano_discriminant = (length(psi)/2)^2  - 1/27
  cardano_u1 = length(psi)/2 + sqrt(cardano_discriminant)
  cardano_u2 = length(psi)/2 - sqrt(cardano_discriminant)
  J = as.integer(round(cardano_u1^(1/3) + cardano_u2^(1/3),0)) # rounding because of floating point weirdness

  # Building Delta Matrix
  Delta = matrix(psi, nrow = J)
  Delta = cbind(Delta, rep(0,J)) # Adding back the reference value
  Delta = t(apply(Delta,1,function(d) exp(d)/sum(exp(d)))) # logit link

  return(c(Delta))

}

# Adding a name as an attribute
attr(model_to_Delta_NP, "name") = "model_to_Delta_NP"

