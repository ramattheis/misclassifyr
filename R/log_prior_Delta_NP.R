#' Evaluates the log of the prior of Delta at model_to_Delta_NP(psi).
#'
#' @param psi A numeric vector parameterizing `Delta` through `model_to_Delta_NP`.
#' @return A numeric value equal to the log of the flat prior of `Delta` at `psi`, re-scaled for the logit transform.
#' @export
log_prior_Delta_NP = function(psi){

  # J is a deterministic function of psi for model_to_Delta_NP
  # Using Cardano's formula to solve for J...
  cardano_discriminant = (length(psi)/2)^2  - 1/27
  cardano_u1 = length(psi)/2 + sqrt(cardano_discriminant)
  cardano_u2 = length(psi)/2 - sqrt(cardano_discriminant)
  J = as.integer(round(cardano_u1^(1/3) + cardano_u2^(1/3),0)) # rounding because of floating point weirdness

  # Building Delta (without the last entry)
  Delta_ = matrix(psi, nrow = J)

  # Computing the log prior in logit space equivalent to a flat prior in probability space
  log_prior_Delta = sum(apply(Delta_,1, logit_link_volume))

  return(log_prior_Delta)

}
