#' Evaluates the log of the prior of Delta at model_to_Delta_NP_ind(psi).
#'
#' @param psi A numeric vector parameterizing `Delta` through `model_to_Delta_NP_ind`.
#' @return A numeric value equal to the log of the flat prior of `Delta` at `psi`, re-scaled for the logit transform.
#' @export
log_prior_Delta_NP_ind = function(psi){

  # J is a deterministic function of psi for model_to_Delta_NP_ind
  J = as.integer((1 + sqrt(1 + 2 * length(psi))) / 2)

  # Building Delta^{(1)} and Delta^{(2)}
  Delta = matrix(psi, nrow = J-1)
  Delta1 = Delta[,1:J]
  Delta2 = Delta[,(J+1):(2*J)]

  # Computing the log prior in logit space equivalent to a flat prior in probability space
  log_prior_Delta1 = sum(apply(Delta1,2, logit_link_volume))
  log_prior_Delta2 = sum(apply(Delta1,2, logit_link_volume))

  return(log_prior_Delta1 + log_prior_Delta2)

}
