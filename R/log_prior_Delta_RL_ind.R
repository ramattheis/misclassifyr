#' Evaluates the log of the prior of Delta at model_to_Delta_RL_ind(psi).
#'
#' @param psi A numeric vector parameterizing `Delta` through `model_to_Delta_RL_ind`.
#' @return A numeric value equal to the log of the flat prior of `Delta` at `psi`, re-scaled for the logit transform.
#' @export
log_prior_Delta_RL_ind = function(psi){

  # J is a deterministic function of psi for model_to_Delta_RL_ind
  J = as.integer((length(psi) + 2)/4)

  # log prior for the row scales for Delta^{(1)} and Delta^{(2)}
  log_prior_Delta_rowscale_1 = logit_link_volume(psi[1:(J-1)])
  log_prior_Delta_rowscale_2 = logit_link_volume(psi[((J-1)+1):(2*(J-1))])

  # log prior for the column scales for Delta^{(1)} and Delta^{(2)}
  psi = psi[-(1:(2*(J-1)))]
  psi1 = psi[1:J]
  psi2 = psi[(J+1):(2*J)]
  col1 = exp(psi1)/(1 + exp(psi1)) # Column margins are probabilities but need not sum to one
  col2 = exp(psi2)/(1 + exp(psi2)) # Column margins are probabilities but need not sum to one
  log_prior_Delta_colscale_1 = sum(log(col1) - log(1-col1))
  log_prior_Delta_colscale_2 = sum(log(col2) - log(1-col2))

  return(log_prior_Delta_rowscale_1 + log_prior_Delta_rowscale_2 +
         log_prior_Delta_colscale_1 + log_prior_Delta_colscale_2 )

}
