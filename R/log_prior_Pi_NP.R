#' Evaluates the log of the prior of Pi at model_to_Pi_NP(phi).
#'
#' @param phi A numeric vector parameterizing `Pi` through `model_to_Pi`.
#' @return A numeric value equal to the log of the flat prior of `Pi` at `phi`, re-scaled for the logit transform.
#' @export
log_prior_Pi_NP = function(phi){

  # For nonparametric Pi, the volume distortion is just logit_link_volume
  log_det_logit_Pi = logit_link_volume(phi)

  return(log_det_logit_Pi)

}
