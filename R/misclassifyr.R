#' misclassifyr
#'
#' @param data Tabulated data or a list of tabulated data split by controls
#' @param model_options Options for the model to be estimated
#' @param estimate_options Options for downstream estimates
#' @return An object that includes estimates and information from the estimation process
#' @export
misclassifyr <- function(data, model_options, estimate_options) {
  # Call ll, delta, and pi as needed
  log_likelihood <- ll(params, data)
  delta_params <- delta(model_params)
  pi_params <- pi(other_params)

  # Additional estimation and inference
  estimates <- estimate_model(log_likelihood, delta_params, pi_params, data)

  # Return results
  return(list(
    estimates = estimates,
    log_likelihood = log_likelihood,
    delta_params = delta_params,
    pi_params = pi_params
  ))
}
