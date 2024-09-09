# This script illustrates the workflow of estimation with misclassifyr

require(misclassifyr)
require(ggplot2)

# Load the data
synthetic = synthetic_data()

# Estimating the misclassification model, making figures, etc
synthetic_results = misclassifyr(
  tab = synthetic$tab,
  J = list(5,5),
  K = list(5,5),
  estimate_beta = T,
  model_to_Pi = model_to_Pi_NP,
  model_to_Delta = model_to_Delta_NP_ind,
  makeplots = T,
  mle = T,
  bayesian = T,
  log_prior_Delta = log_prior_Delta_NP_ind,
  log_prior_Pi = log_prior_Pi_NP,
  X_names = list(1:5,1:5),
  Y_names =  list(1:5,1:5),
  W_names = list(1,2),
  X_col_name = "Synthetic X",
  Y_col_name = "Synthetic Y",
  X_vals =  list(1:5,1:5),
  Y_vals = list(1:5,1:5),
  cores = 2
)


# Checking results
if(class(synthetic_results$Pi_hat_mle) == "list"){

  # Extracting Pi
  Pi_hat_mle = synthetic_results$Pi_hat_mle[[1]]
  Pi_star = c(synthetic$Pi[[1]])
  Pi_out_df = cbind(Pi_hat_mle, Pi_star) |> as.data.frame()
  ggplot(Pi_out_df, aes(x = Pi_star, y = Pi_hat_mle)) + geom_point() +
    geom_abline(slope=1, intercept = 0, linetype = "dotted") +
    theme_minimal()

  # Extracting Delta
  Delta_hat_mle = synthetic_results$Delta_hat_mle[[1]]
  Delta_star = c(synthetic$Delta[[1]])
  Delta_out_df = cbind(Delta_hat_mle, Delta_star) |> as.data.frame()
  ggplot(Delta_out_df, aes(x = Delta_star, y = Delta_hat_mle)) + geom_point() +
    geom_abline(slope=1, intercept = 0, linetype = "dotted") +
    theme_minimal()

  if(!identical(synthetic_results$trace_plots_Pi[[1]], NA)){
    synthetic_results$trace_plots_Pi[[1]]$`Pr(X = 1, Y = 1)`
    synthetic_results$trace_plots_Pi[[1]]$`Pr(X = 3, Y = 1)`
    synthetic_results$trace_plots_Pi[[1]]$`Pr(X = 4, Y = 2)`
  }

  # How much wider is the (consrvative) set-ID robust CI?
  CCTCI_width = synthetic_results$CCTCI[2] - synthetic_results$CCTCI[1]
  MLE_width = 2*1.96*synthetic_results$se_beta_mle
  Bayes_width = 2*1.96*synthetic_results$posterior_beta_sd
  Bayes_width / MLE_width
  CCTCI_width / MLE_width
}





