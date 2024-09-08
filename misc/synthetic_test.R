# This script illustrates the workflow of estimation with misclassifyr

require(misclassifyr)

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

