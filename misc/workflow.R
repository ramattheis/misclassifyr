# This script illustrates the workflow of estimation with misclassifyr

require(misclassifyr)

# Load the data
data("ancienregime")

# arguments for prep_misclassification
data = ancienregime
outcome_1 = "son_income_1780"
outcome_2 = "son_income_1770"
regressor = "father_income_1750"
outcome_1_bin = "son_occupation_1780"
outcome_2_bin = "son_occupation_1770"
regressor_bin = "father_occupation_1750"
controls = c("birthplace","birthyear")
record_vals = T
weights = "linked_weight"
X_names = NA
Y1_names = NA
Y2_names = NA



# Prepare data for estimation
misclassifyr_inputs = prep_misclassification_data(
  data = ancienregime,
  outcome_1 = "son_income_1780",
  outcome_2 = "son_income_1770",
  regressor = "father_income_1750",
  outcome_1_bin = "son_occupation_1780",
  outcome_2_bin = "son_occupation_1770",
  regressor_bin = "father_occupation_1750",
  controls = "birthplace",
  weights = "linked_weight",
  record_vals = T,
  round_vals = 0
)

# Analysis
out = misclassifyr(
  tab = misclassifyr_inputs$tab,
  J = misclassifyr_inputs$J,
  K = misclassifyr_inputs$K,
  model_to_Delta = model_to_Delta_NP_ind,
  X_names = misclassifyr_inputs$X_names,
  Y1_names = misclassifyr_inputs$Y1_names,
  Y2_names = misclassifyr_inputs$Y2_names,
  X_vals = misclassifyr_inputs$X_vals,
  Y1_vals = misclassifyr_inputs$Y1_vals,
  Y2_vals = misclassifyr_inputs$Y2_vals,
  W_names = misclassifyr_inputs$W_names,
  mle = F,
  estimate_beta = T,
  X_col_name = "father_outcome_1750",
  Y_col_name = "son_outcome_1780",
  bayesian = T,
  log_prior_Delta = log_prior_Delta_NP_ind,
  cores = 1
)



# misclassifyr arguments for debugging
tab = misclassifyr_inputs$tab
J = misclassifyr_inputs$J
K = misclassifyr_inputs$K
model_to_Pi = model_to_Pi_NP
model_to_Delta = model_to_Delta_NP_ind
phi_0 = NA
psi_0 = NA
makeplots = T
X_names = misclassifyr_inputs$X_names
Y1_names = misclassifyr_inputs$Y1_names
Y2_names = misclassifyr_inputs$Y2_names
X_vals = misclassifyr_inputs$X_vals
Y_vals = misclassifyr_inputs$Y_vals
W_names = misclassifyr_inputs$W_names
mle = T
check_stability = F
stability_sd = 0.1
optim_maxit = 1e5
optim_tol = 1e-8
estimate_beta = F
X_col_name = "father_income_1750"
Y_col_name = "son_income_1780"
bayesian = F
log_prior_Pi = log_prior_Pi_NP
log_prior_Delta = log_prior_Delta_NP_ind
gibbs_proposal_sd = 0.01
n_mcmc_draws = 1e3
n_burnin = 5e2
thinning_rate = 1
cores = 6
