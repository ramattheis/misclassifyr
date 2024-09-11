# This script illustrates the workflow of estimation with misclassifyr

require(misclassifyr)

# Load the data
data("ancienregime")

# Prepare data for estimation
misclassifyr_inputs = prep_misclassification_data(
  ancienregime,
  "son_income_1780",
  "son_income_1770",
  "father_income_1750",
  c("birthplace","birthyear"),
  "linked_weight",
  record_vals = T)

# Analysis
out = misclassifyr(
  tab = misclassifyr_inputs$tab,
  J = misclassifyr_inputs$J,
  K = misclassifyr_inputs$K,
  model_to_Delta = model_to_Delta_NP_ind,
  X_names = misclassifyr_inputs$X_names,
  Y_names = misclassifyr_inputs$Y_names,
  X_vals = misclassifyr_inputs$X_vals,
  Y_vals = misclassifyr_inputs$Y_vals,
  W_names = misclassifyr_inputs$W_names,
  estimate_beta = T,
  X_col_name = "father_income_1750",
  Y_col_name = "son_income_1780",
  bayesian = T,
  log_prior_Delta = log_prior_Delta_NP,
  cores = 4
)



# (Temporary, unbundling to be the arguments of misclassifyr)
tab = misclassifyr_inputs$tab
J = misclassifyr_inputs$J
K = misclassifyr_inputs$K
estimate_beta = T
makeplots = T
model_to_Pi = model_to_Pi_NP
model_to_Delta = model_to_Delta_NP
phi_0 = NA
psi_0 = NA
X_names = misclassifyr_inputs$X_names
Y_names = misclassifyr_inputs$Y_names
W_names = misclassifyr_inputs$W_names
X_vals = misclassifyr_inputs$X_vals
Y_vals = misclassifyr_inputs$Y_vals
X_col_name = "Father's Income 1750"
Y_col_name = "Son's Income 1780"
mle = T
optim_tol = 1e-8
optim_maxit = 1e5
check_stability = F
stability_sd = 0.1
bayesian = T
log_prior_Delta = log_prior_Delta_NP
log_prior_Pi = log_prior_Pi_NP
n_mcmc_draws = 1e4
n_burnin = 5e3
thinning_rate = 1
gibbs_proposal_sd = 0.02
cores = 2
