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

# (Temporary, unbundling to be the arguments of misclassifyr)
tab = misclassifyr_inputs$tab
J = misclassifyr_inputs$J
K = misclassifyr_inputs$K
estimate_beta = T
model_to_Pi = model_to_Pi_NP
model_to_Delta = model_to_Delta_NP_ind
phi_0 = NA
psi_0 = NA
split_eta = NA
X_names = misclassifyr_inputs$X_names
Y_names = misclassifyr_inputs$Y_names
W_names = misclassifyr_inputs$W_names
X_vals = misclassifyr_inputs$X_vals
Y_vals = misclassifyr_inputs$Y_vals
lambda_pos = NA
lambda_dd = NA
optim_maxit = 1e5
optim_tol = 1e-10
optim_stepsize = NA
check_stability = F
stability_sd = 0.1
cores = 6
