# This script illustrates the workflow of estimation with misclassifyr

require(misclassifyr)

# Load the data
synthetic = synthetic_data()

# (Temporary, unbundling to be the arguments of misclassifyr)
tab = synthetic$tab
J = list(5,5)
K = list(5,5)
estimate_beta = T
model_to_Pi = model_to_Pi_NP
model_to_Delta = model_to_Delta_NP_ind
makeplots = T
phi_0 = NA
psi_0 = NA
X_names = list(1:5,1:5)
Y_names =  list(1:5,1:5)
W_names = list(1,2)
X_col_name = "Synthetic X"
Y_col_name = "Synthetic Y"
X_vals =  list(1:5,1:5)
Y_vals = list(1:5,1:5)
mle = T
optim_tol = 1e-8
optim_maxit = 1e5
check_stability = F
stability_sd = 0.1
bayesian = T
log_prior_Delta = log_prior_Delta_NP_ind
log_prior_Pi = log_prior_Pi_NP
n_mcmc_draws = 2e3
n_burnin = 1e3
thinning_rate = 1
gibbs_proposal_sd = 0.02
cores = 2
