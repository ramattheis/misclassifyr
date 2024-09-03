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
phi_0 = NA
psi_0 = NA
split_eta = NA
X_names = list(1:5,1:5)
Y_names =  list(1:5,1:5)
W_names = list(1,2)
X_vals =  list(1:5,1:5)
Y_vals = list(1:5,1:5)
n_mcmc_draws = 2e4
n_burnin = 1e4
thinning_rate = 2
gibbs_poropsal_sd = 0.1
cores = 1
