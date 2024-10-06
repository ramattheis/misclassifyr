require(misclassifyr)

linked = read.csv("~/Dropbox/research/projects/spuriousmobility/data/temp/validation_example.csv")

tab_ABE = prep_misclassification_data(
  data =  subset(linked, !is.na(incwage_B) &
                   !is.na(incwage_B_inst) &
                   !is.na(yrs_educ_A)), # Removing NAs ,
  outcome_1 = "incwage_B",
  outcome_2 = "incwage_B_inst",
  regressor = "yrs_educ_A",
  outcome_1_bin = "incwage_bin_B",
  outcome_2_bin = "incwage_bin_B_inst",
  regressor_bin = "yrs_educ_bin_A",
  controls = c("bpl","cohort"),
  weights = NA,
  record_vals = T,
  round_vals = 0
)

tab = tab_ABE$tab
J = tab_ABE$J
K = tab_ABE$J
X_names = tab_ABE$X_names
Y1_names = tab_ABE$Y1_names
Y2_names = tab_ABE$Y2_names
W_names = tab_ABE$W_names
model_to_Pi = model_to_Pi_NP
model_to_Delta = model_to_Delta_NP
estimate_beta = T
phi_0 = NA
psi_0 = rep(-2,16)
X_vals = tab_ABE$X_vals
Y_vals = tab_ABE$Y_vals
X_col_name = "Years of Education"
Y_col_name = "Wages"
makeplots = T
misclassification_size = 0.2
estimate_beta = T
mle = T
optim_tol = 1e-8
optim_maxit = 1e5
check_stability = F
stability_sd = 0.1
bayesian = T
log_prior_Pi = log_prior_Pi_NP
log_prior_Delta = log_prior_Delta_NP_ind
n_mcmc_draws = 2e3
n_burnin = 1e3
thinning_rate = 1
gibbs_proposal_sd = 0.1
cores = 1

