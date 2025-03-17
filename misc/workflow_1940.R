require(misclassifyr)

linked_ABE = read.csv("~/Dropbox/research/projects/spuriousmobility/data/temp/linked_075.csv")


# Tabulating
linked_tab = prep_misclassification_data(
  data = subset(linked_ABE, !is.na(incwage_B) &
                  !is.na(incwage_B_inst) &
                  !is.na(yrs_educ_A)), # Removing NAs
  outcome_1 = "incwage_B",
  outcome_2 = "incwage_B_inst",
  regressor = "yrs_educ_A",
  outcome_1_bin = "incwage_bin_B",
  outcome_2_bin = "incwage_bin_B_inst",
  regressor_bin = "yrs_educ_bin_A",
  record_vals = T,
  round_vals = 0
)

# Posterior of the nonparametric model
bayes_NP_ind = misclassifyr(
  tab = linked_tab$tab,
  J = linked_tab$J,
  K = linked_tab$J,
  X_names = linked_tab$X_names,
  Y1_names = linked_tab$Y1_names,
  Y2_names = linked_tab$Y2_names,
  model_to_Delta = model_to_Delta_NP_ind,
  mle = F,
  bayesian = T,
  n_mcmc_draws = 1e3,
  n_burnin = 2e2
)
