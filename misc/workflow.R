# This script illustrates the workflow of estimation with misclassifyr

require(misclassifyr)

# Load the data
data("ancienregime")

# Prepare data for estimation with a continuous outcome and regressor
misclassifyr_inputs_continuous = prep_misclassification_data(
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

# Analysis with regression
out_reg = misclassifyr(
  tab = misclassifyr_inputs_continuous$tab,
  J = misclassifyr_inputs_continuous$J,
  K = misclassifyr_inputs_continuous$K,
  model_to_Delta = model_to_Delta_NP_ind,
  X_names = misclassifyr_inputs_continuous$X_names,
  Y1_names = misclassifyr_inputs_continuous$Y1_names,
  Y2_names = misclassifyr_inputs_continuous$Y2_names,
  X_vals = misclassifyr_inputs_continuous$X_vals,
  Y_vals = misclassifyr_inputs_continuous$Y_vals,
  W_names = misclassifyr_inputs_continuous$W_names,
  mle = T,
  estimate_beta = T,
  X_col_name = "Father's income in 1750",
  Y_col_name = "Son's income in 1780",
  bayesian = T,
  log_prior_Delta = log_prior_Delta_NP_ind,
  cores = 4
)


## Prepare data for estimation with a discrete outcome and regressor
misclassifyr_inputs_categorical = prep_misclassification_data(
  data = ancienregime,
  outcome_1 = "son_occupation_1780",
  outcome_2 = "son_occupation_1770",
  regressor = "father_occupation_1750",
  X_names = c("Vagabond","Metayer","Journalier", "Petit Metiers","Petite Bourgeoisie","Haute Bourgeoisie"),
  Y1_names = c("Vagabond","Metayer","Journalier", "Petit Metiers","Petite Bourgeoisie","Haute Bourgeoisie"),
  Y2_names = c("Vagabond","Metayer","Journalier", "Petit Metiers","Petite Bourgeoisie","Haute Bourgeoisie"),
  controls = "birthplace",
  weights = "linked_weight",
  record_vals = F
)


# Analysis with regression
out_noreg = misclassifyr(
  tab = misclassifyr_inputs_categorical$tab,
  J = misclassifyr_inputs_categorical$J,
  K = misclassifyr_inputs_categorical$K,
  model_to_Delta = model_to_Delta_NP_ind,
  X_names = misclassifyr_inputs_categorical$X_names,
  Y1_names = misclassifyr_inputs_categorical$Y1_names,
  Y2_names = misclassifyr_inputs_categorical$Y2_names,
  W_names = misclassifyr_inputs_categorical$W_names,
  mle = T,
  estimate_beta = F,
  X_col_name = "Father's occupation in 1750",
  Y_col_name = "Son's occupation in 1780",
  bayesian = F,
  log_prior_Delta = log_prior_Delta_NP_ind,
  cores = 4
)
