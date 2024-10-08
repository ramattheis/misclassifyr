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


Delta_info = make_empirical_Delta_RL_common_alpha_mixed_NP(tab_ABE$tab, tab_ABE$J)

mcout_emp = misclassifyr(
  tab = tab_ABE$tab,
  J = tab_ABE$J,
  K = tab_ABE$J,
  X_names = tab_ABE$X_names,
  Y1_names = tab_ABE$Y1_names,
  Y2_names = tab_ABE$Y2_names,
  W_names = tab_ABE$W_names,
  model_to_Pi = model_to_Pi_NP,
  model_to_Delta = Delta_info$model_to_Delta,
  phi_0 = NA,
  psi_0 = Delta_info$psi_0,
  X_col_name = "Years of Education",
  Y_col_name = "Wages",
  cores = 3
)

beta_mle = Pi_to_beta(
  X_vals = tab_ABE$X_vals,
  Y_vals = tab_ABE$Y_vals,
  W_weights = mcout_emp$W_weights,
  W_names = tab_ABE$W_names,
  mle = T,
  Pi_mle = mcout_emp$Pi_hat_mle,
  cov_Pi = mcout_emp$cov_Pi_mle
)

