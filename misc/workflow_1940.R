require(misclassifyr)

linked_ABE = read.csv("~/Dropbox/research/projects/spuriousmobility/data/misc/linked_ABE_s12.csv")

tab_ABE = prep_misclassification_data(
  data = subset(linked_ABE, !is.na(incwage_B)  & !is.na(incwage_B_inst) & !is.na(yrs_educ_A) ),
  outcome_1 = "incwage_B",
  outcome_2 = "incwage_B_inst",
  regressor = "yrs_educ_A",
  outcome_1_bin = "incwage_bin_B",
  outcome_2_bin = "incwage_bin_B_inst",
  regressor_bin = "yrs_educ_bin_A",
  controls = NA,
  weights = NA,
  record_vals = T,
  round_vals = 0
)

out_ABE = misclassifyr(
  tab = tab_ABE$tab,
  J = tab_ABE$J,
  K = tab_ABE$J,
  X_names = tab_ABE$X_names,
  Y1_names = tab_ABE$Y1_names,
  Y2_names = tab_ABE$Y2_names,
  model_to_Delta = model_to_Delta_RL_ind,
  estimate_beta = T,
  X_vals = tab_ABE$X_vals,
  Y_vals = tab_ABE$Y_vals,
  X_col_name = "Years of Education",
  Y_col_name = "Wages"
)

tab_ABE_cat = prep_misclassification_data(
  data = subset(linked_ABE, !is.na(incwage_B)  & !is.na(incwage_B_inst) & !is.na(yrs_educ_A) ),
  outcome_1 = "incwage_bin_B",
  outcome_2 = "incwage_bin_B_inst",
  regressor = "yrs_educ_bin_A",
  X_names = c("0", "4", "7", "10", "12", "14", "16", "18"),
  Y1_names = c("= $0", "< $500", "< $1000", "< $1500", "< $2000", "< $2500", "< $3000", "> $3000"),
  Y2_names = c("= $0", "< $500", "< $1000", "< $1500", "< $2000", "< $2500", "< $3000", "> $3000"),
  controls = NA,
  weights = NA,
  record_vals = F,
  round_vals = 0
)

out_ABE_cat = misclassifyr(
  tab = tab_ABE_cat$tab,
  J = tab_ABE_cat$J,
  K = tab_ABE_cat$J,
  X_names = tab_ABE_cat$X_names,
  Y1_names = tab_ABE_cat$Y1_names,
  Y2_names = tab_ABE_cat$Y2_names,
  model_to_Delta = model_to_Delta_NP_ind,
  estimate_beta = F,
  X_col_name = "Years of Education",
  Y_col_name = "Wages"
)
