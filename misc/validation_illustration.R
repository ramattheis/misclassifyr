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

model_to_Delta = lapply(tab_ABE$tab, function(tb) {

  function(psi) {

    # Transforming phi to return to probabilities (no sum-to-one constraint)
    alpha = exp(psi)/(1+exp(psi))

    # Splitting phi into components associated with Y1 and Y2
    alpha1 = alpha[1:(length(alpha)/2)]
    alpha2 = alpha[(length(alpha)/2 +1):length(alpha)]

    # Computing the marginal distribution of Y1 and Y2
    tabY1 = tb |>
      dplyr::arrange(Y1) |>
      dplyr::group_by(Y1) |>
      dplyr::summarise( n = sum(n)) |>
      as.data.frame()
    FY1 = tabY1$n / sum(tabY1$n)
    tabY2 = tb |>
      dplyr::arrange(Y2) |>
      dplyr::group_by(Y2) |>
      dplyr::summarise( n = sum(n)) |>
      as.data.frame()
    FY2 = tabY2$n / sum(tabY2$n)

    # Computing the misclassification error distribution
    Delta1 = diag(1-alpha1) + FY1 %*% t(alpha1)
    Delta2 = diag(1-alpha2) + FY2 %*% t(alpha2)
    Delta = lapply(1:nrow(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
    Delta = do.call(cbind, Delta)

    return(c(Delta))
  }
})

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

