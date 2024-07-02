#' misclassifyr
#'
#' This package provides a menu of options for estimation and inference of misclassification models in which the analyst has access to two noisy measures, `Y1` and `Y2` of a latent outcome `Y*`, a correctly measured covariate `X`, and discrete controls `W`.
#'
#' @param tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should be numeric with names `Y1`, `Y2`, `X`, and `n` where `Y1` and `Y2` take each value between `1` and `J`, `X` takes each value between `1` and `K`, and
#' @param J An integer corresponding to the number of unique values of `Y1` and `Y2`.
#' @param K An integer less than or equal to `J` corresponding to the number of unique values of `X`.
#' @param model_to_Pi A function mapping the parameters of a model for the joint distribution to the joint distribution \Pi
#' @param model_to_Delta A function mapping the parameters of a model to the conditional distribution Y1, Y2 | Y*, \Delta
#' @param phi_0 A numeric vector providing the starting location for optimization for the argument to model_to_Pi.
#' @param psi_0 A numeric vector providing the starting location for optimization for the argument to model_to_Delta.
#' @param split_eta An integer indicating where to split the vector `eta` in `phi` and `psi`, the arguments to `model_to_Pi` and `model_to_Delta` respectively.
#' @param lambda_pos scales the penalty for violations of positivity (i.e. all probabilities should be positive).
#' @param lambda_dd scales the penalty for violations of diagonal dominance.
#' @param optim_maxit An integer for the maximum number of iterations in numerical optimization, passed to `optim()`
#' @param optim_tol A positive number defining convergence in numerical optimization, passed to `optim()`
#' @param optim_stepsize A positive number for the step size in the numerical gradient, passed to `optim()`
#' @param check_stability A logical value indicating whether to perform a stability test for the numerical optimizer.
#' @param estimation_options Options for downstream estimates.
#' @param formula A regression formula
#' @return An object that includes estimates and information from the estimation process
#' @export
misclassifyr <- function(
    tab,
    J,
    K,
    model_to_Pi,
    model_to_Delta,
    phi_0 = NA,
    psi_0 = NA,
    split_eta = NA,
    lambda_pos = NA,
    lambda_dd = NA,
    optim_maxit = 1e3,
    optim_tol = 1e-9,
    optim_stepsize = NA,
    check_stability = F,
    stability_sd = 0.1,
    estimation_options = NA,
    formula = NA) {

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  # Note, some input errors depend on future quantities and appear later in the code.

  # Throwing an error for argument type
  if(!is.data.frame(tab)){
    stop("Error: `tab` should be a data frame.")}
  if(!is.function(model_to_Pi)){
    stop("Error: `model_to_Pi` should be a function.")}
  if(!is.function(model_to_Delta)){
    stop("Error: `model_to_Delta` should be a function.")}
  if(!is.na(split_eta) & !is.integer(split_eta)){
    stop("Error: If provided, `split_eta` should be an integer.")
  }

  # Throwing errors if `tab` isn't properly formatted
  if(!identical(setdiff(c("Y1","Y2","X","n"), colnames(tab)),character(0))){
    stop("Error: `tab` should have four columns: `Y1`, `Y2`, `X`, and `n`.")
  }
  if( typeof(tab$Y1) != "integer" | max(tab$Y1) > J | min(tab$Y1) < 1 |
      typeof(tab$Y2) != "integer" | max(tab$Y2) > J | min(tab$Y2) < 1 |
      typeof(tab$X)  != "integer" | max(tab$X)  > K | min(tab$X) < 1 ){
    stop("Error: `Y1`, `Y2`, and `X` should take integer values between `1` and `J` or `K`.")
  }
  if( !(typeof(tab$n) %in% c("double","integer") ) | min(tab$n) < 0){
    stop("Error: `n` should take non-negative integer values representing counts of unique values of `X`, `Y1`, and `Y2`")
  }

  # Throwing an error if J is not weakly larger than K
  if( K > J ){
    stop("Error: `K` should be weakly less than `J`.")
  }

  # Throwing an error if no starting location is found for user-defined model_to_Pi or model_to_Delta
  if(is.na(phi_0) & identical(attr(model_to_Pi,"name"),NULL)){
    stop("Error: `phi_0` must be defined if `model_to_Pi` is user-defined.")
  }
  if(is.na(psi_0) & identical(attr(model_to_Delta,"name"),NULL)){
    stop("Error: `psi_0` must be defined if `model_to_Delta` is user-defined.")
  }

  #------------------------------------------------------------
  # Ensuring tab is "balanced"
  #------------------------------------------------------------

  # Making a "balanced" data frame of the correct dimension with n=0
  emptydf = expand.grid(1:J,1:J,1:K)
  colnames(emptydf) = c("X","Y1","Y2")
  emptydf$n = 0

  # Making sure tab is a data.frame and not a data.table or tibble
  tab = as.data.frame(tab)

  # binding and collapsing the empty df with tab
  tab = tab[,c("X","Y1","Y2","n")]
  tab = rbind(tab,emptydf)
  tab = dplyr::group_by(tab,X,Y1,Y2) |>
    dplyr::summarise(n = sum(n)) |>
    as.data.frame()

  # Ordering tab data by X, Y_1, Y_2
  tab = tab[order(tab$Y2, tab$Y1, tab$X, decreasing = F),]

  #------------------------------------------------------------
  # Defining the objective function for estimation
  #------------------------------------------------------------

  # Defining split_eta for default functions
  if(identical("model_to_Pi_NP",
               attr(model_to_Pi,"name"))){
    split_eta = J*K-1
  }

  # Choosing penalty scales to be consistent with the size of the sample
  ntotal = sum(tab$n)
  # lambda_pos scales the penalty for violations of positivity (i.e. all probabilities should be positive)
  if(identical(lambda_pos, NA)){
    lambda_pos = ntotal^2
  }
  # lambda_dd scales the penalty for violations of diagonal dominance, as defined in Mattheis (2024).
  if(identical(lambda_dd, NA)){
    lambda_dd = ntotal^2
  }

  # Defining the objective function
  objective = function(eta){

    # Mapping model arguments to the entries of the joint distribution
    Pi = model_to_Pi(eta[1:split_eta])
    Delta = model_to_Delta(eta[(split_eta+1):length(eta)] )

    # Computing the log likelihood of the data given Pi and Delta
    ll =  loglikelihood(c(c(Pi), c(Delta)))

    # Returning -1x the log likelihood + penalties
    return(-1*ll)

  }

  #------------------------------------------------------------
  # Setting the starting location for optimization
  #------------------------------------------------------------

  if(is.na(phi_0)){
    if(identical(attr(model_to_Pi,"name"),"model_to_Pi_NP")){
      # Default starting location for phi_0 is uniform Pi
      phi_0 = rep(log(1/J^2),J^2-1)
    }
  }

  if(is.na(psi_0)){
    if(identical(attr(model_to_Delta,"name"),"model_to_Delta_RL_ind")){
      # Default starting location for psi_0 is uniform margin and 4/5 diagonal
      row_0 = rep(1/J,J)
      col_0 = rep(1/5,J)
      psi_0 = log(c(row_0[-J],row_0[-J],col_0,col_0))
      rm(list = c("row_0","col_0"))
    }
    if(identical(attr(model_to_Delta,"name"),"model_to_Delta_NP_ind")){
      # Default starting location for psi_0 is RL with uniform margin and 4/5 diagonal
      row_0 = rep(1/J,J)
      col_0 = rep(1/5,J)
      Delta_0 = diag(J)*(1-col_0) + outer(row_0,col_0,"*")
      psi_0 = log(c(c(Delta_0[-J,]),c(Delta_0[-J,])))
      rm(list = c("row_0","col_0","Delta_0"))
    }
    if(identical(attr(model_to_Delta,"name"),"model_to_Delta_NP")){
      # Default starting location for psi_0 is RL with uniform margin and 4/5 diagonal
      row_0 = rep(1/J,J)
      col_0 = rep(1/5,J)
      Delta_0 = diag(J)*(1-col_0) + outer(row_0,col_0,"*")
      # Building the joint distribution Y1, Y2 conditional on Y*
      Delta_0 = lapply(1:ncol(Delta_0), function(j) diag(Delta_0[j,]) %*% t(Delta_0))
      Delta_0 = do.call(cbind, Delta_0)
      psi_0 = log(c(Delta_0[,-J^2]))
      rm(list = c("row_0","col_0","Delta_0"))
    }

  }

  # Defining the initial eta
  eta_0 = c(phi_0,psi_0)

  # Additionally, setting the default step size for the numerical derivative
  if(identical(optim_stepsize,NA)){
    optim_stepsize = rep(1e-6, length(eta_0))
  }

  # Testing whether optim_stepsize has the right length
  if(!identical(length(eta_0),length(optim_stepsize))){
    stop("Error: eta_0 and optim_stepsize should be the same length.")
  }

  #------------------------------------------------------------
  # Passing common objects to the shared environment
  #------------------------------------------------------------

  # Retrieving the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))

  # Passing common objects to the shared environment
  misclassifyr_env$tab = tab
  misclassifyr_env$J = J
  misclassifyr_env$K = K
  misclassifyr_env$lambda_pos = lambda_pos
  misclassifyr_env$lambda_dd = lambda_dd

  #------------------------------------------------------------
  # Optimizing objective w.r.t. eta
  #------------------------------------------------------------

  # Optimizing objective w.r.t. eta
  out = optim(par = eta_0,
              fn = objective,
              method = "BFGS",
              hessian = T,
              control = list(maxit = optim_maxit,
                             reltol = optim_tol,
                             abstol = optim_tol,
                             ndeps = optim_stepsize))


  #------------------------------------------------------------
  # Testing the stability of the optimizer
  #------------------------------------------------------------

  if((check_stability)){

    # Recording initial eta hat
    eta_hat1 = out$par

    # Initializing total inconsistency
    inconsistency = 0

    for(draw in 1:10){

      # Additional estimates on perturbed starting locations
      out2 = optim(par = eta_0 + rnorm(length(eta_0),sd = stability_sd),
                   fn = objective,
                   method = "BFGS",
                   control = list(maxit = optim_maxit,
                                  reltol = optim_tol,
                                  abstol = optim_tol,
                                  ndeps = optim_stepsize))
      eta_hat2 = out2$par

      # Adding the (absolute) difference between e^eta_hat1 and e^eta_hat2
      inconsistency = inconsistency + sum(abs(exp(eta_hat1) - exp(eta_hat2)))
    }

    if(inconsistency < 0.01){
      numerical_stability = paste0("Estimates stable in a ", stability_sd, " SD normal ball around eta_0.")
    } else {
      numerical_stability = paste0("Estimates stable in a ", stability_sd, " SD normal ball around eta_0, the sum of absolute differences is ", round(inconsistency,3) ," over 10 iterations.")
    }

  } else {
    numerical_stability = "Numerical stability not tested."
  }



  #------------------------------------------------------------
  # Returning estimates and other info
  #------------------------------------------------------------

  # Extracting results
  eta_hat = out$par

  # Return results
  return(list(
    eta_hat = eta_hat,
    log_likelihood = -1*out$value,
    convergence = out$convergence,
    numerical_stability = numerical_stability
  ))
}
