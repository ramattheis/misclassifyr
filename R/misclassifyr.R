#' misclassifyr
#'
#' This function provides a menu of options for estimation and inference of misclassification models in which the analyst has access to two noisy measures, `Y1` and `Y2` of a latent outcome `Y*`, a correctly measured covariate `X`, and discrete controls `W`.
#'
#' @import dplyr
#' @import parallel
#' @import numDeriv
#' @import ggplot2
#'
#' @param tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should be numeric with names `Y1`, `Y2`, `X`, and `n` where `Y1` and `Y2` take each value between `1` and `J`, `X` takes each value between `1` and `K`, and
#' @param J An integer or list corresponding to the number of unique values of `Y1` and `Y2`.
#' @param K An integer or list corresponding to the number of unique values of `X`.
#' @param model_to_Pi A function or list of functions mapping the parameters of a model for the joint distribution to the joint distribution `\eqn{\Pi}`.
#' @param model_to_Delta A function or list of functions mapping the parameters of a model to the conditional distribution Y1, Y2 | Y*, `\eqn{\Delta}`.
#' @param makeplots A logical value for whether to make trace plots and plots of \eqn{\Pi} and \eqn{Delta}. Defaults to TRUE.
#' @param phi_0 A numeric vector or list of numeric vectors providing the starting location for optimization for the argument to model_to_Pi.
#' @param psi_0 A numeric vector or list of numeric vectors providing the starting location for optimization for the argument to model_to_Delta.
#' @param split_eta An integer or list indicating where to split the vector `eta` in `phi` and `psi`, the arguments to `model_to_Pi` and `model_to_Delta` respectively.
#' @param X_names A character vector or list corresponding to the values of the regressor X.
#' @param Y_names A character vector or list corresponding to the values of the outcome Y.
#' @param X_col_name A character vector corresponding to the variable of the regressor X, used only for plots.
#' @param Y_col_name A character vector corresponding to the variable of the outcome Y, used only for plots.
#' @param W_names A character vector corresponding to the values of the control W in each cell.
#' @param estimate_beta A logical value indicating whether to regress Y on X.
#' @param X_vals A numeric vector or list of numeric vectors providing the values of X associated with the columns of Pi.
#' @param Y_vals A numeric vector or list of numeric vectors providing the values of Y associated with the rows of Pi.
#' @param mle A logical value indicating whether to estimate Pi and Delta via MLE. Defaults to TRUE.
#' @param optim_tol A numeric value giving the relative tolerance for optimization with the optim.
#' @param optim_maxit An integer giving the maximum number of iterations for optim.
#' @param check_stability A logical value indicating whether to perform a more rigorous stability test for the numerical optimizer.
#' @param stability_sd A numerical value giving the standard deviation of the noise added to the initial parameter value for the stability test of the MLE.
#' @param bayesian A logical value indicating whether or not to compute the posterior of values.
#' @param log_prior_Pi A function or list of functions evaluating the log of the prior of Pi at `phi` (in logs!).
#' @param log_prior_Delta A function or list of functions evaluating the log of the prior of Delta at `psi` (in logs!).
#' @param n_mcmc_draws An integer corresponding to the length of the MCMC chain.
#' @param n_burnin An integer giving the length of the burn-in period for each MCMC chain, must be shorter than `n_mcmc_draws`.
#' @param thinning_rate An integer indicating how frequently to record posterior draws from the MCMC chain -- e.g. a `thinning_rate` of 2 records every other draw.
#' @param gibbs_proposal_sd A numeric value giving the standard deviation for the proposal distribution in each Gibbs step.
#' @param cores An integer for the number of CPUs available for parallel processing.
#' @return A list containing the following components:
#'   - `$Pi_hat_MLE`: The MLE estimate of the joint distribution of \eqn{X} and \eqn{Y^*}, \eqn{\Pi}.
#' @export
misclassifyr <- function(
    tab,
    J,
    K,
    model_to_Pi = model_to_Pi_NP,
    model_to_Delta = model_to_Delta_NP_ind,
    makeplots = T,
    phi_0 = NA,
    psi_0 = NA,
    X_names = NA,
    Y_names = NA,
    W_names = NA,
    estimate_beta = F,
    estimate_betas = F,
    X_vals = NA,
    Y_vals = NA,
    X_col_name = "X",
    Y_col_name = "Y",
    mle = T,
    optim_tol = 1e-8,
    optim_maxit = 1e5,
    check_stability = F,
    stability_sd = 0.1,
    bayesian = F,
    log_prior_Pi = log_prior_Pi_NP,
    log_prior_Delta = log_prior_Delta_NP_ind,
    n_mcmc_draws = 1e4,
    n_burnin = 5e3,
    thinning_rate = 1,
    gibbs_proposal_sd = 0.1,
    cores = 1){

  #------------------------------------------------------------
  # Catching errors in some variables
  # other input errors caught in MisclassMLE()
  #------------------------------------------------------------

  if(class(bayesian) != "logical"){stop("Error: `bayesian` should be TRUE or FALSE.")}
  if(class(mle) != "logical"){stop("Error: `mle` should be TRUE or FALSE.")}
  if(!bayesian&!mle){stop("Error: either `bayesian` or `mle` should be TRUE.")}

  if(!identical(cores%%1,0)){
    stop("Error: `cores` should be an integer.")
  } else if(parallel::detectCores() < cores){
    stop("Error: You requested more cores than appear available on your machine.")
  }

  if(!identical(optim_maxit%%1,0)){ stop("Error: `optim_maxit` should be an integer.") }
  if(identical(as.numeric(gibbs_proposal_sd),NA)){ stop("Error: `gibbs_proposal_sd` should be numeric.") }

  if(!identical(n_mcmc_draws%%1,0)){ stop("Error: `n_mcmc_draws` should be an integer.") }
  if(!identical(n_burnin%%1,0)){ stop("Error: `n_burnin` should be an integer.") }
  if(!identical(thinning_rate%%1,0)){ stop("Error: `thinning_rate` should be an integer.") }
  if(identical(as.numeric(gibbs_proposal_sd),NA)){ stop("Error: `gibbs_proposal_sd` should be numeric.") }

  # Throwing an error if n_burnin or n_mcmc_draws is too small
  if(n_mcmc_draws < 10000){ stop("Error: `n_mcmc_draws` is too small. Choose a value of at least 10000.")}
  if(n_burnin < 5000){ stop("Error: `n_burnin` is too small. Choose a value of at least 5000.")}

  # Throwing an error if the relative size of n_mcmc_draws/ n_burnin /
  if(n_mcmc_draws - n_burnin < 1000){ stop("Error: `n_burnin` is too close to `n_mcmc_draws`. Choose a smaller `n_burnin` or a larger `n_mcmc_draws`.")}
  if((n_mcmc_draws - n_burnin)/thinning_rate < 1000){ stop("Error: `thinning_rate` is too large. Choose a smaller `n_thinning_rate` or a larger gap between `n_burnin` and `n_mcmc_draws`.")}

  #-----------------------------
  # Recording the types of each object
  #-----------------------------
  input_types = c(class(tab),
                  class(J),
                  class(K),
                  class(model_to_Pi),
                  class(model_to_Delta),
                  class(phi_0),
                  class(psi_0),
                  class(X_names),
                  class(Y_names),
                  class(X_vals),
                  class(Y_vals)) |>
    as.data.frame()

  colnames(input_types) = "type"

  input_types$object = c(
    "tab",
    "J",
    "K",
    "model_to_Pi",
    "model_to_Delta",
    "phi_0",
    "psi_0",
    "X_names",
    "Y_names",
    "X_vals",
    "Y_vals")

  #-----------------------------
  # Are any of the objects lists? If so, making all objects lists
  #-----------------------------

  if(any(input_types$type == "list")){
    # Throwing an error if tab is not a list (nothing should be a list if tab is not a list)
    if(subset(input_types, object == "tab")$type != "list"){stop("Error: No arguments should be lists if `tab` is not a list.")}
    if(!all(input_types$type == "list")){
      # Returning a message that some objects will be copied.
      message(paste0("The following objects were not provided as a list and will be copied across covariate cells: ",
                     paste0(subset(input_types, type != "list")$object, collapse = ", ")))
      # Determining the length of each provided list
      list_lengths = sapply(subset(input_types, type == "list")$object,
                            function(list_name) length(get(list_name)))
      # Throwing an error if not all lists are the same length
      if(length(unique(list_lengths)) > 1){stop("Error: Provided lists are not all the same length.")}
      list_length = unique(list_lengths)
      # Making a list of copies for all objects that are not already lists
      invisible(
        lapply(subset(input_types, type != "list")$object,
                 function(obj_name) {
                      obj = get(obj_name) # Grabbing the object from the environment
                      obj_list = replicate(obj, n = list_length, simplify = F) # Making copies
                      assign(obj_name, obj_list, envir = parent.frame(2)) # Assigning back to the misclassifyr()-level environment
                  }
               )
      )
    }
  }


  #-----------------------------
  # Converting inputs to a list (or list of lists)  for estimate_misclassification
  #-----------------------------

  if(class(tab) == "list"){
    # rebundling each list as a list of lists
    misclassification_inputs = lapply(seq_along(tab), function(j) list(
      tab = tab[[j]],
      J = J[[j]],
      K = K[[j]],
      model_to_Pi = model_to_Pi[[j]],
      model_to_Delta = model_to_Delta[[j]],
      phi_0 = phi_0[[j]],
      psi_0 = psi_0[[j]],
      X_names = X_names[[j]],
      Y_names = Y_names[[j]],
      W_names = W_names[j],
      X_vals = X_vals[[j]],
      Y_vals = Y_vals[[j]]
      ))

  } else {
    # Bundling inputs into a list
    misclassification_inputs = list(
      tab = tab,
      J = J,
      K = K,
      model_to_Pi = model_to_Pi,
      model_to_Delta = model_to_Delta,
      phi_0 = phi_0,
      psi_0 = psi_0,
      X_names = X_names,
      Y_names = Y_names,
      W_names = W_names,
      X_vals = X_vals,
      Y_vals = Y_vals
    )
  }

  #------------------------------------------------------------
  # Defining a function estimate_misclassification to estimate Pi and Delta
  #------------------------------------------------------------

  estimate_misclassification = function(misclassification_input){

    #------------------------------------------------------------
    # Unpacking inputs
    #------------------------------------------------------------

    tab = misclassification_input$tab
    J = misclassification_input$J
    K = misclassification_input$K
    model_to_Pi = misclassification_input$model_to_Pi
    model_to_Delta = misclassification_input$model_to_Delta
    phi_0 = misclassification_input$phi_0
    psi_0 = misclassification_input$psi_0
    X_names = misclassification_input$X_names
    Y_names = misclassification_input$Y_names
    W_names = misclassification_input$W_names
    X_vals = misclassification_input$X_vals
    Y_vals = misclassification_input$Y_vals

    #------------------------------------------------------------
    # Catching input errors
    #------------------------------------------------------------

    # Throwing an error for argument type
    if(!is.data.frame(tab)){
      stop("Error: `tab` should be a data frame.")}
    if(!is.function(model_to_Pi)){
      stop("Error: `model_to_Pi` should be a function.")}
    if(!is.function(model_to_Delta)){
      stop("Error: `model_to_Delta` should be a function.")}

    # Throwing errors if `tab` isn't properly formatted
    if(!identical(setdiff(c("Y1","Y2","X","n"), colnames(tab)),character(0))){
      stop("Error: `tab` should have four columns: `Y1`, `Y2`, `X`, and `n`.")
    }
    if( !(class(tab$n) %in% c("double","numeric","integer") ) | min(tab$n) < 0){
      stop("Error: `n` should take non-negative integer values representing counts (or weighted counts) of unique values of `X`, `Y1`, and `Y2`")
    }
    # Throwing an error if no starting location is found for user-defined model_to_Pi or model_to_Delta
    if(is.na(phi_0) & identical(attr(model_to_Pi,"name"),NULL)){
      stop("Error: `phi_0` must be defined if `model_to_Pi` is user-defined.")
    }
    if(is.na(psi_0) & identical(attr(model_to_Delta,"name"),NULL)){
      stop("Error: `psi_0` must be defined if `model_to_Delta` is user-defined.")
    }
    # Throwing an error if tab is not balanced
    if( any(duplicated(tab[,c("X","Y1","Y2")])) | nrow(tab) != J^2*K){
      stop("Error: `tab` appears not to be balanced across X, Y1, and Y2.")
    }
    # Throwing an error if model_to_Pi doesn't accept additional arguments
    if(!("..." %in% names(formals(model_to_Pi)))){
      stop("Error: `model_to_Pi` must accept additional arguments `...`.")
    }


    #------------------------------------------------------------
    # Setting the starting location for optimization and/or MCMC
    #------------------------------------------------------------

    if(identical(phi_0,NA)){
      if(identical(attr(model_to_Pi,"name"),"model_to_Pi_NP")){
        # Default starting location for phi_0 is the empirical distribution of X and Y1
        tab_xy = tab |>
          dplyr::group_by(X,Y1) |>
          dplyr::summarise(n = sum(n),.groups = "drop") |>
          as.data.frame()
        tab_xy = tab_xy[order(tab_xy$Y1, tab_xy$X),]
        tab_xy$p = tab_xy$n/sum(tab_xy$n)
        phi_0 = softlog(tab_xy$p[1:(J*K-1)] /tab_xy$p[J*K])

        # Default starting location for phi_0 is flat
        # phi_0 = softlog(rep(1/(J*K), J*K-1 )/(1/(J*K)) )
      }
    }

    if(identical(psi_0,NA)){
      if(identical(attr(model_to_Delta,"name"),"model_to_Delta_RL_ind")){
        # Default starting location for psi_0 is uniform margin and 4/5 diagonal
        row_0 = rep(1/J,J-1)
        col_0 = rep(1/5,J)
        row_1 = softlog(row_0 / (1/J))
        row_2 = softlog(row_0 / (1/J))
        col_1 = softlog(col_0 / (1-col_0))
        col_2 = softlog(col_0 / (1-col_0))
        psi_0 = c(row_1,row_2,col_1,col_2)
        rm(list = c("row_0","row_1","row_2","col_0","col_1","col_2"))
      }
      if(identical(attr(model_to_Delta,"name"),"model_to_Delta_NP_ind")){
        # Default starting location for psi_0 is RL with uniform margin and 4/5 diagonal
        row_0 = rep(1/J,J)
        col_0 = rep(1/5,J)
        Delta_0 = diag(J)*(1-col_0) + outer(row_0,col_0,"*")
        Delta_1 = apply(Delta_0,2, function(d) softlog(d / d[length(d)]) )
        Delta_2 = apply(Delta_0,2, function(d) softlog(d / d[length(d)]) )
        psi_0 = c(c(Delta_1[-J,]),c(Delta_2[-J,]))
        rm(list = c("row_0","col_0","Delta_0","Delta_1","Delta_2"))
      }
      if(identical(attr(model_to_Delta,"name"),"model_to_Delta_NP")){
        # Default starting location for psi_0 is RL with uniform margin and 4/5 diagonal
        row_0 = rep(1/J,J)
        col_0 = rep(1/5,J)
        Delta_0 = diag(J)*(1-col_0) + outer(row_0,col_0,"*")
        # Building the joint distribution Y1, Y2 conditional on Y*
        Delta_0 = lapply(1:ncol(Delta_0), function(j) diag(Delta_0[j,]) %*% t(Delta_0))
        Delta_0 = do.call(cbind, Delta_0)
        Delta_0 = t(apply(Delta_0,1,function(d) softlog(d / d[length(d)])))
        psi_0 = c(Delta_0[,-J^2])
        rm(list = c("row_0","col_0","Delta_0"))
      }
    }

    # Defining the initial eta
    eta_0 = c(phi_0,psi_0)

    # Recording split_eta (the first position of psi in eta)
    split_eta = length(phi_0) + 1


    #------------------------------------------------------------
    # Defining the objective function for estimation
    #------------------------------------------------------------

    # Setting the penalty for violations of diagonal dominance
    lambda_dd = sum(tab$n)

    # Defining the objective function
    objective = function(eta){

      # Mapping model arguments to the entries of the joint distribution
      Pi_ = model_to_Pi(eta[1:(split_eta-1)], J)
      Delta_ = model_to_Delta(eta[(split_eta):length(eta)])

      # Computing the log likelihood of the data given Pi and Delta
      ll =  loglikelihood(theta = c(c(Pi_), c(Delta_)),
                          tab, J, K, lambda_dd)

      # Returning the log likelihood + penalties
      return(ll)

    }


    #------------------------------------------------------------
    # Estimation via MLE
    #------------------------------------------------------------

    if(mle){

      #------------------------------------------------------------
      # Optimizing objective w.r.t. eta
      #------------------------------------------------------------

      # Optimizing with standard quasi-newton
      mle_out = optim(
        par = eta_0,
        fn = objective,
        method = "BFGS",
        hessian = T,
        control = list(
          fnscale = -1, # Flipping the scale so that objective is maximized
          maxit = optim_maxit,
          reltol = optim_tol
        )
      )

      # Throwing an error if optim did not converge successfully
      if(mle_out$convergence == 1){stop("Error: Maximum iterations reached before numerical convergence, consider increasing optim max iterations or tolerance.")}
      if(mle_out$convergence > 1){stop(paste("Error: Numerical optimization failed with message:",
                                             mle_out$message, "... consider alternative `optim` settings."))}

      #------------------------------------------------------------
      # Testing the stability of the optimizer
      #------------------------------------------------------------

      # If check_stability, increase the number of alternative starting points
      if(check_stability){
        extra_starting_points = 9
      } else {
        extra_starting_points = 1
      }

      # Recording initial theta hat
      eta_hat1 = mle_out$par
      Pi_hat1 = model_to_Pi(eta_hat1[1:(split_eta-1)],J)
      Delta_hat1 = model_to_Delta(eta_hat1[split_eta:length(eta_hat1)])
      theta_hat1 = c(Pi_hat1,Delta_hat1)

      # Initializing total inconsistency
      inconsistency = 0

      for(draw in 1:extra_starting_points){

        # Additional estimates on perturbed starting locations
        mle_out2 = optim(par = eta_0 + rnorm(length(eta_0),sd = stability_sd),
                     fn = objective,
                     method = "BFGS",
                     control = list(fnscale = -1,
                                    maxit = optim_maxit,
                                    reltol = optim_tol))
        eta_hat2 = mle_out2$par
        Pi_hat2 = model_to_Pi(eta_hat2[1:(split_eta-1)],J)
        Delta_hat2 = model_to_Delta(eta_hat2[split_eta:length(eta_hat2)])
        theta_hat2 = c(Pi_hat2,Delta_hat2)

        # Adding the (absolute) difference between theta_hats
        inconsistency = inconsistency + sum(abs( theta_hat1 - theta_hat2 ))
      }

      if(inconsistency > 0.1){
        warning("Optimal eta is inconsitent across starting locations.")
      }

      #------------------------------------------------------------
      # Computing the covariance matrix of Pi
      #------------------------------------------------------------

      # Is the hessian of the fisher information invertible?
      fisher_info = -1*mle_out$hessian[1:(split_eta-1),1:(split_eta-1)]
      if (abs(det(fisher_info)) < 1e-9) {
        fisher_info_err = "Fisher information matrix is not invertible."
        warning("The Fisher information matrix is not invertible; using the Moore-Penrose inverse instead -- analytical variances may be inaccurate. Consider using Bayesian inference instead.")
      } else {
        fisher_info_err = "Fisher information matrix is invertible."
      }

      # Computing the Jacobian of model_to_Pi
      model_to_Pi_wrapper = function(phi){ return(model_to_Pi(phi,J)) }
      model_to_Pi_jacobian = numDeriv::jacobian(model_to_Pi_wrapper, mle_out$par[1:(split_eta-1)])

      # Computing the covariance matrix of Pi
      cov_Pi =  model_to_Pi_jacobian %*% pracma::pinv(fisher_info) %*% t(model_to_Pi_jacobian)

      # Forcing the diagonal of cov_Pi to be non-negative
      diag(cov_Pi) = pmax(diag(cov_Pi),0)

      #------------------------------------------------------------
      # Returning estimates and other info
      #------------------------------------------------------------

      # Objects to return
      Pi_hat_mle = model_to_Pi_wrapper(mle_out$par[1:(split_eta-1)])
      Delta_hat_mle = model_to_Delta(mle_out$par[(split_eta):(length(mle_out$par))])
      cov_Pi_mle = cov_Pi
      eta_hat_mle = mle_out$par
      log_likelihood_mle = mle_out$value
      optim_counts = mle_out$counts
      model_to_Pi_jacobian = model_to_Pi_jacobian
      eta_hessian_mle = mle_out$hessian
      fisher_info_err = fisher_info_err
      inconsistency_mle = inconsistency

    } else {
      # If MLE is not requested, returning NA for each output
      Pi_hat_mle = NA
      Delta_hat_mle = NA
      cov_Pi_mle = NA
      eta_hat_mle = NA
      log_likelihood_mle = NA
      optim_counts = NA
      model_to_Pi_jacobian = NA
      eta_hessian_mle = NA
      fisher_info_err = NA
      inconsistency_mle = NA
    }

    #------------------------------------------------------------
    # Bayesian inference w/ Gibbs sampler
    #------------------------------------------------------------

    if(bayesian){

      # Defining (the log of) the prior
      prior = function(eta){

        # computing the prior for Pi and Delta at phi and psi separately
        prior_Pi_val = log_prior_Pi(eta[1:(split_eta - 1)])
        prior_Delta_val = log_prior_Delta(eta[split_eta:length(eta)])

        # Returning the sum of the prior
        return(prior_Pi_val + prior_Delta_val)
      }

      #------------------------------------------------------------
      # Defining a gibbs sampler
      #------------------------------------------------------------

      # Defining an environment to keep track of the current parameter value and likelihood
      gibbs.env = new.env()
      gibbs.env$counter = 0
      gibbs.env$accepted_proposals = 0
      gibbs.env$eta_current = eta_0
      gibbs.env$ll_current = objective(eta_0) + prior(eta_0)
      gibbs.env$eta_history = list()
      gibbs.env$ll_history = list()

      # Defining a function to draw a proposal for a single parameter and accept/reject it
      gibbs_step = function(index){
        # Drawing the proposal
        eta_proposed = gibbs.env$eta_current
        gibbs_jump = rnorm(1, mean = 0, sd = gibbs_proposal_sd)
        eta_proposed[index] = eta_proposed[index] + gibbs_jump

        # Finding the log posterior likelihood
        ll_proposed = objective(eta_proposed) + prior(eta_proposed)

        # The MH acceptance probability is the difference of the proposed and current
        # log posterior likelihoods minus the gibbs jump
        if(log(runif(1)) <  ll_proposed - gibbs.env$ll_current - gibbs_jump){
          # Accept the proposal
          gibbs.env$eta_current = eta_proposed
          gibbs.env$ll_current = ll_proposed
          gibbs.env$accepted_proposals = gibbs.env$accepted_proposals + 1
        } # Nothing changes if the proposal is rejected
        # Nothing to return, all changes recorded in gibbs.env
      }

      #------------------------------------------------------------
      # Exploring the posterior of Pi and Delta via Gibbs-MCMC
      #------------------------------------------------------------

      while(gibbs.env$counter < n_mcmc_draws){

        # Increasing the step
        gibbs.env$counter = gibbs.env$counter+1

        # Looping through parameters, accepting or rejecting proposals in sequence
        invisible(sapply(seq_along(eta_0), gibbs_step))

        # Recording every Nth value
        if(gibbs.env$counter%%thinning_rate == 0){
          # Only printing progress in development versions...
          cat(paste0("steps = ",gibbs.env$counter,
                     "... accepted proposals = ", gibbs.env$accepted_proposals,
                     "... log likelihood = ", floor(gibbs.env$ll_current),"\n" ))
          gibbs.env$ll_history = append(gibbs.env$ll_history, list(gibbs.env$ll_current))
          gibbs.env$eta_history = append(gibbs.env$eta_history, list(gibbs.env$eta_current))
        }
      }

      #------------------------------------------------------------
      # Output for Bayesian inference
      #------------------------------------------------------------

      # Extracting the posterior for eta
      posterior_eta = do.call(rbind,gibbs.env$eta_history[(n_burnin/thinning_rate+1):length(gibbs.env$eta_history)])

      # Defining functions for extracting Pi and Delta from eta_hat
      eta_hat_to_Pi_hat = function(eta_hat, draw){

        # Extracting Pi from eta hat
        Pi_hat = model_to_Pi(eta_hat[1:(split_eta - 1)],J)

        # Binding with corresponding names of X and Y*
        X_name = sapply(X_names, function(v) rep(v,J)) |> c()
        Y_name = rep(Y_names, K) |> c()

        # Binding as a data frame
        Pi_hat_df = as.data.frame(cbind(Pi_hat, X_name, Y_name))

        if(!identical(X_vals,NA)){
          # Binding with corresponding names of X and Y*
          Pi_hat_df$X_val = sapply(X_vals, function(v) rep(v,J)) |> c()
          Pi_hat_df$Y_val = rep(Y_vals, K) |> c()
        } else {
          Pi_hat_df$X_val = NA
          Pi_hat_df$Y_val = NA
        }

        # Recording posterior draw
        Pi_hat_df$draw = draw

        return(Pi_hat_df)

      }

      eta_hat_to_Delta_hat = function(eta_hat, draw){

        # Extracting Delta from eta hat
        Delta_hat = model_to_Delta(eta_hat[split_eta:length(eta_hat)])

        # Binding to corresponding names of Y*, Y1, Y2
        Y2_name = rep(Y_names, J^2)
        Y1_name = rep(c(sapply(Y_names, function(y) rep(y,J))),J)
        Ys_name = c(sapply(Y_names, function(y) rep(y, J^2)))

        # Binding as a data frame
        Delta_hat_df = as.data.frame(cbind(Delta_hat, Ys_name, Y1_name, Y2_name))

        if(!identical(Y_vals,NA)){
          # Binding with corresponding names of X and Y*
          Delta_hat_df$Y2_val = rep(Y_vals, J^2)
          Delta_hat_df$Y1_val = rep(c(sapply(Y_vals, function(y) rep(y,J))),J)
          Delta_hat_df$Ys_val = c(sapply(Y_vals, function(y) rep(y, J^2)))
        } else {
          Delta_hat_df$Ys_val = NA
          Delta_hat_df$Y1_val = NA
          Delta_hat_df$Y2_val = NA
        }

        # Recording posterior draw
        Delta_hat_df$draw = draw

        return(Delta_hat_df)

      }

      # Building the posterior for Pi
      posterior_Pi = do.call(rbind,lapply(1:nrow(posterior_eta), function(d) eta_hat_to_Pi_hat(posterior_eta[d,],d)))

      # Building the posterior for Delta
      posterior_Delta = do.call(rbind,lapply(1:nrow(posterior_eta), function(d) eta_hat_to_Delta_hat(posterior_eta[d,],d)))

      # Renaming draws to match MCMC history
      posterior_Pi$draw = n_burnin + thinning_rate*posterior_Pi$draw
      posterior_Delta$draw = n_burnin + thinning_rate*posterior_Delta$draw

      # Recording the history of log likelihood
      ll_history = unlist(gibbs.env$ll_history) |> as.data.frame()
      colnames(ll_history) = "ll"
      ll_history$draw = thinning_rate*c(1:nrow(ll_history))

      # Recording
      accepted_proposals = gibbs.env$accepted_proposals

      #------------------------------------------------------------
      # Trace plots
      #------------------------------------------------------------

      if(makeplots){
        # Defining a function to make trace plots
        traceplots = function(posterior_draws, posterior_name){
          # Forming as a data.frame
          trace_df = as.data.frame(posterior_draws)
          trace_df$draw = n_burnin + thinning_rate*c(1:nrow(trace_df))
          # Making the plot
          trace_plot = ggplot2::ggplot(trace_df, ggplot2::aes(x = draw, y = posterior_draws)) +
            ggplot2::geom_line(linewidth = 0.5, col = "#637384") +
            ggplot2::labs(x = "MCMC draw", y = posterior_name) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              axis.ticks = ggplot2::element_blank(),
              axis.title = ggplot2::element_text(size=12)
            )
          return(trace_plot)
        }

        # Trace plots for eta
        trace_plots_eta = lapply(1:ncol(posterior_eta), function(j) traceplots(posterior_eta[,j], paste0("eta_",j)))
        names(trace_plots_eta) = paste0("eta_",1:ncol(posterior_eta))

        # Trace plots for Pi
        X_names_plot = unique(posterior_Pi$X_name)
        Y_names_plot = unique(posterior_Pi$Y_name)
        XY_names = expand.grid(X_names_plot,Y_names_plot)
        colnames(XY_names) = c("X","Y")
        XY_names$label = paste0( "Pr(X = ",XY_names$X, ", Y = ", XY_names$Y, ")")
        trace_plots_Pi = lapply(1:nrow(XY_names), function(j)
          traceplots(subset(posterior_Pi, X_name == XY_names$X[j] & Y_name == XY_names$Y[j])$Pi_hat,
                     XY_names$label[j]))
        names(trace_plots_Pi) = XY_names$label

        # Trace plots for Delta
        YYY_names = expand.grid(Y_names_plot, Y_names_plot, Y_names_plot)
        colnames(YYY_names) = c("Y1","Y2","Y*")
        YYY_names$label = paste0( "Pr(Y1 = ",YYY_names$Y1, ", Y2 = ", YYY_names$Y2, " | Y* = ", YYY_names$`Y*`, ")")
        trace_plots_Delta = lapply(1:nrow(YYY_names), function(j)
          traceplots(subset(posterior_Delta, Y1_name == YYY_names$Y1[j] & Y2_name == YYY_names$Y2[j] & Ys_name == YYY_names$`Y*`[j] )$Delta_hat,
                     YYY_names$label[j]))
        names(trace_plots_Delta) = YYY_names$label

      } else {
        # If plots are not requested, returning NA
        trace_plots_eta = NA
        trace_plots_Pi = NA
        trace_plots_Delta = NA
      }

    } else {
      # If Bayesian inference is not requested, returning NA for all objects
      posterior_eta = NA
      posterior_Pi = NA
      posterior_Delta = NA
      ll_history = NA
      accepted_proposals = NA
      trace_plots_eta = NA
      trace_plots_Pi = NA
      trace_plots_Delta = NA
    }

    #------------------------------------------------------------
    # Returning estimates and other info
    #------------------------------------------------------------

    # Return results
    return(list(
      Pi_hat_mle = Pi_hat_mle,
      Delta_hat_mle = Delta_hat_mle,
      cov_Pi_mle = cov_Pi_mle,
      eta_hat_mle = eta_hat_mle,
      log_likelihood_mle = log_likelihood_mle,
      optim_counts = optim_counts,
      model_to_Pi_jacobian = model_to_Pi_jacobian,
      eta_hessian_mle = eta_hessian_mle,
      fisher_info_err = fisher_info_err,
      inconsistency_mle = inconsistency_mle,
      posterior_Pi = posterior_Pi,
      posterior_Delta = posterior_Delta,
      posterior_eta = posterior_eta,
      ll_history = ll_history,
      accepted_proposals = accepted_proposals,
      trace_plots_eta = trace_plots_eta,
      trace_plots_Pi = trace_plots_Pi,
      trace_plots_Delta = trace_plots_Delta,
      tab = tab,
      J = J,
      K = K,
      X_names = X_names,
      Y_names = Y_names,
      X_vals = X_vals,
      Y_vals = Y_vals,
      W_weight = sum(tab$n),
      phi_0 = phi_0,
      psi_0 = psi_0,
      eta_0 = eta_0
    ))

  }

  #------------------------------------------------------------
  # Estimating Pi and Delta
  #------------------------------------------------------------

  # Is there more than one covariate value?
  if(class(tab)=="list"){

    # Setting up parallel processing
    workers = parallel::makeCluster(cores)

    # Teaching the workers misclassifyr
    parallel::clusterEvalQ(workers, require(misclassifyr)) |> invisible()

    # Exporting common objects to workers
    common_objects = c("estimate_misclassification", "makeplots","mle","optim_maxit","optim_tol","check_stability","stability_sd",
                       "bayesian","n_mcmc_draws", "n_burnin","thinning_rate", "gibbs_proposal_sd",
                       "log_prior_Pi", "log_prior_Delta")
    for (obj in common_objects) { # exporting workers from the local environment one by one
      parallel::clusterExport(workers, varlist = obj, envir = environment()) |> invisible()
    }

    # Drawing from the posterior of Pi and Delta within covariate cells
    misclassification_out = pbapply::pblapply(
      misclassification_inputs,
      estimate_misclassification,
      cl = workers
    )

    #------------------------------------------------------------
    # Restructuring misclassification_out into lists of each output
    #------------------------------------------------------------

    # Recording the names of each type of output
    output_names = names(misclassification_out[[1]])

    # Defining a quick function to rebundle output into one list across covariate cells
    rebundle_output = function(index){
      output = lapply(misclassification_out, function(out) out[[index]])
      names(output) = unlist(W_names) # Adding covariate cell names
      return(output)
    }

    # Rebundling across covariate cells
    misclassification_output = lapply(seq_along(output_names), rebundle_output)

    # Adding back output names
    names(misclassification_output) = output_names

  } else {

    # Estimating Pi and Delta for the full population
    misclassification_output = estimate_misclassification(misclassification_inputs)

  }



  #------------------------------------------------------------
  # Generating Pi and Delta figures
  #------------------------------------------------------------

  if(mle){ # Are there MLE estimates?

    if(makeplots){ # Should plots be generated?

      if(class(misclassification_output$Pi_hat_mle) == "list"){ # Aggregating across lists of Pi_hat and Delta_hat

        # If the dimension of Pi is inconsistent across covariate cells, not plotting
        if(length(unique(paste(unlist(misclassification_output$J),
                               unlist(misclassification_output$K)))) == 1){

          # Normalizing the weight of each cell
          W_weights = unlist(misclassification_output$W_weight)
          W_weights = W_weights / sum(W_weights)

          # Averaging Pi across covariate cells
          Pi_hat_mle_plot_df = do.call(cbind, misclassification_output$Pi_hat_mle) %*% W_weights |> c() |> as.data.frame()
          colnames(Pi_hat_mle_plot_df) = c("Pi_hat")
          Pi_hat_mle_plot_df$X = c(sapply(misclassification_output$X_names[[1]], function(xn) rep(xn, misclassification_output$J[[1]])))
          Pi_hat_mle_plot_df$Y = rep(misclassification_output$Y_names[[1]], misclassification_output$J[[1]])

          # Averaging Delta across covariate cells
          Delta_hat_mle_plot_df = do.call(
            cbind,
            lapply(misclassification_output$Delta_hat_mle,
                   function(Delta) apply(
                     matrix(Delta, nrow = misclassification_output$J[[1]]),
                     2, sum))
            ) %*% W_weights |>c() |> as.data.frame()
          colnames(Delta_hat_mle_plot_df) = c("Delta_hat")
          Delta_hat_mle_plot_df$Ys = c(sapply(misclassification_output$Y_names[[1]], function(xn) rep(xn, misclassification_output$J[[1]])))
          Delta_hat_mle_plot_df$Y1 = rep(misclassification_output$Y_names[[1]], misclassification_output$J[[1]])

          # Plotting the joint distribution of X and Y*
          Pi_hat_mle_plot = ggplot2::ggplot(Pi_hat_mle_plot_df,
                                            ggplot2::aes(x = X,y = Y, fill = Pi_hat ) ) +
            ggplot2::geom_raster(show.legend = F) +
            ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676") +
            ggplot2::geom_text(ggplot2::aes(label = round(Pi_hat,2),
                          color = Pi_hat > median(Pi_hat) )) +
            ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
            ggplot2::labs(x = X_col_name, y = Y_col_name) +
            ggplot2::theme_minimal() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           axis.text = ggplot2::element_text(size = 11),
                           axis.title = ggplot2::element_text(size = 12))

          # Plotting the conditional distribution Y1 | Y*
          Delta_hat_mle_plot = ggplot2::ggplot(Delta_hat_mle_plot_df,
                                               ggplot2::aes(x = Ys,y = Y1, fill = Delta_hat)) +
            ggplot2::geom_raster(show.legend = F) +
            ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676", limits= c(0,1)) +
            ggplot2::geom_text(ggplot2::aes(label = round(Delta_hat,2),
                                            color = Delta_hat > 0.5 )) +
            ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
            ggplot2::labs(x = paste("Latent", Y_col_name), y = paste("Observed", Y_col_name)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           axis.text = ggplot2::element_text(size = 11),
                           axis.title = ggplot2::element_text(size = 12))

        } else {
          Pi_hat_mle_plot = NA
          Delta_hat_mle_plot = NA
        }

      } else {

        # Averaging Pi across covariate cells
        Pi_hat_mle_plot_df = c(misclassification_output$Pi_hat_mle) |> as.data.frame()
        colnames(Pi_hat_mle_plot_df) = c("Pi_hat")
        Pi_hat_mle_plot_df$X = c(sapply(misclassification_output$X_names, function(xn) rep(xn, misclassification_output$J)))
        Pi_hat_mle_plot_df$Y = rep(misclassification_output$Y_names, misclassification_output$J)

        # Averaging Delta across covariate cells
        Delta_hat_mle_plot_df = apply( matrix(misclassification_output$Delta_hat_mle, nrow = misclassification_output$J),2, sum) |> c() |> as.data.frame()
        colnames(Delta_hat_mle_plot_df) = c("Delta_hat")
        Delta_hat_mle_plot_df$Ys = c(sapply(misclassification_output$Y_names, function(xn) rep(xn, misclassification_output$J)))
        Delta_hat_mle_plot_df$Y1 = rep(misclassification_output$Y_names, misclassification_output$J)

        # Plotting the joint distribution of X and Y*
        Pi_hat_mle_plot = ggplot2::ggplot(Pi_hat_mle_plot_df,
                                          ggplot2::aes(x = X,y = Y, fill = Pi_hat ) ) +
          ggplot2::geom_raster(show.legend = F) +
          ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676") +
          ggplot2::geom_text(ggplot2::aes(label = round(Pi_hat,2),
                                          color = Pi_hat > median(Pi_hat) )) +
          ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
          ggplot2::labs(x = X_col_name, y = Y_col_name) +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         axis.text = ggplot2::element_text(size = 11),
                         axis.title = ggplot2::element_text(size = 12))

        # Plotting the conditional distribution Y1 | Y*
        Delta_hat_mle_plot = ggplot2::ggplot(Delta_hat_mle_plot_df,
                                             ggplot2::aes(x = Ys,y = Y1, fill = Delta_hat)) +
          ggplot2::geom_raster(show.legend = F) +
          ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676", limits= c(0,1)) +
          ggplot2::geom_text(ggplot2::aes(label = round(Delta_hat,2),
                                          color = Delta_hat > 0.5 )) +
          ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
          ggplot2::labs(x = paste("Latent", Y_col_name), y = paste("Observed", Y_col_name)) +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         axis.text = ggplot2::element_text(size = 11),
                         axis.title = ggplot2::element_text(size = 12))
      }

    }

  } else {

      Pi_hat_mle_plot = NA
      Delta_hat_mle_plot = NA

    }

  if(bayesian){ # Are there Bayesian estimates?

    if(makeplots){ # Should plots be generated?

      if(class(misclassification_output$posterior_Pi) == "list"){ # Aggregating across lists of Pi_hat and Delta_hat

        # If the dimension of Pi is inconsistent across covariate cells, not plotting
        if(length(unique(paste(unlist(misclassification_output$J),
                               unlist(misclassification_output$K)))) == 1){

          # Normalizing the weight of each cell
          W_weights = unlist(misclassification_output$W_weight)
          W_weights = W_weights / sum(W_weights)

          # Finding posterior mean Pi within covariate cells
          posterior_Pi_mean = lapply(misclassification_output$posterior_Pi,
                                    function(Pi_draws) Pi_draws |>
                                      dplyr::group_by(X_name, Y_name) |>
                                      dplyr::summarise(Pi_hat = mean(Pi_hat),
                                                       .groups = "drop") |>
                                      as.data.frame())
          posterior_Pi_mean_plot_df = posterior_Pi_mean[[1]][,c("X_name","Y_name")]
          # Averaging Pi across covariate cells
          posterior_Pi_mean_plot_df$Pi_hat = do.call(cbind,lapply(posterior_Pi_mean, function(Pi_mean) Pi_mean$Pi_hat)) %*% W_weights


          # Finding posterior mean Delta within covariate cells
          posterior_Delta_mean = lapply(
            lapply(misclassification_output$posterior_Delta,
                   function(Delta_draws) Delta_draws |>
                     dplyr::group_by(Ys_name, Y1_name, draw) |>
                     dplyr::summarise(Delta_hat = sum(Delta_hat),
                                      .groups = "drop") |>
                     as.data.frame()),
            function(Delta_draws) Delta_draws |>
              dplyr::group_by(Ys_name, Y1_name) |>
              dplyr::summarise(Delta_hat = mean(Delta_hat),
                               .groups = "drop") |>
              as.data.frame())

          # Recording margin names
          posterior_Delta_mean_plot_df = posterior_Delta_mean[[1]][,c("Ys_name","Y1_name")]

          # Averaging Delta across covariate cells
          posterior_Delta_mean_plot_df$Delta_hat = do.call(cbind,lapply(posterior_Delta_mean, function(Delta_mean) Delta_mean$Delta_hat)) %*% W_weights

          # Plotting the joint distribution of X and Y*
          Pi_hat_posterior_plot = ggplot2::ggplot(posterior_Pi_mean_plot_df,
                                            ggplot2::aes(x = X_name,y = Y_name, fill = Pi_hat ) ) +
            ggplot2::geom_raster(show.legend = F) +
            ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676") +
            ggplot2::geom_text(ggplot2::aes(label = round(Pi_hat,2),
                                            color = Pi_hat > median(Pi_hat) )) +
            ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
            ggplot2::labs(x = X_col_name, y = Y_col_name) +
            ggplot2::theme_minimal() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           axis.text = ggplot2::element_text(size = 11),
                           axis.title = ggplot2::element_text(size = 12))

          # Plotting the conditional distribution Y1 | Y*
          Delta_hat_posterior_plot = ggplot2::ggplot(posterior_Delta_mean_plot_df,
                                               ggplot2::aes(x = Ys_name,y = Y1_name, fill = Delta_hat)) +
            ggplot2::geom_raster(show.legend = F) +
            ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676", limits= c(0,1)) +
            ggplot2::geom_text(ggplot2::aes(label = round(Delta_hat,2),
                                            color = Delta_hat > 0.5 )) +
            ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
            ggplot2::labs(x = paste("Latent", Y_col_name), y = paste("Observed", Y_col_name)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           axis.text = ggplot2::element_text(size = 11),
                           axis.title = ggplot2::element_text(size = 12))

        } else {
          Pi_hat_posterior_plot = NA
          Delta_hat_posterior_plot = NA
        }

      } else {

        # Finding posterior mean Pi within covariate cells
        posterior_Pi_mean_plot_df = misclassification_output$posterior_Pi |>
          dplyr::group_by(X_name, Y_name) |>
          dplyr::summarise(Pi_hat = mean(Pi_hat), .groups = "drop") |>
          as.data.frame()

        # Finding posterior mean Delta within covariate cells
        posterior_Delta_mean_plot_df = lapply(
          misclassification_output$posterior_Delta,
          function(Delta_draws) Delta_draws |>
            dplyr::group_by(Ys_name, Y1_name, draw) |>
            dplyr::summarise(Delta_hat = sum(Delta_hat), .groups = "drop")) |>
            dplyr::group_by(Ys_name, Y1_name) |>
            dplyr::summarise(Delta_hat = mean(Delta_hat),
                             .groups = "drop") |>
            as.data.frame()

        # Plotting the joint distribution of X and Y*
        Pi_hat_posterior_plot = ggplot2::ggplot(posterior_Pi_mean_plot_df,
                                                ggplot2::aes(x = X_name,y = Y_name, fill = Pi_hat ) ) +
          ggplot2::geom_raster(show.legend = F) +
          ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676") +
          ggplot2::geom_text(ggplot2::aes(label = round(Pi_hat,2),
                                          color = Pi_hat > median(Pi_hat) )) +
          ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
          ggplot2::labs(x = X_col_name, y = Y_col_name) +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         axis.text = ggplot2::element_text(size = 11),
                         axis.title = ggplot2::element_text(size = 12))

        # Plotting the conditional distribution Y1 | Y*
        Delta_hat_posterior_plot = ggplot2::ggplot(posterior_Delta_mean_plot_df,
                                                   ggplot2::aes(x = Ys_name,y = Y1_name, fill = Delta_hat)) +
          ggplot2::geom_raster(show.legend = F) +
          ggplot2::scale_fill_gradient(low="#EDF3FF", high="#002676", limits= c(0,1)) +
          ggplot2::geom_text(ggplot2::aes(label = round(Delta_hat,2),
                                          color = Delta_hat > 0.5 )) +
          ggplot2::scale_color_manual(guide = "none", values = c("black", "white")) +
          ggplot2::labs(x = paste("Latent", Y_col_name), y = paste("Observed", Y_col_name)) +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         axis.text = ggplot2::element_text(size = 11),
                         axis.title = ggplot2::element_text(size = 12))

      }

    } else {
      # If plots are not requested, returning NA
      Pi_hat_posterior_plot = NA
      Delta_hat_posterior_plot = NA
    }

  } else {
    # If Bayesian inference is not requested, returning NA
    Pi_hat_posterior_plot = NA
    Delta_hat_posterior_plot = NA

  }

  #------------------------------------------------------------
  # Estimating beta and SE
  #------------------------------------------------------------

  if(estimate_beta){

    # Computing beta from Pi for MLE
    if(mle){

      # Extracting W_weights and normalizing
      W_weights = unlist(misclassification_output$W_weight)
      W_weights = W_weights/sum(W_weights)

      # Estimating beta given Pi
      beta_hat_mle = Pi_to_beta(misclassification_output$Pi_hat_mle, X_vals, Y_vals, W_weights)
      names(beta_hat_mle) = "MLE beta"

      # Estimating the SE via the delta method
      se_beta_mle = se_beta_deltamethod(
        Pi = misclassification_output$Pi_hat_mle,
        cov_Pi = misclassification_output$cov_Pi_mle,
        X_vals,
        Y_vals,
        W_weights
      )
      names(se_beta_mle) = "SE for MLE beta"

      # If controls are used, returning a list of betas within each cell
      if(class(tab) == "list"){

        # Estimating beta within each control cell
        betas_hat_mle = sapply(seq_along(W_names), function(j)
          Pi_to_beta(misclassification_output$Pi_hat_mle[[j]], X_vals[[j]], Y_vals[[j]], 1))
        names(betas_hat_mle) = paste("MLE beta in cell", W_names)

        # Estimating the SE in each control cell
        se_betas_mle = sapply(seq_along(W_names), function(j)
          se_beta_deltamethod(misclassification_output$Pi_hat_mle[[j]], misclassification_output$cov_Pi_mle[[j]], X_vals[[j]], Y_vals[[j]], 1))
        names(se_betas_mle) = paste("SE for MLE beta in cell", W_names)

      } else {
        # If there are no controls, returning NA
        betas_mle = NA
        se_betas_mle = NA

      }

    } else {
      # If MLE is not computed, returning NAs
      beta_hat_mle = NA
      se_beta_mle = NA
      betas_hat_mle = NA
      se_betas_mle = NA
    }

    if(bayesian){

      # Aggregating across control cells
      if(class(tab) == "list"){

        # Extracting W_weights and normalizing
        W_weights = unlist(misclassification_output$W_weight)
        W_weights = W_weights/sum(W_weights)

        # Defining a quick function to aggregate the posterior across control cells
        posterior_agg = function(d){
          # Grabbing the dth draw of the posterior
          posterior_df_list = lapply(seq_along(W_names), function(w)
            subset(misclassification_output$posterior_Pi[[w]], draw == d))

          # Scaling weights to reflect control cell size
          posterior_df_list = Map(function(post_df, W_weight) {
            post_df$Pi_hat = post_df$Pi_hat * W_weight
            return(post_df)
          }, posterior_df_list, W_weights)

          # Binding together and returning
          posterior_df = do.call(rbind, posterior_df_list)

          return(posterior_df)
        }

        posterior_beta = sapply(unique(misclassification_output$posterior_Pi[[1]]$draw), function(d)
          lm(Y_val ~ X_val,
             data = posterior_agg(d),
             weight = posterior_agg(d)$Pi_hat
          )$coefficients[2] |> unname())


        # Recording the posterior of beta within covariate cells
        posterior_betas = lapply(seq_along(W_names), function(j)
          sapply(unique(misclassification_output$posterior_Pi[[j]]$draw ), function(d)
            lm(Y_val ~ X_val,
               data = subset(misclassification_output$posterior_Pi[[j]], draw == d),
               weight = subset(misclassification_output$posterior_Pi[[j]], draw == d)$Pi_hat
            )$coefficients[2] |> unname()
          )
        )
        names(posterior_betas) = W_names

        # Computing Chen, Christensen, and Tamer Partial ID-robust CI

        # Summing the likelihood history across covariate cells
        ll_history = sapply(unique(misclassification_output$posterior_Pi[[1]]$draw), function(d) {
          # For each draw 'd', loop over all control cells 'w' and sum the likelihood
          sum(sapply(seq_along(W_names), function(w) {
            # Subset the ll_history for the current control cell 'w' and draw 'd'
            subset(misclassification_output$ll_history[[w]], draw == d)$ll
          }))
        })

        # sort ll_history to find top 95%
        ll_sorted = sort(ll_history, decreasing = F)
        ll_critical = ll_sorted[floor(length(ll_sorted)*0.05) ]
        CCT_draws = unique(misclassification_output$posterior_Pi[[1]]$draw[ll_history > ll_critical])
        CCTCI = c(min(posterior_beta[ll_history > ll_critical]), max(posterior_beta[ll_history > ll_critical]) )

        # Extracting the median and the sd
        posterior_beta_med = median(posterior_beta)
        posterior_beta_sd = sd(posterior_beta)
        posterior_betas_med = lapply(posterior_betas, median)
        posterior_betas_sd = lapply(posterior_betas, sd)

      } else {

        # Estimating beta for each draw of the posterior
        posterior_beta = sapply(unique(misclassification_output$posterior_Pi$draw ), function(d)
          lm(Y_val ~ X_val,
             data = subset(misclassification_output$posterior_Pi, draw == d),
             weight = subset(misclassification_output$posterior_Pi, draw == d)$Pi_hat
             )$coefficients[2] |> unname()
          )

        # Computing Chen, Christensen, and Tamer Partial ID-robust CI

        # Summing the likelihood history across covariate cells
        ll_history = sapply(unique(misclassification_output$posterior_Pi$draw), function(d) {
          # For each draw 'd' sum the likelihood
          sum(subset(misclassification_output$ll_history, draw == d)$ll)
        })

        # sort ll_history to find top 95%
        ll_sorted = sort(ll_history, decreasing = F)
        ll_critical = ll_sorted[floor(length(ll_sorted)*0.05) ]
        CCT_draws = unique(misclassification_output$posterior_Pi$draw[ll_history > ll_critical])
        CCTCI = c(min(posterior_beta[ll_history > ll_critical]), max(posterior_beta[ll_history > ll_critical]) )

        # Recording the median and SD
        posterior_beta_med = median(posterior_beta)
        posterior_beta_sd = sd(posterior_beta)

        # Returning NA if there are not multiple control cells
        posterior_betas = NA
        posterior_betas_med = NA
        posterior_betas_sd = NA

      }


    } else {

      posterior_beta = NA
      posterior_beta_med = NA
      posterior_beta_sd = NA
      posterior_betas = NA
      posterior_betas_med = NA
      posterior_betas_sd = NA
      CCT_draws = NA
      CCTCI = NA

    }


  } else {

    # If beta is not computed, returning NAs
    beta_hat_mle = NA
    se_beta_mle = NA
    betas_hat_mle = NA
    se_betas_mle = NA
    posterior_beta = NA
    posterior_beta_med = NA
    posterior_beta_sd = NA
    posterior_betas = NA
    posterior_betas_med = NA
    posterior_betas_sd = NA
    CCT_draws = NA
    CCTCI = NA

  }


  #------------------------------------------------------------
  # Returning a list of all outputs
  #------------------------------------------------------------

  return(
    append(
      misclassification_output,
      list(
        Pi_hat_mle_plot = Pi_hat_mle_plot,
        Delta_hat_mle_plot = Delta_hat_mle_plot,
        Pi_hat_posterior_plot = Pi_hat_posterior_plot,
        Delta_hat_posterior_plot = Delta_hat_posterior_plot,
        beta_hat_mle = beta_hat_mle,
        se_beta_mle = se_beta_mle,
        betas_hat_mle = betas_hat_mle,
        se_betas_mle = se_betas_mle,
        posterior_beta = posterior_beta,
        posterior_beta_med = posterior_beta_med,
        posterior_beta_sd = posterior_beta_sd,
        posterior_betas = posterior_betas,
        posterior_betas_med = posterior_betas_med,
        posterior_betas_sd = posterior_betas_sd,
        CCT_draws = CCT_draws,
        CCTCI = CCTCI
      )
    )
  )
}
