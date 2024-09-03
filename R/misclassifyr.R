#' misclassifyr
#'
#' This function provides a menu of options for estimation and inference of misclassification models in which the analyst has access to two noisy measures, `Y1` and `Y2` of a latent outcome `Y*`, a correctly measured covariate `X`, and discrete controls `W`.
#'
#' @import parallel
#'
#' @param tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should be numeric with names `Y1`, `Y2`, `X`, and `n` where `Y1` and `Y2` take each value between `1` and `J`, `X` takes each value between `1` and `K`, and
#' @param J An integer or list corresponding to the number of unique values of `Y1` and `Y2`.
#' @param K An integer or list corresponding to the number of unique values of `X`.
#' @param model_to_Pi A function or list of functions mapping the parameters of a model for the joint distribution to the joint distribution \Pi
#' @param model_to_Delta A function or list of functions mapping the parameters of a model to the conditional distribution Y1, Y2 | Y*, \Delta
#' @param phi_0 A numeric vector or list of numeric vectors providing the starting location for optimization for the argument to model_to_Pi.
#' @param psi_0 A numeric vector or list of numeric vectors providing the starting location for optimization for the argument to model_to_Delta.
#' @param split_eta An integer or list indicating where to split the vector `eta` in `phi` and `psi`, the arguments to `model_to_Pi` and `model_to_Delta` respectively.
#' @param X_names A character vector or list corresponding to the values of the regressor X.
#' @param Y_names A character vector or list corresponding to the values of the outcome Y.
#' @param W_names A character vector corresponding to the values of the control W in each cell.
#' @param estimate_beta A logical value indicating whether to regress Y on X.
#' @param estimate_betas A logical value indicating whether to regress Y on X within covariate cells.
#' @param X_vals A numeric vector or list of numeric vectors providing the values of X associated with the columns of Pi.
#' @param Y_vals A numeric vector or list of numeric vectors providing the values of Y associated with the rows of Pi.
#' @param n_mcmc_draws An integer corresponding to the length of the MCMC chain.
#' @param n_burnin An integer giving the length of the burn-in period for each MCMC chain, must be shorter than `n_mcmc_draws`.
#' @param thinning_rate An integer indicating how frequently to record posterior draws from the MCMC chain -- e.g. a `thinning_rate` of 2 records every other draw.
#' @param cores An integer for the number of CPUs available for parallel processing.
#' @return An object that includes estimates and information from the estimation process
#' @export
misclassifyr <- function(
    tab,
    J,
    K,
    model_to_Pi = "model_to_Pi_NP",
    model_to_Delta = "model_to_Delta_NP_ind",
    phi_0 = NA,
    psi_0 = NA,
    X_names = NA,
    Y_names = NA,
    W_names = NA,
    estimate_beta = F,
    estimate_betas = F,
    X_vals = NA,
    Y_vals = NA,
    n_mcmc_draws = 1e5,
    n_burnin = 5e4,
    thinning_rate = 5,
    cores = 1) {

  #------------------------------------------------------------
  # Catching errors in some variables
  # other input errors caught in MisclassMLE()
  #------------------------------------------------------------

  if(!(is.integer(cores)|(is.numeric(cores) & abs(cores - floor(cores)) < 1e-16))){
    stop("Error: `cores` should be an integer.")
  } else if(parallel::detectCores() < cores){
    stop("Error: You requested more cores than appear available on your machine.")
  }

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
      # Returning a warning that some objects will be copied.
      warning(paste0("The following objects were not provided as a list and will be copied across covariate cells: ",
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
  # Converting inputs to a list (or list of lists)  for MisclassMLE
  #-----------------------------

  if(class(tab) == "list"){
    # rebundling each list as a list of lists
    MisclassMLE_inputs = lapply(seq_along(tab), function(j) list(
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
    # Bundling inputs into a list for MisclassMLE
    MisclassMLE_inputs = list(
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
  # Defining a function MisclassMLE to estimate Pi and Delta
  #------------------------------------------------------------

  MisclassMLE = function(MisclassMLE_input){

    #------------------------------------------------------------
    # Unpacking inputs
    #------------------------------------------------------------

    tab = MisclassMLE_input$tab
    J = MisclassMLE_input$J
    K = MisclassMLE_input$K
    model_to_Pi = MisclassMLE_input$model_to_Pi
    model_to_Delta = MisclassMLE_input$model_to_Delta
    phi_0 = MisclassMLE_input$phi_0
    psi_0 = MisclassMLE_input$psi_0
    X_names = MisclassMLE_input$X_names
    Y_names = MisclassMLE_input$Y_names
    W_names = MisclassMLE_input$W_names
    X_vals = MisclassMLE_input$X_vals
    Y_vals = MisclassMLE_input$Y_vals

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
    # Setting the starting location for MCMC
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
        phi_0 = softlog(tab_xy$p[1:(J*K-1)])

        # Default starting location for phi_0 is flat
        #phi_0 = softlog(rep(1/(J*K), J*K-1 ))
      }
    }

    if(identical(psi_0,NA)){
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

    # Recording split_eta (the first position of psi in eta)
    split_eta = length(phi_0) + 1


    #------------------------------------------------------------
    # Defining the objective function for estimation
    #------------------------------------------------------------

    # Defining the objective function
    objective = function(eta){

      # Mapping model arguments to the entries of the joint distribution
      Pi_ = model_to_Pi(eta[1:(split_eta-1)], J, K)
      Delta_ = model_to_Delta(eta[(split_eta):length(eta)])

      # Computing the log likelihood of the data given Pi and Delta
      ll =  loglikelihood(theta = c(c(Pi_), c(Delta_)),
                          tab, J, K, gibbs.env$lambda_pos, gibbs.env$lambda_dd)

      # Returning the log likelihood + penalties
      return(ll)

    }

    #------------------------------------------------------------
    # Defining the prior density
    #------------------------------------------------------------

    prior = function(eta){

      # Mapping model arguments to the entries of the joint distribution
      Pi_ = model_to_Pi(eta[1:(split_eta-1)], J, K)
      Delta_ = model_to_Delta(eta[(split_eta):length(eta)])

      # A flat prior in levels implies an exponential density in logs,
      # and a flat density in logs (so just returning the sum of logs)
      return(sum(softlog(c(Pi_,Delta_))))

    }

    #------------------------------------------------------------
    # Defining a gibbs sampler
    #------------------------------------------------------------

    # Defining an environment to keep track of the current parameter value and likelihood
    gibbs.env = new.env()
    gibbs.env$counter = 0
    gibbs.env$lambda_pos_0 = sum(tab$n)
    gibbs.env$lambda_dd_0 = sum(tab$n)
    gibbs.env$lambda_pos_T = sum(tab$n)^1.5
    gibbs.env$lambda_dd_T = sum(tab$n)^1.5
    gibbs.env$lambda_pos = sum(tab$n)
    gibbs.env$lambda_dd = sum(tab$n)
    gibbs.env$accepted_proposals = 0
    gibbs.env$eta_current = eta_0
    gibbs.env$eta_history = list()
    gibbs.env$ll_history = list()

    # Defining a function to draw a proposal for a single parameter and accept/reject it
    gibbs_step = function(index){
      # Drawing the proposal
      eta_proposed = gibbs.env$eta_current
      gibbs_jump = rnorm(1,0,0.1)
      eta_proposed[index] = eta_proposed[index] + gibbs_jump

      # Finding the log posterior likelihood
      ll_proposed = objective(eta_proposed) + prior(eta_proposed)
      ll_current = objective(gibbs.env$eta_current) + prior(gibbs.env$eta_current)

      # The MH acceptance probability is the difference of the proposed and current
      # log posterior likelihoods minus the gibbs jump
      if(log(runif(1)) <  ll_proposed - ll_current - gibbs_jump){
        # Accept the proposal
        gibbs.env$eta_current = eta_proposed
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

      # Increasing the scale of penalties
      penalty_size = 1 / (1 + exp(-1*(100/n_mcmc_draws)*(gibbs.env$counter - n_mcmc_draws/10)))
      gibbs.env$lambda_pos = (1-penalty_size)*gibbs.env$lambda_pos_0 + penalty_size*gibbs.env$lambda_pos_T
      gibbs.env$lambda_dd = (1-penalty_size)*gibbs.env$lambda_dd_0 + penalty_size*gibbs.env$lambda_dd_T

      # Looping through parameters, accepting or rejecting proposals in sequence
      invisible(sapply(seq_along(eta_0), gibbs_step))

      # Recording new LL evaluated at the max penalty
      gibbs.env$lambda_pos = gibbs.env$lambda_pos_T
      gibbs.env$lambda_dd = gibbs.env$lambda_dd_T
      ll_current = objective(gibbs.env$eta_current) + prior(gibbs.env$eta_current)

      # Recording every Nth value
      if(gibbs.env$counter%%thinning_rate==1){
        cat(paste0("steps = ",gibbs.env$counter,
                   "... accepted proposals = ", gibbs.env$accepted_proposals,
                   "... log likelihood = ", floor(ll_current),"\n" ))
        gibbs.env$ll_history = append(gibbs.env$ll_history, list(ll_current))
        gibbs.env$eta_history = append(gibbs.env$eta_history, list(gibbs.env$eta_current))
      }
    }

    # Extracting the posterior for eta
    posterior_eta = do.call(rbind,gibbs.env$eta_history[(n_burnin/thinning_rate):length(gibbs.env$eta_history)])

    # Defining functions for extracting Pi and Delta from eta_hat
    eta_hat_to_Pi_hat = function(eta_hat, draw){

      # Extracting Pi from eta hat
      Pi_hat = model_to_Pi(eta_hat[1:(split_eta - 1)],J,K)

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


    # Posterior for Pi
    posterior_Pi = do.call(rbind,lapply(1:nrow(posterior_eta), function(d) eta_hat_to_Pi_hat(posterior_eta[d,1:(split_eta - 1)],d)))

    # Posterior for Delta
    posterior_Delta = do.call(rbind,lapply(1:nrow(posterior_eta), function(d) eta_hat_to_Pi_hat(posterior_eta[d,1:(split_eta - 1)],d)))

    #------------------------------------------------------------
    # Returning estimates and other info
    #------------------------------------------------------------

    # Return results
    return(list(
      posterior_Pi = posterior_Pi,
      posterior_Delta = posterior_Delta,
      posterior_eta = posterior_eta,
      ll_history = unlist(gibbs.env$ll_history),
      accepted_proposals = gibbs.env$accepted_proposals
    ))

  }


  #------------------------------------------------------------
  # Estimating Pi and Delta
  #------------------------------------------------------------

  # Is there more than one covariate value?
  if(class(tab)=="list"){

    # Determining the weight of each covariate cell
    W_weights = sapply(tab, function(tab_) sum(tab_$n)) |> unname()
    W_weights = W_weights / sum(W_weights)

    # Setting up parallel processing
    workers = parallel::makeCluster(cores)

    # Teaching the workers MisclassMLE
    parallel::clusterEvalQ(workers, require(misclassifyr)) |> invisible()
    parallel::clusterExport(workers, "MisclassMLE") |> invisible()

    # Drawing from the posterior of Pi and Delta within covariate cells
    MisclassMLE_out = pbapply::pblapply(
      MisclassMLE_inputs,
      MisclassMLE,
      cl = workers
    )

    # Restructuring MisclassMLE_out into lists of each output
    posterior_Pi = lapply(MisclassMLE_out, "[[", 1)
    posterior_Delta = lapply(MisclassMLE_out, "[[", 2)
    posterior_eta = lapply(MisclassMLE_out, "[[", 3)
    ll_history = lapply(MisclassMLE_out, "[[", 4)
    accepted_proposals = lapply(MisclassMLE_out, "[[", 5)

    MisclassMLE_out = list(
      posterior_Pi,
      posterior_Delta,
      posterior_eta,
      ll_history,
      accepted_proposals
    )

  } else {

    # Weight is one
    W_weights = 1

    # Estimating Pi and Delta for the full population
    MisclassMLE_out = MisclassMLE(MisclassMLE_inputs)

  }

  #------------------------------------------------------------
  # Generating figures
  #------------------------------------------------------------



  #------------------------------------------------------------
  # Estimating beta and SE
  #------------------------------------------------------------

  if(estimate_beta){

    # Estimating beta given Pi
    beta = Pi_to_beta(MisclassMLE_out$Pi_hat,X_vals, Y_vals, W_weights)

    # Estimating the SE via the delta method


  } else {

    return(MisclassMLE_out)

  }

}
