#' misclassifyr
#'
#' This function provides a menu of options for estimation and inference of misclassification models in which the analyst has access to two noisy measures, `Y1` and `Y2` of a latent outcome `Y*`, a correctly measured covariate `X`, and discrete controls `W`.
#'
#' @import parallel numDeriv pracma
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
#' @param optim_maxit An integer for the maximum number of iterations in numerical optimization, passed to `optim()`
#' @param optim_tol A positive number defining convergence in numerical optimization, passed to `optim()`
#' @param check_stability A logical value indicating whether to perform a more rigorous stability test for the numerical optimizer.
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
    optim_maxit = 1e5,
    optim_tol = 1e-8,
    check_stability = F,
    stability_sd = 0.1,
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
                  class(Y_vals),
                  class(optim_maxit),
                  class(optim_tol),
                  class(check_stability),
                  class(stability_sd)) |>
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
    "Y_vals",
    "optim_maxit",
    "optim_tol",
    "check_stability",
    "stability_sd")

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
      Y_vals = Y_vals[[j]],
      optim_maxit = optim_maxit[[j]],
      optim_tol = optim_tol[[j]],
      check_stability = check_stability[[j]],
      stability_sd = stability_sd[[j]]
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
      Y_vals = Y_vals,
      optim_maxit = optim_maxit,
      optim_tol = optim_tol,
      check_stability = check_stability,
      stability_sd = stability_sd
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
    lambda_pos = MisclassMLE_input$lambda_pos
    lambda_dd = MisclassMLE_input$lambda_dd
    optim_maxit = MisclassMLE_input$optim_maxit
    optim_tol = MisclassMLE_input$optim_tol
    check_stability = MisclassMLE_input$check_stability
    stability_sd = MisclassMLE_input$stability_sd

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
    # Setting the starting location for optimization
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
        phi_0 = log(tab_xy$p[1:(J*K-1)])
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

    # Creating an environment to keep track of penalty scaling parameters
    # and other parameters across draws
    optim.env = new.env()
    optim.env$tab = tab
    optim.env$J = J
    optim.env$K = K
    optim.env$lambda_pos = sum(tab$n)^2
    optim.env$lambda_dd = sum(tab$n)^2

    # Defining the objective function
    objective = function(eta){

      # Mapping model arguments to the entries of the joint distribution
      Pi_ = model_to_Pi(eta[1:(split_eta-1)], J, K)
      Delta_ = model_to_Delta(eta[(split_eta):length(eta)])

      # Computing the log likelihood of the data given Pi and Delta
      ll =  loglikelihood(
        tab = optim.env$tab,
        theta = c(c(Pi_), c(Delta_)),
        J = optim.env$J,
        K = optim.env$K,
        lambda_pos = optim.env$lambda_pos,
        lambda_dd = optim.env$lambda_dd
        )

      # Returning -1x the log likelihood + penalties
      return(-1*ll)

    }

    #------------------------------------------------------------
    # Optimizing objective w.r.t. eta
    #------------------------------------------------------------

    # Optimizing objective w.r.t. eta
    out = optim(par = eta_0,
                fn = objective,
                method = "BFGS",
                hessian = T,
                control = list(maxit = optim_maxit,
                               reltol = optim_tol))

    # Throwing an error if optim did not converge successfully
    if(out$convergence == 1){stop("Error: Maximum iterations reached before numerical convergence, consider increasing optim max iterations or tolerance.")}
    if(out$convergence > 1){stop(paste("Error: Numerical optimization failed with message:",
                                        out$message,
                                        "... consider alternative `optim` settings."))}

    # Raising a warning if estimates are close to the initial value
    if(sum((exp(out$par) - exp(eta_0))^2) < 0.1){warning("Solution appears close to the intial value. Consider choosing a more conservative tolerance for optim and/or proceed with caution.")}

    #------------------------------------------------------------
    # Testing the stability of the optimizer
    #------------------------------------------------------------

    # If check_stability, increase the number of alternative starting points
    if(check_stability){
      extra_starting_points = 9
    } else {
      extra_starting_points = 1
    }

    # Recording initial eta hat
    eta_hat1 = out$par

    # Initializing total inconsistency
    inconsistency = 0

    for(draw in 1:extra_starting_points){

      # Additional estimates on perturbed starting locations
      out2 = optim(par = eta_0 + rnorm(length(eta_0),sd = stability_sd),
                   fn = objective,
                   method = "BFGS",
                   control = list(maxit = optim_maxit,
                                  reltol = optim_tol))
      eta_hat2 = out2$par

      # Adding the (absolute) difference between e^eta_hat1 and e^eta_hat2
      inconsistency = inconsistency + sum(abs(exp(eta_hat1) - exp(eta_hat2)))
    }

    if(inconsistency > 0.01){
      warning("Optimal eta is inconsitent across starting locations.")
    }

    #------------------------------------------------------------
    # Computing the covariance matrix of Pi
    #------------------------------------------------------------

    # Is the hessian of the fisher information invertible?
    # Note that the objective function is flipped, so we don't need to multiply by -1 again
    fisher_info = out$hessian[1:(split_eta-1),1:(split_eta-1)]
    if (abs(det(fisher_info)) < 1e-9) {
      fisher_info_err = "Fisher information matrix is not invertible."
      warning("The Fisher information matrix is not invertible; using the Moore-Penrose inverse instead -- analytical variances may be inaccurate.")
    } else {
      fisher_info_err = "Fisher information matrix is invertible."
    }

    # Computing the Jacobian of model_to_Pi
    model_to_Pi_wrapper = function(phi){
      return(model_to_Pi(phi,J,K))
    }

    model_to_Pi_jacobian = numDeriv::jacobian(model_to_Pi_wrapper, out$par[1:(split_eta-1)])

    # Computing the covariance matrix of Pi
    cov_Pi =  model_to_Pi_jacobian %*% pracma::pinv(fisher_info) %*% t(model_to_Pi_jacobian)

    # Forcing the diagonal of cov_Pi to be non-negative
    diag(cov_Pi) = pmax(diag(cov_Pi),0)

    #------------------------------------------------------------
    # Returning estimates and other info
    #------------------------------------------------------------

    # Return results
    return(list(
      Pi_hat = model_to_Pi_wrapper(out$par[1:(split_eta-1)]),
      Delta_hat = model_to_Delta(out$par[(split_eta):(length(out$par))]),
      cov_Pi = cov_Pi,
      eta_hat = out$par,
      log_likelihood = -1*out$value,
      optim_counts = out$counts,
      model_to_Pi_jacobian = model_to_Pi_jacobian,
      eta_hessian = out$hessian,
      fisher_info_err = fisher_info_err,
      inconsistency = inconsistency
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

    # Estimating Pi and Delta within covariate cells
    MisclassMLE_out = pbapply::pblapply(
      MisclassMLE_inputs,
      MisclassMLE,
      cl = workers
    )

    # Restructuring MisclassMLE_out into lists of each output
    Pi_hat = lapply(MisclassMLE_out, "[[", 1)
    Delta_hat = lapply(MisclassMLE_out, "[[", 2)
    cov_Pi = lapply(MisclassMLE_out, "[[", 3)
    eta_hat = lapply(MisclassMLE_out, "[[", 4)
    log_likelihood = lapply(MisclassMLE_out, "[[", 5)
    optim_counts = lapply(MisclassMLE_out, "[[", 6)
    model_to_Pi_jacobian = lapply(MisclassMLE_out, "[[", 7)
    eta_hessian = lapply(MisclassMLE_out, "[[", 8)
    fisher_info_err = lapply(MisclassMLE_out, "[[", 9)
    inconsistency = lapply(MisclassMLE_out, "[[", 10)

    MisclassMLE_out = list(
      Pi_hat = Pi_hat,
      Delta_hat = Delta_hat,
      cov_Pi = cov_Pi,
      eta_hat = eta_hat,
      log_likelihood = log_likelihood,
      optim_counts = optim_counts,
      model_to_Pi_jacobian = model_to_Pi_jacobian,
      eta_hessian = eta_hessian,
      fisher_info_err = fisher_info_err,
      inconsistency = inconsistency
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
