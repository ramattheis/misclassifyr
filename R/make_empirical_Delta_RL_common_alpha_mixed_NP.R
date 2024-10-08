#' Creates a function for model_to_Delta based on a combination of RL structure errors and the empirical distribution of Y1 assuming RL errors are independent of Y* and non-parametric errors in Y2
#'
#' @param tab tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should have names `Y1`, `Y2`, `X`, and `n` where `n` is a non-negative numeric vector corresponding to the counts of `Y1`,`Y2`, and `X`. The rows should be ordered according to `order(Y2,Y1,X)`.
#' @return A list including 1. a function or a list of functions that map model parameters `psi` to the misclassification matrix `Delta`, 2. a vector or list of vectors corresponding to initial values of psi , and 3. a function or list of functions for the log prior of Delta in this model.
#' @export
make_empirical_Delta_RL_common_alpha_mixed_NP = function(tab,J){

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  # If (and only if) tab is a list, J should be a list of the same length
  if(class(tab) == "list"){
    if(class(J) == "list"){
      if(length(J) != length(tab)){
        stop("If `tab` is a list, the list `J` should have the same length.")
      }
    } else {
      stop("If `tab` is a list, `J` should also be a list.")
    }
  } else {
    if(class(J) == "list"){ stop("`J` should not be a list if `tab` is not a list.") }
  }


  #------------------------------------------------------------
  # Constructing functions for model_to_Delta for each cell
  #------------------------------------------------------------

  if(class(tab) == "list"){

    # Computing the empirical frequency of Y1 and Y2 within covariate cells
    FY1s = lapply(tab, function(tb) {
      tabY1 = tb |>
        dplyr::arrange(Y1) |>
        dplyr::group_by(Y1) |>
        dplyr::summarise( n = sum(n)) |>
        as.data.frame()
      tabY1$n / sum(tabY1$n)
    })

    # Defining a function for model_to_Delta for each covariate cell based on the
    # empirical frequency of Y1 and non-parametric model of Y2
    model_to_Delta = lapply(FY1s, function(FY1) {

      function(psi) {

        # Extracting alpha (the RL error rate in Y1)
        alpha = exp(psi[1])/(1+exp(psi[1]))
        psi = psi[-1]

        # Assembling misclassification error matrix for Y1
        Delta1 = diag(rep(1-alpha,length(FY1))) + FY1 %*% t(rep(alpha,length(FY1)))

        # Building non-parametric Delta^{(2)}
        Delta2 = matrix(psi, nrow = length(FY1)-1)
        Delta2 = rbind(Delta2, rep(0, ncol(Delta2))) # Adding back the reference value
        Delta2 = apply(Delta2,2,function(d) exp(d)/sum(exp(d)))

        # Computing the misclassification error distribution
        Delta = lapply(1:nrow(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
        Delta = do.call(cbind, Delta)

        return(c(Delta))
      }
    })

  } else {

    # Computing the empirical frequency of Y1
    tabY1 = tab |>
      dplyr::arrange(Y1) |>
      dplyr::group_by(Y1) |>
      dplyr::summarise( n = sum(n)) |>
      as.data.frame()
    FY1 = tabY1$n / sum(tabY1$n)

    model_to_Delta <- local({

      function(psi) {

        # Extracting alpha (the RL error rate in Y1)
        alpha = exp(psi[1])/(1+exp(psi[1]))
        psi = psi[-1]

        # Assembling misclassification error matrix for Y1
        Delta1 = diag(rep(1-alpha,length(FY1))) + FY1 %*% t(rep(alpha,length(FY1)))

        # Building non-parametric Delta^{(2)}
        Delta2 = matrix(psi, nrow = length(FY1)-1)
        Delta2 = rbind(Delta2, rep(0, ncol(Delta2))) # Adding back the reference value
        Delta2 = apply(Delta2,2,function(d) exp(d)/sum(exp(d)))

        # Computing the misclassification error distribution
        Delta = lapply(1:nrow(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
        Delta = do.call(cbind, Delta)

        return(c(Delta))

      }

    })
  }

  #------------------------------------------------------------
  # Defining the initial value of psi_0
  #------------------------------------------------------------

  if(class(tab) == "list"){

    psi_0 = lapply(seq_along(tab), function(i) {

      # Initial value of Delta2 is RL errors with 20% linkage error
      Delta2 =  diag(J[[i]])*(1 - 0.2) + outer(rep(1 / J[[i]], J[[i]]), rep(0.2, J[[i]]), "*")

      # Transform to probabilities to logs
      Delta2 = apply(Delta2, 2, function(col) { softlog(col / col[length(col)]) })

      # Removing the reference column and flattening to a vector
      psi2 = c(Delta2[-J[[i]],])

      # Combine with initial value for the common-alpha RL error portion of the model
      return(c(-2,psi2))

    })

  } else {

    # Initial value of Delta2 is RL errors with 20% linkage error
    Delta2 =  diag(J)*(1 - 0.2) + outer(rep(1 / J, J), rep(0.2, J), "*")

    # Transform to probabilities to logs
    Delta2 = c(apply(Delta2, 2, function(col) { softlog(col / col[length(col)]) }))

    # Removing the reference column and flattening to a vector
    psi2 = c(Delta2[-J,])

    # Combine with initial value for the common-alpha RL error portion of the model
    psi_0 = c(-2,psi2)

  }

  #------------------------------------------------------------
  # Defining the prior
  #------------------------------------------------------------

  if(class(tab) == "list"){
    # The prior for the Delta2 portion is can be written as 1/2 the log prior for Delta1, Delta2 in the
    # non-parameteric model assuming independence, since the function is symmertic in the componets of psi
    # associated with Delta1/Delta2, and the log prior is the sum of these two components.

    log_prior_Delta = replicate(length(tab), function(psi) {

       # Prior for the common-alpha RL model for Y1
      alpha_component =  dlogis(psi[1], log = T)

      # Non-parametric misclassificaiton distribution for Y2
      Delta2_component = 1/2*log_prior_Delta_NP_ind(c(psi[-1],psi[-1]))

      # Log prior is the sum of each component
      return(alpha_component + Delta2_component)

    })
  } else {

    log_prior_Delta = function(psi) {

      # Prior for the common-alpha RL model for Y1
      alpha_component =  dlogis(psi[1], log = T)

      # Non-parametric misclassificaiton distribution for Y2
      Delta2_component = 1/2*log_prior_Delta_NP_ind(c(psi[-1],psi[-1]))

      # Log prior is the sum of each component
      return(alpha_component + Delta2_component)
    }
  }

  #------------------------------------------------------------
  # Returning a list with each object
  #------------------------------------------------------------

  # Returning a list
  return(list(
    model_to_Delta = model_to_Delta,
    psi_0 = psi_0,
    log_prior_Delta = log_prior_Delta
  ))

}
