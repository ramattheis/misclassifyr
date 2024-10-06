#' Creates a function for model_to_Delta based on RL structure errors and the empirical distribution of Y1 and Y2
#'
#' @param tab tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should have names `Y1`, `Y2`, `X`, and `n` where `n` is a non-negative numeric vector corresponding to the counts of `Y1`,`Y2`, and `X`. The rows should be ordered according to `order(Y2,Y1,X)`.
#' @param J An integer or list corresponding to the number of unique values of `Y1` and `Y2`.
#' @return A list including 1. a function or a list of functions that map model parameters `psi` to the misclassification matrix `Delta`, 2. a vector or list of vectors corresponding to initial values of psi , and 3. a function or list of functions for the log prior of Delta in this model.
#' @export
make_empirical_Delta_RL = function(tab,J){

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

    FY2s = lapply(tab, function(tb) {
      tabY2 = tb |>
        dplyr::arrange(Y2) |>
        dplyr::group_by(Y2) |>
        dplyr::summarise( n = sum(n)) |>
        as.data.frame()
      tabY2$n / sum(tabY2$n)
    })

    # Pasting together the empirical frequency of Y1 and Y2
    FY12s = Map(c,FY1s,FY2s)

    # Defining a function for model_to_Delta for each covariate cell based on the
    # empirical frequency of Y1 and Y2
    model_to_Delta = lapply(FY12s, function(FY12) {

      function(psi) {

        # Transforming phi to return to probabilities (no sum-to-one constraint)
        alpha = exp(psi)/(1+exp(psi))

        # Splitting phi into components associated with Y1 and Y2
        alpha1 = alpha[1:(length(alpha)/2)]
        alpha2 = alpha[(length(alpha)/2 +1):length(alpha)]

        # Recording the marginal distribution of Y1 and Y2
        FY1 = FY12[1:(length(FY12)/2)]
        FY2 = FY12[(length(FY12)/2 + 1):length(FY12)]

        # Computing the misclassification error distribution
        Delta1 = diag(1-alpha1) + FY1 %*% t(alpha1)
        Delta2 = diag(1-alpha2) + FY2 %*% t(alpha2)
        Delta = lapply(1:nrow(Delta2), function(j) diag(Delta2[j,]) %*% t(Delta1))
        Delta = do.call(cbind, Delta)

        return(c(Delta))
      }
    })

  } else {

    # Computing the empirical frequency of Y1 and Y2 within covariate cells
    tabY1 = tab |>
      dplyr::arrange(Y1) |>
      dplyr::group_by(Y1) |>
      dplyr::summarise( n = sum(n)) |>
      as.data.frame()
    FY1 = tabY1$n / sum(tabY1$n)

    tabY2 = tab |>
      dplyr::arrange(Y2) |>
      dplyr::group_by(Y2) |>
      dplyr::summarise( n = sum(n)) |>
      as.data.frame()
    tabY2$n / sum(tabY2$n)
    FY2 = tabY2$n / sum(tabY2$n)

    model_to_Delta <- local({

      function(psi) {

        # Transforming phi to return to probabilities (no sum-to-one constraint)
        alpha = exp(psi)/(1+exp(psi))

        # Splitting phi into components associated with Y1 and Y2
        alpha1 = alpha[1:(length(alpha)/2)]
        alpha2 = alpha[(length(alpha)/2 +1):length(alpha)]

        # Computing the misclassification error distribution
        Delta1 = diag(1-alpha1) + FY1 %*% t(alpha1)
        Delta2 = diag(1-alpha2) + FY2 %*% t(alpha2)
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
    psi_0 = lapply(seq_along(tab), function(i) rep(-2,2*J))
  } else {
    psi_0 = rep(-2,2*J)
  }

  #------------------------------------------------------------
  # Defining the prior
  #------------------------------------------------------------

  if(class(tab) == "list"){
    log_prior_Delta = replicate(length(tab),dlogis)
  } else {
    log_prior_Delta = dlogis
  }

  #------------------------------------------------------------
  # Returning a list with each object
  #------------------------------------------------------------

  # Returning a list
  return(list(
    model_to_Delta,
    psi_0,
    log_prior_Delta
  ))

}

