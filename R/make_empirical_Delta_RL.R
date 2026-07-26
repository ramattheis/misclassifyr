#' Creates a function for model_to_Delta based on RL structure errors and the empirical distribution of Y1 and Y2
#'
#' @param tab tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should have names `Y1`, `Y2`, `X`, and `n` where `n` is a non-negative numeric vector corresponding to the counts of `Y1`,`Y2`, and `X`. The rows should be ordered according to `order(Y2,Y1,X)`.
#' @param J An integer or list corresponding to the number of unique values of `Y1` and `Y2`.
#' @return A list with three elements ready to be handed to
#'   [misclassifyr()]: `model_to_Delta`, a function (or list of functions,
#'   one per control cell) mapping model parameters `psi` to the
#'   misclassification matrix `Delta`; `psi_0`, a vector (or list of vectors)
#'   of starting values; and `log_prior_Delta`, a function (or list of
#'   functions) giving the log prior of `Delta` in this model.
#' @examples
#' set.seed(1)
#' syn <- synthetic_data(J = 3, K = 3, I = 1, sample_size = 5000,
#'                       dgp_delta = "Record Linkage, independent, 10 - 30%")
#' design <- make_empirical_Delta_RL(syn$tab[[1]], J = 3)
#'
#' # 2 * J false-link rates (one per category, per measure); psi_0 starts them
#' # all near plogis(-2) = 0.12
#' length(design$psi_0)
#'
#' # Verify the design: J x J^2 with rows summing to one
#' D <- matrix(design$model_to_Delta(design$psi_0), nrow = 3)
#' dim(D)
#' rowSums(D)
#' design$log_prior_Delta(design$psi_0)
#' @export
make_empirical_Delta_RL = function(tab,J){

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  # If (and only if) tab is a list, J should be a list of the same length
  if(is_cell_list(tab)){
    if(is_cell_list(J)){
      if(length(J) != length(tab)){
        stop("If `tab` is a list, the list `J` should have the same length.")
      }
    } else {
      stop("If `tab` is a list, `J` should also be a list.")
    }
  } else {
    if(is_cell_list(J)){ stop("`J` should not be a list if `tab` is not a list.") }
  }


  #------------------------------------------------------------
  # Constructing functions for model_to_Delta for each cell
  #------------------------------------------------------------

  if(is_cell_list(tab)){

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

        # Transforming psi to return to probabilities (no sum-to-one constraint)
        alpha = exp(psi)/(1+exp(psi))

        # Splitting psi into components associated with Y1 and Y2
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

        # Transforming psi to return to probabilities (no sum-to-one constraint)
        alpha = exp(psi)/(1+exp(psi))

        # Splitting psi into components associated with Y1 and Y2
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

  if(is_cell_list(tab)){
    psi_0 = lapply(seq_along(tab), function(i) rep(-2,2*J[[i]]))
  } else {
    psi_0 = rep(-2,2*J)
  }

  #------------------------------------------------------------
  # Defining the prior
  #------------------------------------------------------------

  if(is_cell_list(tab)){
    log_prior_Delta = replicate(length(tab), function(psi) sum(dlogis(psi, log = T)))
  } else {
    log_prior_Delta = function(psi) {return(sum(dlogis(psi, log = T)))}
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

