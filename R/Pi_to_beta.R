#' Estimate a linear regression coefficient \eqn{\beta} from the joint distribution Pi
#'
#' This function computes \eqn{\beta} in the linear regression \eqn{Y^* = \alpha + \beta X + \epsilon} given the joint distribution of \eqn{X} and \eqn{Y^*}, \eqn{\Pi}, and a set of scalars associated with each value of \eqn{X} and \eqn{Y^*}.
#'
#' @param X_vals A numeric vector or list of numeric vectors providing the values of X associated with the columns of Pi.
#' @param Y_vals A numeric vector or list of numeric vectors providing the values of Y associated with the rows of Pi.
#' @param W_weights A numeric vector representing the sample size of each control cell.
#' @param W_names A character vector corresponding to the values of the control \eqn{W} in each cell.
#' @param mle Logical value indicating whether maximum likelihood estimates of \eqn{\Pi} have been provided. Defaults to `TRUE`.
#' @param bayesian Logical value indicating whether posterior draws of \eqn{\Pi} have been provided. Defaults to `FALSE`.
#' @param Pi_mle A numeric vector or list of numeric vectors containing the elements of Pi.
#' @param cov_Pi A numeric vector or a list of numeric vectors representing the covariance of estimates of the elements of Pi.
#' @param posterior_Pi A data frame of posterior draws of Pi (or a list of such data frames split by control cell), as returned in `$posterior_Pi` by `misclassifyr(bayesian = TRUE)`.
#' @param ll_history A data frame with columns `ll` and `draw` (or a list of such data frames split by control cell), as returned in `$ll_history` by `misclassifyr(bayesian = TRUE)`. Required for the Chen-Christensen-Tamer HPD confidence set; if omitted, `HPD_draws` and `HPDCI` are returned as `NA`.
#' @export
Pi_to_beta = function(
    X_vals,
    Y_vals,
    W_weights = NA,
    W_names = NA,
    mle = T,
    bayesian = F,
    Pi_mle = NA,
    cov_Pi = NA,
    posterior_Pi = NA,
    ll_history = NA
    ){

  #------------------------------------------------------------
  # Catching input errors
  #------------------------------------------------------------

  if(mle){
    # Throwing an error if MLE esimates and variances of Pi aren't provided
    if((identical(Pi_mle, NA) | identical(cov_Pi,NA))){
      stop("If `mle == TRUE`, `Pi_mle` and `cov_Pi` should be provided.")
    }
    # Throwing an error if not all objects are lists / all objects are not lists
    if(bayesian && (is_cell_list(Pi_mle) || is_cell_list(posterior_Pi)) &&
       (!(is_cell_list(Pi_mle) && is_cell_list(posterior_Pi)) ||
        length(Pi_mle) != length(posterior_Pi))){
      stop("If `Pi_mle` or `posterior_Pi` is a list, both should be lists of the same length.  ")
    }
    # Throwing an error if W_weights is not provided
    if(is_cell_list(Pi_mle) & identical(W_weights, NA)){
      stop("If `Pi_mle` is a list, `W_weights` should be provided.")
    }
  }

  if(bayesian){
    if(is.null(ll_history)){ ll_history = NA }
    if(identical(posterior_Pi, NA)){
      stop("If `bayesian == TRUE`, `posterior_Pi` should be provided.")
    }
    # Throwing an error if W_weights is not provided
    if(is_cell_list(posterior_Pi) & identical(W_weights, NA)){
      stop("If `posterior_Pi` is a list, `W_weights` should be provided.")
    }
    # Throwing an error if ll_history's structure doesn't match posterior_Pi
    if(!identical(ll_history, NA)){
      if(is_cell_list(posterior_Pi) != is_cell_list(ll_history)){
        stop("`ll_history` should have the same structure as `posterior_Pi` (both single data frames, or both lists split by control cell).")
      }
      if(is_cell_list(ll_history) && length(ll_history) != length(posterior_Pi)){
        stop("`ll_history` and `posterior_Pi` should be lists of the same length.")
      }
    }
  }

  if(identical(W_names, NA) &  !identical(W_weights, NA)){
    # Naming control cells if no names are provided
    W_names = paste0("control cell ",seq_along(W_weights))
  }

  #------------------------------------------------------------
  # Computing beta and the SE for MLE
  #------------------------------------------------------------

  # Computing beta from Pi for MLE
  if(mle){

    # Extracting W_weights and normalizing
    W_weights = unlist(W_weights)
    W_weights = W_weights/sum(W_weights)

    # Estimating beta given Pi
    beta_hat_mle = Pi_to_beta_inner(Pi_mle, X_vals, Y_vals, W_weights)
    names(beta_hat_mle) = "MLE beta"

    # Estimating the SE via the delta method
    se_beta_mle = se_beta_deltamethod(
      Pi = Pi_mle,
      cov_Pi = cov_Pi,
      X_vals,
      Y_vals,
      W_weights
    )
    names(se_beta_mle) = "SE for MLE beta"

    # If controls are used, returning a list of betas within each cell
    if(is_cell_list(Pi_mle)){

      # Estimating beta within each control cell
      betas_hat_mle = sapply(seq_along(W_names), function(j)
        Pi_to_beta_inner(Pi_mle[[j]], X_vals[[j]], Y_vals[[j]], 1))
      names(betas_hat_mle) = paste("MLE beta in cell", W_names)

      # Estimating the SE in each control cell
      se_betas_mle = sapply(seq_along(W_names), function(j)
        se_beta_deltamethod(Pi_mle[[j]], cov_Pi[[j]], X_vals[[j]], Y_vals[[j]], 1))
      names(se_betas_mle) = paste("SE for MLE beta in cell", W_names)

    } else {
      # If there are no controls, returning NA
      betas_hat_mle = NA
      se_betas_mle = NA
    }

  } else {
    # If MLE is not computed, returning NAs
    beta_hat_mle = NA
    se_beta_mle = NA
    betas_hat_mle = NA
    se_betas_mle = NA
  }


  #------------------------------------------------------------
  # Computing the posterior of beta given posterior draws of Pi
  #------------------------------------------------------------

  if(bayesian){

    # Aggregating across control cells
    if(is_cell_list(posterior_Pi)){

      # Extracting W_weights and normalizing
      W_weights = unlist(W_weights)
      W_weights = W_weights/sum(W_weights)

      # The posterior draw index (common across control cells)
      draws = unique(posterior_Pi[[1]]$draw)

      # Defining a quick function to aggregate the posterior across control cells
      posterior_agg = function(d){
        # Grabbing the dth draw of the posterior
        posterior_df_list = lapply(seq_along(W_names), function(w)
          subset(posterior_Pi[[w]], draw == d))

        # Scaling weights to reflect control cell size
        posterior_df_list = Map(function(post_df, W_weight) {
          post_df$Pi_hat = post_df$Pi_hat * W_weight
          return(post_df)
        }, posterior_df_list, W_weights)

        # Binding together and returning
        posterior_df = do.call(rbind, posterior_df_list)

        return(posterior_df)
      }

      posterior_beta = sapply(draws, function(d){
        posterior_df = posterior_agg(d)
        lm(Y_val ~ X_val,
           data = posterior_df,
           weight = posterior_df$Pi_hat
        )$coefficients[2] |> unname()
      })

      # Recording the posterior of beta within covariate cells
      posterior_betas = lapply(seq_along(W_names), function(j)
        sapply(unique(posterior_Pi[[j]]$draw), function(d){
          posterior_df = subset(posterior_Pi[[j]], draw == d)
          lm(Y_val ~ X_val,
             data = posterior_df,
             weight = posterior_df$Pi_hat
          )$coefficients[2] |> unname()
        })
      )
      names(posterior_betas) = W_names

      # Computing Chen, Christensen, and Tamer partial-ID-robust CI
      # (requires the MCMC likelihood history)
      if(!identical(ll_history, NA)){

        # Summing the likelihood history across covariate cells
        ll_draws = sapply(draws, function(d) {
          sum(sapply(seq_along(W_names), function(w) {
            subset(ll_history[[w]], draw == d)$ll
          }))
        })

        # sort to find top 95% (with fewer than 20 draws the 5% cutoff index
        # is 0, so clamp to the minimum: all draws are kept)
        ll_sorted = sort(ll_draws, decreasing = F)
        ll_critical = ll_sorted[max(1, floor(length(ll_sorted)*0.05))]
        HPD_draws = draws[ll_draws >= ll_critical]
        HPDCI = c(min(posterior_beta[ll_draws >= ll_critical]),
                  max(posterior_beta[ll_draws >= ll_critical]))
      } else {
        HPD_draws = NA
        HPDCI = NA
      }

      # Extracting the median and the sd
      posterior_beta_med = median(posterior_beta)
      posterior_beta_sd = sd(posterior_beta)
      posterior_betas_med = lapply(posterior_betas, median)
      posterior_betas_sd = lapply(posterior_betas, sd)

    } else {

      # The posterior draw index
      draws = unique(posterior_Pi$draw)

      # Estimating beta for each draw of the posterior
      posterior_beta = sapply(draws, function(d){
        posterior_df = subset(posterior_Pi, draw == d)
        lm(Y_val ~ X_val,
           data = posterior_df,
           weight = posterior_df$Pi_hat
        )$coefficients[2] |> unname()
      })

      # Computing Chen, Christensen, and Tamer partial-ID-robust CI
      # (requires the MCMC likelihood history)
      if(!identical(ll_history, NA)){

        # Summing the likelihood over each draw
        ll_draws = sapply(draws, function(d) {
          sum(subset(ll_history, draw == d)$ll)
        })

        # sort to find top 95% (with fewer than 20 draws the 5% cutoff index
        # is 0, so clamp to the minimum: all draws are kept)
        ll_sorted = sort(ll_draws, decreasing = F)
        ll_critical = ll_sorted[max(1, floor(length(ll_sorted)*0.05))]
        HPD_draws = draws[ll_draws >= ll_critical]
        HPDCI = c(min(posterior_beta[ll_draws >= ll_critical]),
                  max(posterior_beta[ll_draws >= ll_critical]))
      } else {
        HPD_draws = NA
        HPDCI = NA
      }

      # Recording the median and SD
      posterior_beta_med = median(posterior_beta)
      posterior_beta_sd = sd(posterior_beta)

      # Returning NA if there are not multiple control cells
      posterior_betas = NA
      posterior_betas_med = NA
      posterior_betas_sd = NA

    }


  } else {

    # Returning NAs if the posterior of Pi was not provided
    posterior_beta = NA
    posterior_beta_med = NA
    posterior_beta_sd = NA
    posterior_betas = NA
    posterior_betas_med = NA
    posterior_betas_sd = NA
    HPD_draws = NA
    HPDCI = NA

  }

  #------------------------------------------------------------
  # Binding and returning results
  #------------------------------------------------------------

  return(list(
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
    HPD_draws = HPD_draws,
    HPDCI = HPDCI
  ))
}
