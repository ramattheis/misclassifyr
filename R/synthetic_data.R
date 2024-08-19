#' Generates synthetic misclassification data.
#'
#' Longer description of what it does...
#'
#' @param J An integer indicating the dimension of Y.
#' @param K An integer indicating the dimension of X.
#' @param I An integer indicating the dimension of W.
#' @param sample_size An integer denoting the number of synthetic observations.
#' @param dgp_delta A character string indicating the data generating process for the synthetic noise
#' @param dgp_pi A character string indicating the data generating process for the joint distribution of X and Y*
#' @return A list including tabulated data `tab` and matrices `Pi`, `Delta`
#' @details some details...
#' @examples
#' \dontrun{
#' some example code # Should return something
#' }
#' @export
synthetic_data = function(J=5,
                          K=5,
                          I=2,
                          sample_size = 1e6,
                          dgp_delta = "Nonparametric, independent, strong diagonal",
                          dgp_pi = "Exponential"){


  #------------------------------------------------------------
  # (Inner loop) Defining a function to draw Pi, Delta, and
  # tab conditional on the control, W
  #------------------------------------------------------------

  within_control = function(){

    #------------------------------------------------------------
    # Drawing Pi
    #------------------------------------------------------------

    if(dgp_pi == "Exponential"){
      pi = do.call(rbind,lapply(1:J, function(j) exp(-1/2*(1:K - j)^2)))
      pi = pi/sum(pi)
    }

    #------------------------------------------------------------
    # Drawing Delta
    #------------------------------------------------------------

    if(dgp_delta == "Nonparametric, independent, strong diagonal"){

      delta1 = apply(diag(1/2,J) + matrix(exp(rnorm(J^2, mean = -5, sd = 1) ), J,J),
                     2, function(v) v/sum(v))
      delta2 = apply(diag(1/2,J) + matrix(exp(rnorm(J^2, mean = -5, sd = 1) ), J,J),
                     2, function(v) v/sum(v))

      delta = lapply(1:ncol(delta2), function(j) diag(delta2[j,]) %*% delta1)
      delta = do.call(cbind, delta)

    }

    if(dgp_delta == "Nonparametric, independent, medium diagonal"){

      delta1 = apply(diag(1/2,J) + matrix(exp(rnorm(J^2, mean = -4, sd = 1) ), J,J),
                     2, function(v) v/sum(v))
      delta2 = apply(diag(1/2,J) + matrix(exp(rnorm(J^2, mean = -4, sd = 1) ), J,J),
                     2, function(v) v/sum(v))

      delta = lapply(1:ncol(delta2), function(j) diag(delta2[j,]) %*% delta1)
      delta = do.call(cbind, delta)

    }

    if(dgp_delta == "Nonparametric, independent, weak diagonal"){

      delta1 = apply(diag(1/2,J) + matrix(exp(rnorm(J^2, mean = -3, sd = 1) ), J,J),
                     2, function(v) v/sum(v))
      delta2 = apply(diag(1/2,J) + matrix(exp(rnorm(J^2, mean = -3, sd = 1) ), J,J),
                     2, function(v) v/sum(v))

      delta = lapply(1:ncol(delta2), function(j) diag(delta2[j,]) %*% delta1)
      delta = do.call(cbind, delta)

    }


    if(dgp_delta == "Record Linkage, independent, 10 - 30%"){

      # Building Delta^{(1)} and Delta^{(2)}
      col1 = runif(J,0.1,0.3)
      col2 = runif(J,0.1,0.3)
      row1 = rexp(J)
      row1 = row1/sum(row1)
      row2 = rexp(J)
      row2 = row2/sum(row2)

      delta1 = diag(J)*(1-col1) + outer(row1,col1,"*")
      delta2 = diag(J)*(1-col2) + outer(row2,col2,"*")

      # Building the joint distribution Y1, Y2 conditional on Y*
      delta = lapply(1:ncol(delta2), function(j) diag(delta2[j,]) %*% t(delta1))
      delta = do.call(cbind, delta)

    }

    #------------------------------------------------------------
    # Drawing the data conditional on Pi and Delta
    #------------------------------------------------------------

    # Generate some tabulated data
    likelihood = t(pi) %*% delta
    X = rep(1:K,J^2)
    Y1 = rep(sapply(1:J, function(i) rep(i,K)),J)
    Y2 = c(sapply(1:J, function(i) rep(i,J*K)))
    tab = cbind(X,Y1,Y2) |> as.data.frame()
    tab$n = rmultinom(1, sample_size, prob = c(likelihood))

    #------------------------------------------------------------
    # Returning
    #------------------------------------------------------------

    return(list(
      tab = tab,
      Pi = pi,
      Delta = delta
    ))

  }

  #------------------------------------------------------------
  # (Outer loop) Drawing Pi, Delta, tab for each control cell
  #------------------------------------------------------------

  out_Pi = list()
  out_Delta = list()
  out_tab = list()

  for(i in 1:I){
    out_cell = within_control() # As opposed to an incel
    out_Pi = append(out_Pi, list(out_cell$Pi))
    out_Delta = append(out_Delta, list(out_cell$Delta))
    out_tab = append(out_tab, list(out_cell$tab))
  }

  #------------------------------------------------------------
  # Returning the list of objects
  #------------------------------------------------------------

  return(list(
    Pi = out_Pi,
    Delta = out_Delta,
    tab = out_tab
  ))

}
