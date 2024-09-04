#' Returns the log likelihood of the data.
#'
#' @param theta A numeric vector of length  Jx(J^2+K) describing the joint distribution of the data.
#' @param tab A dataframe or a list of dataframes containing tabulated data or a list of tabulated data split by controls. The columns should be numeric with names `Y1`, `Y2`, `X`, and `n` where `Y1` and `Y2` take each value between 1 and `J`, `X` takes each value between `1` and `K`, and
#' @param J An integer or list corresponding to the number of unique values of `Y1` and `Y2`.
#' @param K An integer or list corresponding to the number of unique values of `X`.
#' @param lambda_pos A numeric value scaling a violations of positivity for the entries of Pi and Delta.
#' @param lambda_dd A numeric value scaling a violations of diagonal dominance for Delta.
#' @return the log likelihood of the data given theta, i.e. Pi and Delta.
#' @keywords internal
loglikelihood = function(theta,tab,J,K,lambda_pos,lambda_dd){

  # Building Pi
  Pi = theta[1:(J*K)]       # Extracting all but the last entry of Pi
  Pi = matrix(Pi,nrow = J)  # Converting to matrix, Y* rows, X cols

  # Building Delta
  Delta = matrix(theta[(J*K+1):length(theta)], nrow = J)

  # Adding a penalty for diagonal dominance
  penalty_dd = -1*lambda_dd*sum(
    sapply(1:J^2, function(i)
      pmax(0,Delta[ifelse(i>=J & i%%J==0,J,i%%J),
                   (floor((i-1)/J)*J+1):(floor((i-1)/J+1)*J)] -
             Delta[ifelse(i>=J & i%%J==0,J,i%%J),i])^2))

  # Computing the log likelihood
  llsum = sum(tab$n * softlog(c(t(Pi) %*% Delta)))

  return(llsum + penalty_dd)

}

