#' Returns the log likelihood of the data.
#'
#' @param theta A numeric vector of length  Jx(J^2+K) describing the joint distribution of the data.
#' @param tab
#' @param J An integer or list corresponding to the number of unique values of `Y1` and `Y2`.
#' @param K An integer or list corresponding to the number of unique values of `X`.
#' @param lambda_pos A numeric value scaling a violations of positivity for the entries of Pi and Delta.
#' @param lambda_dd A numeric value scaling a violations of diagonal dominance for Delta.
#' @return the log likelihood of the data given theta, i.e. Pi and Delta.
#' @keywords internal
#' @export
loglikelihood = function(theta,tab,J,K,lambda_pos,lambda_dd){

  # Building Pi
  Pi = theta[1:(J*K)]       # Extracting all but the last entry of Pi
  Pi = matrix(Pi,nrow = J)  # Converting to matrix, Y* rows, X cols

  # Building Delta
  Delta = matrix(theta[(J*K+1):length(theta)], nrow = J)

  # Adding a penalty for positivity
  penalty_pos = -1*lambda_pos*sum(sapply(theta, function(z) min(z,0))^2)

  # Adding a penalty for diagonal dominance
  penalty_dd = -1*lambda_dd*sum(
    sapply(1:J^2, function(i)
      pmax(0,Delta[ifelse(i>=J & i%%J==0,J,i%%J),
                   (floor((i-1)/J)*J+1):(floor((i-1)/J+1)*J)] -
             Delta[ifelse(i>=J & i%%J==0,J,i%%J),i])^2))

  # Computing the log likelihood
  llsum = sum(tab$n * softlog(c(t(Pi) %*% Delta)))

  return(llsum + penalty_dd + penalty_pos)

}

