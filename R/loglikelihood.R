#' Returns the log likelihood of the data.
#'
#' Longer description of what it does...
#' It's assumed that`tab`, `J`, and `K` are defined in the environment.
#'
#' @param theta A numeric vector of length  Jx(J^2+K) describing the joint distribution of the data.
#' @return the log likelihood of the data given theta, i.e. Pi and Delta.
#' @details some details.
#' @examples
#' \dontrun{
#' some example code # Should return something
#' }
#' @export
loglikelihood = function(theta){

  # Importing tab, J, K, and lambda from the shared environment
  misclassifyr_env = get(".misclassifyr_env", envir = asNamespace("misclassifyr"))
  tab = misclassifyr_env$tab
  J = misclassifyr_env$J
  K = misclassifyr_env$K
  lambda_pos = misclassifyr_env$lambda_pos
  lambda_dd = misclassifyr_env$lambda_dd

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

