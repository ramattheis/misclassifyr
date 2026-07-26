#' Maps model parameters, psi, to Delta, the fully non-parametric  distribution of Y1, Y2 | Y*
#'
#' @param psi A numeric vector of length`J`^2*(`J`-1) containing all but the last row of `Delta`.
#' @return A numeric vector of length `J`^3 corresponding to the values of the
#'   `J`x`J`^2 matrix `Delta`. The rows of the matrix index the latent outcome
#'   `Y*`; the columns index the pair `(Y1, Y2)` with `Y1` varying fastest, so
#'   that column `(y2 - 1) * J + y1` holds `Pr(Y1 = y1, Y2 = y2 | Y* = y*)`.
#' @seealso [model_to_Delta_NP_ind()] for the conditionally independent
#'   version, [model_to_Delta_RL_ind()] for the record-linkage version.
#' @examples
#' # J = 2: 2 x 4 matrix, 2^2 * (2 - 1) = 4 free parameters
#' D <- matrix(model_to_Delta_NP(rep(0, 4)), nrow = 2)
#' dim(D)          # 2 x 4
#' rowSums(D)      # each row is a distribution over (Y1, Y2): all ones
#' @export
model_to_Delta_NP = function(psi){

  # J is a deterministic function of psi for model_to_Delta_NP
  # Using Cardano's formula to solve for J...
  cardano_discriminant = (length(psi)/2)^2  - 1/27
  cardano_u1 = length(psi)/2 + sqrt(cardano_discriminant)
  cardano_u2 = length(psi)/2 - sqrt(cardano_discriminant)
  J = as.integer(round(cardano_u1^(1/3) + cardano_u2^(1/3),0)) # rounding because of floating point weirdness

  # Building Delta Matrix
  Delta = matrix(psi, nrow = J)
  Delta = cbind(Delta, rep(0,J)) # Adding back the reference value
  Delta = t(apply(Delta,1,function(d) exp(d)/sum(exp(d)))) # logit link

  return(c(Delta))

}

# Adding a name as an attribute
attr(model_to_Delta_NP, "name") = "model_to_Delta_NP"

