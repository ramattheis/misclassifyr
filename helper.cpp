#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ll(NumericVector theta,
          DataFrame tab,
          double lambda_pos = 1e6,
          double lambda_dd = 1e4) {

  // Determining J
  int J = (2 + sqrt(4 + 12 * (theta.size() + 1))) / 6;

  // Extracting rho
  NumericVector rho = theta[Range(0, J - 2)];
  rho.push_back(1 - sum(rho)); // enforcing that parameters sum to one

  // Extracting Pi_X
  NumericMatrix Pi_x(J, J);
  for (int j = 0; j < J; ++j) {
    NumericVector pi_x = theta[Range((j+1)*(J-1), (j+2)*(J-1)-1)];
    pi_x.push_back(1 - sum(pi_x)); // enforcing parameters sum to one
    Pi_x(_, j) = pi_x;
  }

  // Extracting Delta_1 and Delta_2
  NumericMatrix Delta1(J, J), Delta2(J, J);
  for (int j = 0; j < J; ++j) {
    NumericVector delta1 = theta[Range((J+j+1)*(J-1), (J+j+2)*(J-1)-1)];
    delta1.push_back(1 - sum(delta1)); // enforcing parameters sum to one
    Delta1(_, j) = delta1;
    NumericVector delta2 = theta[Range((2*J+j+1)*(J-1), (2*J+j+2)*(J-1)-1)];
    delta2.push_back(1 - sum(delta2)); // enforcing parameters sum to one
    Delta2(_, j) = delta2;
  }

  // Adding a penalty for violating positivity
  double pp_rho = sum(pmin(rho,0)*pmin(rho,0));
  NumericVector pp_D1(J), pp_Pi(J), pp_D2(J);
  for (int j = 0; j<J; ++j) {
    pp_D1[j] = sum(pmin(Delta1(_, j), 0) * pmin(Delta1(_, j), 0));
    pp_D2[j] = sum(pmin(Delta2(_, j), 0) * pmin(Delta2(_, j), 0));
    pp_Pi[j] = sum(pmin(Pi_x(_, j), 0) * pmin(Pi_x(_, j), 0));
  }
  double penalty_pos = 0.5 * lambda_pos * (pp_rho + sum(pp_D1) + sum(pp_D2) + sum(pp_Pi));

  // Adding a penalty for violating diagonal dominance
  NumericVector ddy1(J), ddy2(J);
  for (int j = 0; j < J; ++j) {
    ddy1[j] = sum(pmax(Delta1(_, j) - Delta1(j, j), 0) * pmax(Delta1(_, j) - Delta1(j, j), 0));
    ddy2[j] = sum(pmax(Delta2(_, j) - Delta2(j, j), 0) * pmax(Delta2(_, j) - Delta2(j, j), 0));
  }
  double penalty_dd = 0.5 * lambda_dd * (sum(ddy1) + sum(ddy2));

  // Computing the likelihood
  double ll_sum = 0;
  IntegerVector x = tab["x"];
  IntegerVector y1 = tab["y1"];
  IntegerVector y2 = tab["y2"];
  IntegerVector n = tab["n"];
  for (int q = 0; q < J; ++q) {
    for (int r = 0; r < J; ++r) {
      for (int s = 0; s < J; ++s) {
        double p_sum = 0;
        for (int t = 0; t < J; ++t) {
          p_sum += Delta1(r, t) * Delta2(s, t) * Pi_x(t, q) * rho[q];
        }
        // Find n_qrs
        LogicalVector index_qrs = (x == q+1) & (y1 == r+1) & (y2 == s+1);
        int nqrs = 0;
        for(size_t i = 0; i < index_qrs.size(); i++) {
          if(index_qrs[i]) {
            nqrs += n[i];
          }
        }
        // Check if p_sum is zero to avoid NaN
        ll_sum -= nqrs * log(p_sum);
      }
    }
  }

  return ll_sum + penalty_pos + penalty_dd;

}








// [[Rcpp::export]]
NumericVector gradll(NumericVector theta,
                     DataFrame tab,
                     double lambda_pos = 1e6,
                     double lambda_dd = 1e4) {


  // Defining an empty vector for the final gradient
  NumericVector gradient(theta.size());

  // Determining J
  int J = (2 + sqrt(4 + 12 * (theta.size() + 1))) / 6;

  // Extracting rho
  NumericVector rho = theta[Range(0, J - 2)];
  rho.push_back(1 - sum(rho)); // enforcing that parameters sum to one

  // Extracting Pi_X
  NumericMatrix Pi_x(J, J);
  for (int j = 0; j < J; ++j) {
    NumericVector pi_x = theta[Range((j+1)*(J-1), (j+2)*(J-1)-1)];
    pi_x.push_back(1 - sum(pi_x)); // enforcing parameters sum to one
    Pi_x(_, j) = pi_x;
  }

  // Extracting Delta_1 and Delta_2
  NumericMatrix Delta1(J, J), Delta2(J, J);
  for (int j = 0; j < J; ++j) {
    NumericVector delta1 = theta[Range((J+j+1)*(J-1), (J+j+2)*(J-1)-1)];
    delta1.push_back(1 - sum(delta1)); // enforcing parameters sum to one
    Delta1(_, j) = delta1;
    NumericVector delta2 = theta[Range((2*J+j+1)*(J-1), (2*J+j+2)*(J-1)-1)];
    delta2.push_back(1 - sum(delta2)); // enforcing parameters sum to one
    Delta2(_, j) = delta2;
  }

  // Computing the likelihood p for each qrs
  IntegerVector x = tab["x"];
  IntegerVector y1 = tab["y1"];
  IntegerVector y2 = tab["y2"];
  IntegerVector n = tab["n"];
  NumericVector p(J*J*J);
  for (int q = 0; q < J; ++q) {
    for (int r = 0; r < J; ++r) {
      for (int s = 0; s < J; ++s) {
        double p_sum = 0;
        for (int t = 0; t < J; ++t) {
          p_sum += Delta1(r, t) * Delta2(s, t) * Pi_x(t, q) * rho[q];
        }
        // Find indices where the conditions are met
        LogicalVector index_qrs = (x == q+1) & (y1 == r+1) & (y2 == s+1);
        // Record the value of p
        for(size_t i = 0; i < index_qrs.size(); i++) {
          if(index_qrs[i]) {
            p[i]=p_sum;
          }
        }
      }
    }
  }

  // Computing the gradient with respect to rho
  // Looping through X to compute gradient w.r.t. rho
  for (int q = 0; q < J - 1; ++q) {
    // Gradient of the likelihood
    double grad_ll = 0;
    double grad_pos = 0;
    for (int i = 0; i < n.size(); ++i) {
      if (x[i] == q + 1) {
        grad_ll += n[i] / rho[q];
      }
      if (x[i] == J) {
        grad_ll -= n[i] / rho[J - 1];
      }
    }

    // Gradient of the positivity penalty
    grad_pos = -lambda_pos * (rho[q] * (rho[q] < 0) - rho[J - 1] * (rho[J - 1] < 0));

    // Recording the value of the gradient
    gradient[q] = grad_ll + grad_pos;
  }

  // Looping through Y* and X to compute gradient w.r.t. Pi_x
  for (int q = 0; q < J; ++q) {
    for (int t = 0; t < J - 1; ++t) {
      double grad_ll = 0;
      for (int r = 0; r < J; ++r) {
        for (int s = 0; s < J; ++s) {
          for (int i = 0; i < n.size(); ++i) {
            if (x[i] == q + 1 && y1[i] == r + 1 && y2[i] == s + 1) {
              grad_ll += rho[q] * n[i] / p[i] * (Delta1(r, t) * Delta2(s, t) - Delta1(r, J - 1) * Delta2(s, J - 1));
            }
          }
        }
      }

      double grad_pos = -lambda_pos * (Pi_x(t, q) * (Pi_x(t, q) < 0) - Pi_x(J - 1, q) * (Pi_x(J - 1, q) < 0));
      gradient[(q+1)*(J - 1) + t] = grad_ll + grad_pos;
    }
  }

  // Looping through Y1 and Y* to compute gradient w.r.t. Delta1
  for(int t = 0; t < J; t++) {
    for(int r = 0; r < J-1; r++) {
      double grad_ll = 0;
      for(int q = 0; q < J; q++) {
        for(int s = 0; s < J; s++) {
          // Filtering the data
          LogicalVector condition = (x == q+1) & (y1 == r+1) & (y2 == s+1);
          LogicalVector conditionJ = (x == q+1) & (y1 == J) & (y2 == s+1);
          double n_sum = 0, p_sum = 0, nJ_sum = 0, pJ_sum = 0;
          for(size_t i = 0; i < condition.size(); i++) {
            if(condition[i]) {
              n_sum += n[i];
              p_sum += p[i];
            }
            if(conditionJ[i]) {
              nJ_sum += n[i];
              pJ_sum += p[i];
            }
          }
          double n_val = n_sum / p_sum;
          double nJ_val = nJ_sum / pJ_sum;
          grad_ll += rho[q] * Delta2(s, t) * Pi_x(t, q) * (n_val - nJ_val);
        }
      }

      double grad_pos = -lambda_pos * (Delta1(r, t) * (Delta1(r, t) < 0) -
                                       Delta1(J-1, t) * (Delta1(r, t) < 0));

      double grad_dd = -lambda_dd * ((Delta1(r, t) - Delta1(t, t)) * (Delta1(r, t) > Delta1(t, t)) -
                                     (Delta1(J-1, t) - Delta1(t, t)) * (Delta1(J-1, t) > Delta1(t, t)));

      gradient[(J + t + 1)*(J - 1) + r ] = grad_ll + grad_pos + grad_dd;
    }
  }

  // Looping through Y2 and Y* to compute gradient w.r.t. Delta2
  for(int t = 0; t < J; t++) {
    for(int s = 0; s < J-1; s++) {
      double grad_ll = 0;
      for(int q = 0; q < J; q++) {
        for(int r = 0; r < J; r++) {
          // Filtering the data
          LogicalVector condition = (x == q+1) & (y1 == r+1) & (y2 == s+1);
          LogicalVector conditionJ = (x == q+1) & (y1 == r+1) & (y2 == J);
          double n_sum = 0, p_sum = 0, nJ_sum = 0, pJ_sum = 0;
          for(size_t i = 0; i < condition.size(); i++) {
            if(condition[i]) {
              n_sum += n[i];
              p_sum += p[i];
            }
            if(conditionJ[i]) {
              nJ_sum += n[i];
              pJ_sum += p[i];
            }
          }
          double n_val = n_sum / p_sum;
          double nJ_val = nJ_sum / pJ_sum;
          grad_ll += rho[q] * Delta1(r, t) * Pi_x(t, q) * (n_val - nJ_val);
        }
      }

      double grad_pos = -lambda_pos * (Delta2(s, t) * (Delta2(s, t) < 0) -
                                       Delta2(J-1, t) * (Delta2(s, t) < 0));

      double grad_dd = -lambda_dd * ((Delta2(s, t) - Delta2(t, t)) * (Delta2(s, t) > Delta2(t, t)) -
                                     (Delta2(J-1, t) - Delta2(t, t)) * (Delta2(J-1, t) > Delta2(t, t)));

      gradient[(2*J + t + 1)*(J - 1) + s] = grad_ll + grad_pos + grad_dd;
    }
  }

  // Flipping the sign of the gradient
  gradient = -1*gradient;

  return gradient;

}

