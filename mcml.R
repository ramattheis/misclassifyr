#---------------------------------------------------------------------
# summary: This script estimates parameters in a misclassification 
#      problem with two noisy measures of a latent outcome and one
#      accurate regressor, all discrete with dimension J, via ML, 
#      using custom C++ functions defined in helper.cpp
# author: Ross Mattheis
# created: February 14, 2024
# last modified: February 29, 2024
# notes:
#---------------------------------------------------------------------

# Making sure data.table is loaded 
require(data.table)

# Loading a C++ functions to quickly evaluate the likelihood and gradient
Rcpp::sourceCpp(file.path(projPath,"code","helper.cpp"))

#------------------------------------------------------
# Estimation via ML
#------------------------------------------------------

mcml = function(data, 
                maximum_iterations = 1e4, 
                tolerance = 1e-6){
  
  # Making sure data is a data.table 
  data = as.data.table(data)
  
  ## Determining J
  J = as.integer(round(nrow(data)^(1/3),0))
  
  ## Making sure that the elements of the tabulation are positive 
  if(any(data$n == 0)){
    data$n = data$n + 1
  }
  
  ## Starting rho at the marginal distribution of X 
  data = data[order(data$x),]
  ppx = data[,.(n=sum(n)),by=.(x)]$n
  ppx = ppx/sum(ppx)
  
  ## Initializing the search at an "informed" (diagonal-dominant) point
  ppys_x = diag(rep(1,J)) + 1/sqrt(J)
  ppys_x = apply(ppys_x, 2, function(p) p/sum(p) )
  ppy1_ys = diag(rep(1,J)) + 1/J
  ppy1_ys = apply(ppy1_ys, 2, function(p) p/sum(p) )
  ppy2_ys = diag(rep(1,J)) + 1/J
  ppy2_ys = apply(ppy2_ys, 2, function(p) p/sum(p) )
  
  theta0 = c(ppx[1:(J-1)],
             c(ppys_x[-J,]),
             c(ppy1_ys[-J,]),
             c(ppy2_ys[-J,]))
  
  # quasi-Newton w/ numerical approximation to the gradient
  out = optim(par = theta0,
              fn = ll,
              gr = gradll,
              tab = data,
              method = "BFGS",
              control = list(maxit = maximum_iterations,
                             reltol = tolerance))
  
  # Checking whether convergence is achieved
  if(out$convergence != 0){
    return("Convergence failed, try increasing maximum iterations or tolerance.")
  } else if(out$counts[1] < 200){
    return("Convergence failed, numerical optimization problem.")
  }
  
  #------------------------------------------------------
  # Recording parameter values
  #------------------------------------------------------
  
  # Extracting parameter values
  theta_h = out$par
  
  # Determining J
  J = as.integer((2+sqrt(4+12*(length(theta_h)+1)))/6)
  
  # Extracting rho
  rho = theta_h[1:(J-1)]
  rho[J] = 1 - sum(rho) # enforcing that parameters sum to one
  theta_h = theta_h[-c(1:(J-1))]
  
  # Extracting Pi
  Pi_x = matrix(nrow = J,ncol = 0) # making an empty column
  for(j in 1:J){ # Looping through the rows
    pi_x = theta_h[1:(J-1)] # extracting Pi_X
    pi_x[J] = 1 - sum(pi_x) # enforcing parameters sum to one
    theta_h = theta_h[-c(1:(J-1))] # popping off values from theta
    Pi_x = cbind(Pi_x, pi_x) # stacking the matrix
  }
  
  # Extracting Delta_1
  Delta1 = matrix(nrow = J,ncol = 0) # making an empty column
  for(j in 1:J){ # Looping through the columns
    delta1 = theta_h[1:(J-1)] # extracting M_Y1|Y*=j
    delta1[J] = 1 - sum(delta1) # enforcing parameters sum to one
    theta_h = theta_h[-c(1:(J-1))] # popping off values from theta
    Delta1 = cbind(Delta1, delta1) # stacking the matrix
  }
  
  # Extracting Delta_2
  Delta2 = matrix(nrow = J,ncol = 0) # making an empty column
  for(j in 1:J){ # Looping through the columns
    delta2 = theta_h[1:(J-1)] # extracting M_Y2|Y*=j
    delta2[J] = 1 - sum(delta2) # enforcing parameters sum to one
    theta_h = theta_h[-c(1:(J-1))] # popping off values from theta
    Delta2 = cbind(Delta2, delta2) # stacking the matrix
  }
  
  #--- Transforming parameter values into a data frame
  
  # rho
  out_rho = as.data.frame(rho)
  out_rho$param_name = paste0("Pr(x = ", 1:J,")")
  colnames(out_rho) = c("ml_estimate","param_name")
  
  # Pi
  Pi_x = unname(Pi_x)
  out_pix = reshape2::melt(Pi_x)
  out_pix$param_name = paste0("Pr(Y* = ",out_pix$Var1,
                              "| X = ", out_pix$Var2 ,")" )
  out_pix = out_pix[,c("value","param_name")]
  colnames(out_pix) = c("ml_estimate","param_name")

  # Delta_1
  Delta1 = unname(Delta1)
  out_d1 = reshape2::melt(Delta1)
  out_d1$param_name = paste0("Pr(Y1 = ",out_d1$Var1,
                             "| Y* = ", out_d1$Var2 ,")" )
  out_d1 = out_d1[,c("value","param_name")]
  colnames(out_d1) = c("ml_estimate","param_name")

  # Delta_2
  Delta2 = unname(Delta2)
  out_d2 = reshape2::melt(Delta2)
  out_d2$param_name = paste0("Pr(Y2 = ",out_d2$Var1,
                             "| Y* = ", out_d2$Var2 ,")" )
  out_d2 = out_d2[,c("value","param_name")]
  colnames(out_d2) = c("ml_estimate","param_name")

  # Binding and returning
  out = rbind(out_rho,out_pix, out_d1,out_d2)
  out = out[,c("param_name","ml_estimate")]
  
  return(out)
  
}

