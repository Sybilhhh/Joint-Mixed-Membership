#set.seed(9397)
library("mvtnorm")
library("numDeriv")
library("Rcpp")
library("RcppArmadillo")
library("RcppEigen")
library("LaplacesDemon")

sourceCpp("mcmc.cpp")

####################################### Simulation setting and control ###################################################
REP= 1    ## number of replication
TT= 6000 ## number of iteration of mcmc
S= 3000 #T-S is the number of burn-in iteration
BASELINE_DEGREE = 0 # Baseline hazard degree: 0 - constant; 1 - linear; 2 - quadratic
IF_QUADRATIC = 0 # 1 - quadratic polynomial; 0 - natural cubic spline; 
IF_MMM = 1 # 1 - the model has MMM part; 0 - the model does not have MMM part
IF_COX = 1 # 1 - the model has COX part; 0 - the model does not have COX part
IF_AIC = 1 # 1 - compute AIC; 0 - not compute AIC and will be faster
CR = 30 # censoring rate setting:  50 - 50%; others - 30%

source("generate_data.r")

#######################################################################################################################
# MCMC
start_time = Sys.time()
for(iter in 1:REP){
  num_list = 11
  ####### initial value #######
  if(IF_COX == 0){
    lambda_cpp = rep(1,J_len-1)
    c_cpp = c(1,-1)
    beta_cpp = c(1,-1,1)
  }else{
    lambda_cpp = rep(1,J_len-1)
    c_cpp = rep(0.5,NY)
    beta_cpp = rep(0,NZ)
  }
  
  # --------------- MMM ---------------#
  if(IF_MMM == 0){
    b_cpp = as.numeric(matrix(c(-5,-1,4,
              -2,-1,3,
              -1,-1,2,
              -2,1,4,
              -3,1,2,
              -2,1,2), nrow = NP*NY, ncol = M_mcmc, byrow = T))
    psi_epscpp = rep(0.3,NY)
    delta_0_cpp = 0.3
    xi_cpp = c(0.3,0.3,0.4)
    g_cpp = as.numeric(true_g[[iter]][1:N_train,])
  }else{
    b_cpp =as.numeric(matrix(c(-1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2,
                               -1,-1,2), nrow = NP*NY, ncol = M_mcmc, byrow = T))
    psi_epscpp = rep(0.1,NY)
    delta_0_cpp = 0.5
    xi_cpp = rep(1,M_mcmc)/M_mcmc
    g_cpp = as.numeric(rdirichlet(N_train,delta_0_cpp*xi_cpp))
  }
  
  
  ####### prior distribution ########
  #---------- Cox --------------------#
  # lambda
  alpha1 = 2
  alpha2 = 0.01
  # beta
  beta0 = rep(0,NZ)
  Sigma_beta0 = diag(NZ)
  # c
  c0 = 0
  sigma_c0 = 1
  # --------------- MMM ---------------#
  # b
  b0 = rep(0,NP)
  Sigma_b0 = diag(NP)*100
  # delta
  a0_delta = 1
  b0_delta = 5
  # g
  d0 = 1
  # psi_eps
  alpha_eps0 = 9
  beta_eps0 = 4
  
  prior = list()
  # Cox
  prior[[1]] = alpha1
  prior[[2]] = alpha2
  prior[[3]] = beta0
  prior[[4]] = Sigma_beta0
  prior[[5]] = c0
  prior[[6]] = sigma_c0
  # MMM
  prior[[7]] = b0
  prior[[8]] = Sigma_b0
  prior[[9]] = a0_delta
  prior[[10]] = b0_delta
  prior[[11]] = d0
  prior[[12]] = alpha_eps0
  prior[[13]] = beta_eps0
  
  
  para = list()
  # Cox
  para[[1]] = lambda_cpp
  para[[2]] = beta_cpp
  para[[3]] = c_cpp
  # MMM
  para[[4]] = b_cpp
  para[[5]] = g_cpp
  para[[6]] = psi_epscpp
  para[[7]] = delta_0_cpp
  para[[8]] = xi_cpp
  
  result = mcmc(observed_data[((iter-1)*num_list+1):((iter-1)*num_list+11)], para, mh, prior, K, iter, REP, TT, S, pi, IF_MMM, IF_COX, IF_AIC, IF_QUADRATIC)
  AIC[iter] = mean(result$AIC)
  
   for(i in 1:S){
    b_mcmc[i,,] = matrix(result$b_str[((i-1)*NP*NY*M_mcmc+1):(i*NP*NY*M_mcmc)],nrow=NP*NY,byrow=F)
    g_mcmc[i,,] = matrix(result$g_str[((i-1)*N_train*M_mcmc+1):(i*N_train*M_mcmc)],nrow=N_train,byrow=F)
  }
  gb_mcmc = matrix(result$gb,ncol=NY*NP,byrow=F)
  cgb_mcmc = matrix(result$cgb,ncol=NP,byrow=F)
  delta_0_mcmc = result$delta_0_str
  xi_mcmc = matrix(result$xi_str,ncol=M_mcmc,byrow=T)
  psi_epsmcmc = matrix(result$psieps_str,nrow=S,byrow=T) 
  
  b_rep[iter,,] = apply(b_mcmc,2:3,mean)
  g_rep[iter,,] = apply(g_mcmc,2:3,mean)
  xi_rep[iter,] = apply(xi_mcmc,2,mean)
  delta_0_rep[iter] = mean(delta_0_mcmc)
  psi_epsrep[iter,] = apply(psi_epsmcmc,2,mean)
  
  # Cox
  c_mcmc = matrix(result$c_str,nrow=S,byrow=T) 
  beta_mcmc = matrix(result$beta_str,nrow=S,byrow=T) 
  lambda_mcmc = matrix(result$lambda_str,nrow=S,byrow=T)
  
  c_rep[iter,] = apply(c_mcmc,2,mean)
  beta_rep[iter,] = apply(beta_mcmc,2,mean)
  lambda_rep[iter,] = apply(lambda_mcmc,2,mean)
}

Iter_time = Sys.time() - start_time
sprintf("Need %.3f mins to iter %d replication", Iter_time, REP)

####################################### Result analysis ###################################################
apply(b_rep,2:3,mean)
mean(delta_0_rep)
apply(xi_rep,2,mean)
apply(psi_epsrep,2,mean)

apply(c_rep,2,mean)
apply(beta_rep,2,mean)
apply(lambda_rep,2,mean)

mean(AIC) 

result$at_b/TT
result$at_delta/TT
result$at_c/TT
result$at_beta/TT
mean(result$at_g/TT)

write.table(c_rep,"c_rep.txt")
write.table(beta_rep,"beta_rep.txt")
write.table(lambda_rep,"lambda_rep.txt")
write.table(b_rep,"b_rep.txt")
write.table(delta_0_rep,"delta_0_rep.txt")
write.table(xi_rep,"xi_rep.txt")
write.table(psi_epsrep,"psi_epsrep.txt")
write.table(g_rep,"g_rep.txt")
write.table(AIC,"AIC.txt")


save.image(file = "1_const_NP6_N200.RData")


