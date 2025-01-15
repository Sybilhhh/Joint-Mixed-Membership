N= 300 # number of patients
N_test = 100
N_train = N-N_test
M = 3 # number of extreme profile
MT = 12 # maximum number of repeated measurement of each patient
NY=2 # dimension of y
NZ=3 # number of covariates
NP = 6 # number of spline basis, if quadratic, then must be NP = 3
J_len = 6
M_mcmc = 3
NP_true = 3

#################### 1.True value and Parameters #########################################################
#---------- MMM ---------------
# The coefficient of basis, each column represent three groups, 
# each row represent coefficient of each basis
# the 1st to NP-th row is the first response, (NP+1) to (2*NP) is second
b <- matrix(c(-5,-1,4,
              -2,-1,3,
              -1,-1,2,
              -2,1,4,
              -3,1,2,
              -2,1,1), nrow = NP_true*NY, ncol = M, byrow = T)
delta_0 = 0.3 # spread parameter
xi = c(0.3,0.3,0.4) # average proportion
delta = delta_0 * xi # concentration parameter
psi_eps = diag(rep(0.3,NY+1))[1:NY,1:NY]
cut = 0.02


#----------------- Cox ---------------------
beta = c(1,-1,1) # coefficient of covariates
c = c(1,-1) # coefficient of f(t)*u_i
#c = -1
K = 10 # number of grids used in quadrature method

#----------------- Save data ----------------------#
#---------- Mixed Membership model---------------#
b_mcmc = array(0, dim = c(S,NP*NY,M_mcmc))
g_mcmc = array(0,dim=c(S,N_train,M_mcmc))
delta_0_mcmc = rep(0,S)
xi_mcmc = array(0, dim = c(S,M_mcmc))
psi_epsmcmc = array(0,dim=c(S,NY))

b_rep = array(0,dim = c(REP,NP*NY,M_mcmc))
g_rep = array(0,dim = c(REP,N_train,M_mcmc))
delta_0_rep = rep(0,REP)
xi_rep = array(0,dim = c(REP,M_mcmc))
psi_epsrep = array(0,dim=c(REP,NY))

#----------------- Cox model ---------------------#
c_mcmc = array(0,dim=c(S,NY))
beta_mcmc = array(0,dim=c(S,NZ))
lambda_mcmc = array(0,dim=c(S,J_len-1))

c_rep = array(0,dim=c(REP,NY))
beta_rep = array(0,dim=c(REP,NZ))
lambda_rep = array(0,dim=c(REP,J_len-1))

AIC = rep(0,REP)

c_beta = 5.5
c_c = 0.1
c_delta = 0.08
c_b = 5
c_g = 0.1

mh = list()
mh[[1]] = c_beta
mh[[2]] = c_c
mh[[3]] = c_delta
mh[[4]] = c_b
mh[[5]] = c_g



######################## generate data ##################################################
observed_data = readRDS("observed_data.RData")

