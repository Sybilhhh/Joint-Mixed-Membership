N=300 # number of patients
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
'b <- matrix(c(-5,-1,4,
              -2,-1,3,
              -2,1,4,
              -3,1,2), nrow = NP*NY, ncol = M, byrow = T)'
'b <- matrix(c(-5,-1,4,
              -2,-1,3,
              -1,-1,2,
              -2,1,1,
              -1,-1,2,
              -2,1,1,
              
              -2,1,4,
              -3,1,2,
              -2,1,1,
              -1,-1,2,
              -2,1,1,
              -1,-1,2), nrow = NP*NY, ncol = M, byrow = T)'
b <- matrix(c(-5,-1,4,
              -2,-1,3,
              -1,-1,2,
              -2,1,4,
              -3,1,2,
              -2,1,1), nrow = NP_true*NY, ncol = M, byrow = T)
'b <- matrix(c(-2,1,4,
              -3,1,2,
              -2,1,2), nrow = NP*NY, ncol = M, byrow = T)'
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
observed_data = list()
true_g = list()
test_data = list()


'# Find the mode of the time
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Define knots using prior computing
temp_kns <- quantile(OT, probs = seq(0,1,1/NP))
kns <- rep(NA,NP)
for (p in 1:NP) {
  # Find most frequent time between quantiles
  kns[p] <- Mode(OT[OT > temp_kns[p] & OT <= temp_kns[p+1]])
}
'
for(iter in 1:REP){
  
  ############  2. Covariates  #############
  #---------- MMM ---------------
  g = rdirichlet(N,delta)
  gb = matrix(NA, nrow = N, ncol = NP_true*NY)
  cgb = matrix(0, nrow = N, ncol = NP_true)
  NT = numeric(N)
  
  #----------------- Cox ---------------------
  FT = numeric(N) # failure time
  CT = numeric(N) # censoring time 
  OT = numeric(N) # observed time
  DELTA = numeric(N) # failure indicator
  Z = matrix(NA, nrow = N, ncol = NZ) # covariates
  Z[,1] = rexp(N,1)-1
  Z[,2] = rnorm(N,0,1)
  Z[,3] = rt(N,5)
  
  Y = c()
  ft = c()
  Time = c() # measurement time
  Queue = c(0)
  
  ############  3. Sample data  #############
  for(i in 1:N){
    for (k in 1:(NY*NP_true)) {
      gb[i,k] = g[i,] %*% b[k,] # useful in computing
    }
    
    for (p in 1:NP_true) {
      for (k in 1:NY) {
        cgb[i,p] = cgb[i,p] + c[k] * gb[i,((k-1)*NP_true+p)]
      }
    }
    
    #----------------- Cox time ---------------------
    temp_ratio = -log(1-runif(1))/exp(Z[i,]%*%beta)[1,1]
    rand = seq(0,10,0.001)
    
    f_rand = 0 # h(T_{ij})
    'f_rand = f_rand + cgb[i,1]
    if(NP>1){
      f_rand = f_rand + cgb[i,2]*rand
    }
    if(NP>2){
      if(IF_QUADRATIC==1){
        f_rand = f_rand + cgb[i,3]*rand^2
      }else{
        for(p in 3:NP){
          f_rand = f_rand + cgb[i,p] *(-1)* ((abs(rand-kns[NP])^3+(rand-kns[NP])^3-(abs(rand-kns[p-2])^3+(rand-kns[p-2])^3))/(2*(kns[NP]-kns[p-2]))-
                                               (abs(rand-kns[NP])^3+(rand-kns[NP])^3-(abs(rand-kns[NP-1])^3+(rand-kns[NP-1])^3))/(2*(kns[NP]-kns[NP-1])))
        }
      }
    }'
    
    f_rand = f_rand + cgb[i,1]
    if(NP_true>1){
      f_rand = f_rand + cgb[i,2]*rand
    }
    if(NP_true>2){
      f_rand = f_rand + cgb[i,3]*rand^3
    }
    
    # FT - true failure time; CT - censoring time; OT - observed time
    if(BASELINE_DEGREE == 0){ # baseline harzard function lambda_0 = 1
      integration = cumsum(exp(f_rand)*0.001) 
      FT[i] = which(temp_ratio<integration)[1]*0.001
      if(CR == 50){
        CT[i] = runif(1,1,1.8) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
      }else{
        CT[i] = runif(1,1.6,2) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
      }
    }
    if(BASELINE_DEGREE == 1){ # baseline harzard function lambda_0 = 1 + t
      integration = cumsum((rand + 1)*exp(f_rand)*0.001) 
      FT[i] = which(temp_ratio<integration)[1]*0.001
      if(CR == 50){
        CT[i] = runif(1,0.8,1.8) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
      }else{
        CT[i] = runif(1,1.6,2) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
      }
    }
    if(BASELINE_DEGREE == 2){ # baseline harzard function lambda_0 = 1 + t^2
      integration = cumsum((0.5*rand^2 + 1)*exp(f_rand)*0.001) 
      FT[i] = which(temp_ratio<integration)[1]*0.001
      if(CR == 50){
        CT[i] = runif(1,0.8,1.8) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
      }else{
        CT[i] = runif(1,1.6,2) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
      }
    }
    
    if(is.na(FT[i])){
      DELTA[i] = 0
      OT[i] = CT[i]
    }else if(FT[i]<CT[i]){
      DELTA[i] = 1
      OT[i] = FT[i]
    }else{
      DELTA[i] = 0
      OT[i] = CT[i]
    }# DELTA - failure indicator
    
    #---------- MMM time ---------------
    # generate measurement time
    'if(DELTA[i]==1){
      NT[i] = sample(6:MT,1)
      #Time_temp = c(sort(sample(setdiff(seq(cut,OT[i],cut)),NT[i]-2)),OT[i])
      Time_temp = sort(c(OT[i],sample(setdiff(seq(cut,trunc(CT[i]/cut)*cut,cut),OT[i]),NT[i]-2)))
      #Time_temp = sort(c(OT[i],sample(setdiff(seq(0.25,2,0.2),OT[i]),NT[i]-1)))
      Time = c(Time,0,Time_temp)
    }else{
      NT[i] = sample(6:MT,1)
      Time_temp = c(sort(sample(seq(cut,OT[i],cut),NT[i]-2)),OT[i])
      #Time_temp = sort(c(OT[i],sample(setdiff(seq(0.25,2,0.2),OT[i]),NT[i]-1)))
      Time = c(Time,0,Time_temp)
    }
    '
    NT[i] = sample(6:MT,1)
    Time_temp = sort(c(OT[i],sample(setdiff(seq(0.1,2,0.1),OT[i]),NT[i]-2)))
    Time = c(Time,0,Time_temp)
    
    # generate response variable
    Queue=c(Queue,sum(NT))  # No. of the first measurement of each patient (used in Rcpp)
    Y_temp = matrix(NA,nrow = NT[i], ncol = NY)
    for(j in 1:NT[i]){
      eps = as.vector(rmvnorm(1,rep(0,NY),as.matrix(psi_eps)))
      
      'for(k in 1:NY){
        h_temp = gb[i,(1+(k-1)*NP)]
        if(NP>1){
          h_temp = h_temp + gb[i,(2+(k-1)*NP)]*Time[Queue[i]+j]
        }
        if(NP>2){
          if(IF_QUADRATIC==1){
            h_temp = h_temp + gb[i,(3+(k-1)*NP)]*Time[Queue[i]+j]^2
          }else{
            for(p in 3:NP){
              h_temp = h_temp + gb[i,(p+(k-1)*NP)] *(-1)* ((abs(Time[Queue[i]+j]-kns[NP])^3+(Time[Queue[i]+j]-kns[NP])^3-(abs(Time[Queue[i]+j]-kns[p-2])^3+(Time[Queue[i]+j]-kns[p-2])^3))/(2*(kns[NP]-kns[p-2]))-
                                                             (abs(Time[Queue[i]+j]-kns[NP])^3+(Time[Queue[i]+j]-kns[NP])^3-(abs(Time[Queue[i]+j]-kns[NP-1])^3+(Time[Queue[i]+j]-kns[NP-1])^3))/(2*(kns[NP]-kns[NP-1])))
            }
          }
        }
        Y_temp[j,k] = h_temp + eps[k]
      }'
      
      for(k in 1:NY){
        h_temp = gb[i,(1+(k-1)*NP_true)]
        if(NP_true>1){
          h_temp = h_temp + gb[i,(2+(k-1)*NP_true)]*Time[Queue[i]+j]
        }
        if(NP_true>2){
          h_temp = h_temp + gb[i,(3+(k-1)*NP_true)]*Time[Queue[i]+j]^3
        }
        Y_temp[j,k] = h_temp + eps[k]
      }
    }
    Y = rbind(Y,Y_temp)
  }
  
  cr = 1-sum(DELTA)/N
  sprintf("censor rate: %.3f", cr)
  SNT = sum(NT)
  Queue=Queue[-(N+1)]  # No. of the first measurement of each patient (used in Rcpp)
  
  kns = quantile(Time, probs = seq(0,1,1/(NP+1)))[2:(NP+1)]
  # Basis spline for each subject i in time j
  B = array(0,dim=c(SNT,NP))
  B[,1] = 1
  if(NP>1){
    B[,2] = Time
  }
  if(NP>2){
    if(IF_QUADRATIC==1){
      B[,3] = Time^2
    }else{
      for(p in 3:NP){
        B[,p] = (-1)*((abs(Time-kns[NP])^3+(Time-kns[NP])^3-(abs(Time-kns[p-2])^3+(Time-kns[p-2])^3))/(2*(kns[NP]-kns[p-2]))-
                        (abs(Time-kns[NP])^3+(Time-kns[NP])^3-(abs(Time-kns[NP-1])^3+(Time-kns[NP-1])^3))/(2*(kns[NP]-kns[NP-1])))
      }
    }
  }
  
  B_train_cpp = as.vector(B[1:Queue[N_train+1],]) # 1:SNT is first time 1, (SNT+1):2SNT is second time t
  
  # labels of survival time in measurement time, used in computing
  OT_label = rep(NA,N)
  for(k in 1:N){
    OT_label[k] = which(Time[(Queue[k]+1):(Queue[k]+NT[k])]==OT[k])-1
  }
  
  # intervals of baseline harzard function
  J_cpp = c(0,as.vector(quantile(OT,seq(0,1,0.2)))[-c(1,6)],as.vector(quantile(OT,seq(0,1,0.2)))[6]+0.001) # grids lambda_h
  Fail_label = numeric(N) # v_ih
  for(i in 1:N){
    Fail_label[i] = which(!(OT[i]>J_cpp))[1]-1
  }
  y_train_cpp = as.vector(Y[1:Queue[N_train+1],]) # 1:SNT is first response y_ij1, (SNT+1):2SNT is second response y_ij2
  z_train_cpp = as.vector(Z[1:N_train,]) # 1:N_train is first covariates Z1, (N_train+1):2N_train is second covariates Z2
  
  y_test_cpp = Y[-(1:Queue[N_train+1]),]
  z_test_cpp = Z[-(1:N_train),]
  
  ############  4. Save data  #############
  num_list = 11
  
  observed_data[[(iter-1)*num_list+1]] = OT_label[1:N_train]
  observed_data[[(iter-1)*num_list+2]] = B_train_cpp
  observed_data[[(iter-1)*num_list+3]] = J_cpp
  observed_data[[(iter-1)*num_list+4]] = Fail_label[1:N_train]
  observed_data[[(iter-1)*num_list+5]] = z_train_cpp
  observed_data[[(iter-1)*num_list+6]] = OT[1:N_train]
  observed_data[[(iter-1)*num_list+7]] = DELTA[1:N_train]
  observed_data[[(iter-1)*num_list+8]] = Queue[1:N_train]
  observed_data[[(iter-1)*num_list+9]] = y_train_cpp
  observed_data[[(iter-1)*num_list+10]] = NT[1:N_train]
  observed_data[[(iter-1)*num_list+11]] = kns
  
  true_g[[iter]] = g
  
  num_list_test = 9
  
  test_data[[(iter-1)*num_list_test+1]] = OT_label[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+2]] = t(B[-(1:Queue[N_train+1]),])
  test_data[[(iter-1)*num_list_test+3]] = Fail_label[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+4]] = t(z_test_cpp)
  test_data[[(iter-1)*num_list_test+5]] = OT[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+6]] = DELTA[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+7]] = Queue[-(1:N_train)]-Queue[-(1:N_train)][1]
  test_data[[(iter-1)*num_list_test+8]] = t(y_test_cpp)
  test_data[[(iter-1)*num_list_test+9]] = NT[-(1:N_train)]
  
}

saveRDS(observed_data, file="observed_data.RData")
saveRDS(test_data, file="test_data.RData")
saveRDS(true_g, file="true_g.RData")

