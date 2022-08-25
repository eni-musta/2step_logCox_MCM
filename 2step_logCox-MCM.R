####################################################################################
# A 2-step approach approach for estimation in the logistic/Cox mixture cure model #
####################################################################################

# logistic/Cox model
# Incidence:  P('cured'|X)=1-phi(gamma_0,X),   phi(gamma,X)=1/(1+exp(-gamma*X)),
#             where gamma_0 is the vector of unknown parameters and X is the vector of covariates
# Latency:  for the uncured subjects, S_u(t|Z)=S_0(t)^exp(beta*Z),
#           where beta_0 is a vectorof unknown parameters and S_0 is the baseline survival function


# Libraries

library(smcure)
library(np)
library(survival)

# Predefined functions

source("functions_2step_logCox_MCM")

# Load the data into the dataframe "Data" and remove the observations with missing values. 
# The dataframe consists of the following columns: 
# Y - follow-up time
# Delta - censoring status:  censored (Delta=0), noncensored (Delta=1)
# Treatment - a binary covariate (0=control, 1=treatment)
# Age - a continuous covariate centered to the mean
# Gender - a binary covariate (0=male, 1=female)

data("e1684")
Data = e1684   
Data = na.omit(Data) 
colnames(Data) = c("Treatment","Y","Delta","Age","Sex")


# Plot the Kaplan Meier estimator of the survival function

plot(survfit(Surv(Data$Y,Data$Delta)~1),xlab = 'Time (years)',ylab = 'Survival probability',main='Kaplan-Meier estimator')


# Compute the maximum likelihood estimators with the package smcure

MCM = smcure(Surv(Y,Delta)~Age+Treatment+Sex,cureform=~Age+Treatment+Sex,data=Data,model="ph",Var=TRUE,nboot=500)

#smcure estimators
gamma_smcure = MCM$b
beta_smcure = MCM$beta
s_smcure = MCM$s


###################################################
# 2-step approach for estimation of the incidence #
###################################################


# Step 1. Compute a nonparametric estimator conditional on a one-dimensional covariate 
#         constructed using the initial estimator 
    

# data_log is a dataframe containing the follow-up times, censoring status and the covariates used to model the incidence
# In this case all covariates are used for both incidence and latency 
# The observations (rows) are ordered according to the follow-up time.

  data_log = data.frame(cbind(Data$Y,Data$Delta,Data$Age,Data$Treatment,Data$Sex))
  colnames(data_log) = c('Y','Delta','Age','Treatment','Sex')
  ord = order(data_log[,1],1-data_log[,2])
  data_log = data_log[ord,]
  
  Y_r = max(data_log[which(data_log[,2]==1),1])    # the largest observed event time
  
  gamma=gamma_smcure # initial estimator
  m=length(gamma)
  n=dim(data_log)[1]
  z=gamma[1]+c(gamma[2:m]%*%t(data_log[,3:(m+1)])) # new 1-dim covariate
  z=(z-mean(z))/sd(z) # The obtained 1-dim covariate is standardized 

  hopt= npcdistbw(formula=data_log[,1]~z,gydat=data_log[which(data_log[,1]<=Y_r),1],cxkertype="epanechnikov")$xbw  # compute bandwidth for estimation of  H(t|z)
  
  
  w = weight_1dim(gamma,data_log[,(3:(m+1))],hopt)  
  pnonp = beran(data_log,w)  
  
# Step 2. Project the nonparametric estimator to the desired parametric class (logistic) using a Bernoulli-type likelihood
  
  gamma_1dim=optim(par=gamma,Lik_logit_1dim,V=data_log,pnonp=pnonp,method='BFGS',control=list(fnscale=-1))$par

#########################################################################
  
# Latency estimation.
# Estimate the regression coefficients and the baseline hazard of the Cox component 
  
  
  n = dim(Data)[1]  #sample size
  X = as.matrix(cbind(rep(1, n), Data$Age, Data$Treatment, Data$Sex)) # covariates for the incidence including the intercept 
  Z = as.matrix(cbind(Data$Age, Data$Treatment, Data$Sex)) #covariates for the latency
  Time = Data$Y          # Follow-up time
  Status = Data$Delta    # Censoring status
  
  beta_new=em2(Time, Status, X, Z,gamma=gamma_1dim,beta=beta_smcure,emmax=100,eps=1e-7)$beta
  
########################################################################  
# Estimate the variance via bootstrap 
  
  nboot = 500        # number of bootstrap samples
  nbeta = dim(Z)[2]  # number of covariates for the latency
  ngamma = dim(X)[2] # number of covariates for the incidence
  log_var=c('Intercept','Age','Treatment','Sex')
  cox_var=c('Age','Treatment','Sex')
    
  gamma_boot = matrix(rep(0, nboot * ngamma), nrow = nboot)
  beta_boot = matrix(rep(0, nboot * nbeta), nrow = nboot)
  gamma_boot_EM = matrix(rep(0, nboot * ngamma), nrow = nboot)
  beta_boot_EM = matrix(rep(0, nboot * nbeta), nrow = nboot)
  
    iter = matrix(rep(0, nboot), ncol = 1)
   
  tempdata = cbind(Time, Status, X, Z)
  data1 = subset(tempdata, Status == 1)
  data0 = subset(tempdata, Status == 0)
  n1 = nrow(data1)
  n0 = nrow(data0)
    
  i = 1
  while (i <= nboot) {
    
    # generate the bootstrap sample
    
      id1 = sample(1:n1, n1, replace = TRUE)
      id0 = sample(1:n0, n0, replace = TRUE)
      bootdata = rbind(data1[id1, ], data0[id0, ])
      bootZ = bootdata[,(3+ngamma):(2+ngamma+nbeta)]
      bootX = bootdata[,3:(2+ngamma)]
      sd_age=sd(bootX[,2]) 
      bootX[,2]=bootX[,2]/sd_age
      bootZ[,1]=bootZ[,1]/sd_age
      bootdata[,4]=bootdata[,4]/sd_age
      bootdata[,7]=bootdata[,7]/sd_age
      
    # estimate the parameters of the incidence  
      
      data_log.boot = data.frame(cbind(bootdata[,1],bootdata[,2],bootX[,-1]))
      colnames(data_log.boot) = c('Y','Delta',log_var[-1])
      ord.boot = order(data_log.boot[,1],1-data_log.boot[,2])
      data_log.boot = data_log.boot[ord.boot,]
      
      Y_r.boot = max(data_log.boot[which(data_log.boot[,2]==1),1])    # the largest observed event time
      
      
      beta.hat.boot <- coxph(Surv(Y,Delta)~Age+Treatment+Sex,subset =Delta!=0,data =  data_log.boot ,method = "breslow")$coef 
      gamma.hat.boot <- glm(Delta~Age+Treatment+Sex,family=binomial(link=logit),data= data_log.boot)$coefficients
      
      MCM.boot=em1(bootdata[,1],bootdata[,2],bootX,bootZ,gamma.hat.boot,beta.hat.boot,emmax=100,eps=1e-7)
      gamma_boot_EM[i,]=MCM.boot$b/c(1,sd_age,1,1)
      
      convergence_EM.boot=as.numeric(MCM.boot$tau<1e-7)
      
      if (convergence_EM.boot==1){
        
      gamma=MCM.boot$b #gamma_smcure
      m=length(gamma)
      z.boot=gamma[1]+c(gamma[2:m]%*%t(data_log.boot[,3:(m+1)]))
      z.boot=(z.boot-mean(z.boot))/sd(z.boot)
      
      hopt.boot=hopt 
      w.boot = weight_1dim(gamma,data_log.boot[,(3:(m+1))],min(hopt.boot,2))  
      pnonp.boot = beran(data_log.boot,w.boot)  
      gamma_1dim.boot=optim(par=gamma,Lik_logit_1dim,V=data_log.boot,pnonp=pnonp.boot,method='BFGS',control=list(fnscale=-1))$par
      
      #estimate latency
      bootfit=em2(bootdata[, 1], bootdata[, 2], bootX, bootZ, gamma=gamma_1dim.boot,beta=beta_smcure,emmax=100,eps=1e-7)
      
      gamma_boot[i, ] = gamma_1dim.boot/c(1,sd_age,1,1)
      beta_boot[i, ] = bootfit$beta/c(sd_age,1,1)
      
      i = i + 1}
    }
    
    gamma_var = apply(gamma_boot, 2, var)
    beta_var = apply(beta_boot, 2, var)
    gamma_sd = sqrt(gamma_var)
    beta_sd = sqrt(beta_var)
  
  fit = list()
  fit$beta = beta_new
  fit$gamma = gamma_new
  fit$gamma_var = gamma_var
  fit$gamma_sd = gamma_sd
  fit$gamma_zvalue = gamma_new/gamma_sd
  fit$gamma_pvalue = (1 - pnorm(abs(fit$gamma_zvalue)))*2
  fit$beta_var = beta_var
  fit$beta_sd = beta_sd
  fit$beta_zvalue = beta_new/beta_sd
  fit$beta_pvalue = (1 - pnorm(abs(fit$beta_zvalue)))*2
  fit$gammanm =  log_var
  fit$betanm = cox_var
  print_logCox_MCM(fit)  
