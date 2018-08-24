#### G Magombedze @imperial college 2015
#############################################################################################################################
################ Likelihood Functions For fitting individual data sets or multiple data sets
### 1. Gaussian or LSE (least square estimates)
### 2. Gaussian Likelihood
### 3. Poison Likelihood


##############################################################################################################################################################
##################################### Functions for the model without aestivation ############################################################################
###############################################################################################################################################################

############################################################################################################################################################### 
 #LSE method (fitting active mosquitoes only):: This is similar to using the gaussian likelihood
 LSE_2qn_residuals<-function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step<-length(T_R$Rain)*30.5  #number of time steps/days in a year
  V_0 = Vint   
  
  Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  return(Data$Adult-model_out$Adult) 
  }
  
   
  Gaussian_Log_LKL_2qn<-function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step<-length(T_R$Rain)*30.5  #number of time steps/days in a year
  V_0 = Vint   
  
  Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  Log_LLK<-sum((Data$Adult-model_out$Adult)^2) 
  return(-2*Log_LLK)
  }
  
  
  # Poisson log likelihood:: Fitting active mosquitoes only
Posisson_Log_LKL_2Eqn<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
    
  Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  
  Log_LLK<- sum(Data$Adult*log(model_out$Adult)-model_out$Adult-log(factorial(Data$Adult)))
  #Log_LLK<- sum(Data$Adult*log(model_out$Adult)-model_out$Adult)
  return(-1*Log_LLK)
  }
  
 Gaussian_Aestivation_2Eqn_pop<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
  
  Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
    
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  Var1=mean((Data$Adult1-mean(Data$Adult1))^2)
  Var2=mean((Data$Adult2-mean(Data$Adult2))^2)
  #Var1=mean((Data$Adult1)^2)-(mean(Data$Adult1))^2
  #Var2=mean((Data$Adult2)^2)-(mean(Data$Adult2))^2
  Const1=0.5*length(Data$Adult1)*log(2*pi*Var1)
  Const2=0.5*length(Data$Adult2)*log(2*pi*Var2)
  Gaussian_LLK<-Const1+Const2+sum((Data$Adult1-model_out$Adult)^2)/(2*Var1)+sum((Data$Adult2-model_out$Adult)^2)/(2*Var2)
  return((Gaussian_LLK))
  } 
  
 Log_LKL_poisson_2Eqn_Pop<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))
  
  t_step =14*30.5  #number of time steps/days in a year
  V_0 = Vint 
  
  Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
   
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  Log_LLK_1<- sum(Data$Adult1*log(model_out$Adult)-model_out$Adult-log(factorial(Data$Adult1)))
  Log_LLK_2<- sum(Data$Adult2*log(model_out$Adult)-model_out$Adult-log(factorial(Data$Adult2)))
  return(-1*(Log_LLK_1+Log_LLK_2))
 }
 
 LSE_2qn_residuals_pop<-function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))
  Data1<-as.data.frame(cbind(time=(AG_Villages$Time),Adult=AG_Villages$ControlV1))
  Data2<-as.data.frame(cbind(time=(AG_Villages$Time),Adult=AG_Villages$ControlV2))
  
  t_step<-length(T_R$Rain)*30.5  #number of time steps/days in a year
  V_0 = Vint   
  
  Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
   
  cost1<-modCost(obs=Data1,model=model_out)
  cost2<-modCost(obs=Data2,model=model_out,cost=cost1)
  return(cost2)
  }
  
##############################################################################################################################################################
##################################### Functions for the model with aestivation ###############################################################################
###############################################################################################################################################################

############################################################################################################################################################### 
 #LSE method (fitting active mosquitoes only):: This is similar to using the gaussian likelihood
 LSE_3qn_residuals_active<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
    
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=(Vector_proj[,2])))
    
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  return(Data$Adult-model_out$Adult)
  }

############################################################################################################################################################
#LSE method (fitting combined active and dormant mosquitoes):: This is similar to using the gaussian likelihood

LSE_3qn_Residuals_active_and_dormant_combined<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
    
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=(Vector_proj[,2]+Vector_proj[,3]))) #combine wet and dry season mosquito populations for fitting
  
  ii <- which (out_V$time %in% c(Data$time)) #match data points and simulated time points
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  return(Data$Adult-model_out$Adult) # return residuals
 }
 
  ########################################################################################################################################
  # Poisson log likelihood:: Fitting active mosquitoes only
Posisson_Log_LKL_3Eqn_active<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
    
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=(Vector_proj[,2])))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  min_val=1.0e-7
  
  Log_LLK<- sum(Data$Adult*log(model_out$Adult+min_val)-model_out$Adult-log(factorial(Data$Adult)))
  #Log_LLK<- sum(Data$Adult*log(model_out$Adult)-model_out$Adult)
  return(-1*Log_LLK)
  }
 
   
  ########################################################################################################################################
  # Poisson log likelihood:: Fitting combined active and dormant mosquitoes
Posisson_Log_LKL_3Eqn_Aestivation<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
    
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=(Vector_proj[,2]+Vector_proj[,3])))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  Log_LLK<- sum(Data$Adult*log(model_out$Adult)-model_out$Adult-log(factorial(Data$Adult)))
  #Log_LLK<- sum(Data$Adult*log(model_out$Adult)-model_out$Adult)
  return(-1*Log_LLK)
  }

Gaussian_Log_LKL_3Eqn_Aestivation<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG$Time),Adult=AG$Control))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
    
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=(Vector_proj[,2]+Vector_proj[,3])))
  
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  Var=mean((Data$Adult-mean(Data$Adult))^2)
  Const=0.5*length(Data$Adult)*log(2*pi*Var)
  Gaussian_LLK<-Const+sum((Data$Adult-model_out$Adult)^2)/(2*Var)
  return((Gaussian_LLK))
 }
  
  
Log_LKL_Aestivation_3Eqn_pop<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
  
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=(Vector_proj[,2]+Vector_proj[,3])))
    
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
 
  
  Log_LLK_1<- sum(Data$Adult1*log(model_out$Adult)-model_out$Adult-log(factorial(Data$Adult1)))
  Log_LLK_2<- sum(Data$Adult2*log(model_out$Adult)-model_out$Adult-log(factorial(Data$Adult2)))
  return(-1*(Log_LLK_1+Log_LLK_2))
  }
  
  Log_LKL_Active_3Eqn_pop<- function(p_to_est, parset = names(p_to_est)) {
  p_vary[parset] <- p_to_est
  Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))
  
  t_step = 14*30.5  #number of time steps/days in a year
  V_0=Vint
  
  Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary)
  tout<-seq(0,t_step,by=h_step)
  out_V<-as.data.frame(cbind('time'=tout,'Adult'=Vector_proj[,2]))
    
  ii <- which (out_V$time %in% c(Data$time))
  model_out <- as.data.frame(cbind(time = out_V$time[ii],Adult = out_V$Adult[ii]))
  
  min_val=1.0e-7
  Log_LLK_1<- sum(Data$Adult1*log(model_out$Adult+min_val)-model_out$Adult-log(factorial(Data$Adult1)))
  Log_LLK_2<- sum(Data$Adult2*log(model_out$Adult+min_val)-model_out$Adult-log(factorial(Data$Adult2)))
  return(-1*(Log_LLK_1+Log_LLK_2))
  }
  
  
Model_DIC <- function(par_posteror_sample,lik_posterior_sample,lik_func) {
  D.bar <- 2*mean(lik_posterior_sample)
  theta.bar <- summary(par_posteror_sample)$statistics[,"Mean"]
  D.hat <- 2*lik_func(theta.bar)
  pD <- D.bar - D.hat
  pV <- var(2*lik_posterior_sample)/2
  return(list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat))
}

Model_DIC_Gaus <- function(par_posteror_sample,lik_posterior_sample,lik_func) {
  D.bar <- -mean(lik_posterior_sample)
  theta.bar <- summary(par_posteror_sample)$statistics[,"Mean"]
  D.hat <- sum(lik_func(theta.bar))
  pD <- D.bar - D.hat
  pV <- var(2*lik_posterior_sample)/2
  return(list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat))
}

Model_DIC_LSE <- function(par_posteror_sample,lik_posterior_sample,lik_func) {
  D.bar <- 2*mean(lik_posterior_sample)
  theta.bar <- summary(par_posteror_sample)$statistics[,"Mean"]
  D.hat <- -2*sum((lik_func(theta.bar))^2)
  pD <- D.bar - D.hat
  pV <- var(2*lik_posterior_sample)/2
  return(list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat))
}

  
  