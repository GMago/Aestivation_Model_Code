#### G Magombedze @imperial college 2016
#############################################################################################################################
################ Model fitting with LSE method and MCMC method
### 1. Fit data with the 2eqn model using LSE and  MCMC methds and save results
### 2. Fit data with the 3eqn model using LSE and  MCMC methds and save results

#Load the main library for model fitting
library('FME')
#library('compiler')
#Import/Load rainfall patterns/data

Rainfall_Boy_Rain <- read.csv("~/Documemts/Data/Boyila_Temp_Rain.csv")

T_R<-Rainfall_Boy_Rain 

#Load mosquito aestivation data
Mosq_data <- read.csv("~/Documemts/Data/AgambiaeSI_data.csv")

AG<-as.data.frame(cbind('Time'=round(Mosq_data$Time),'Control'=Mosq_data$Control,'Treated'=Mosq_data$Treated))
#AG<-as.data.frame(cbind('Time'=round(Mosq_data$Time-Mosq_data$Time[1]),'Control'=Mosq_data$Control,'Treated'=Mosq_data$Treated))

#Load model functions:: 1.Rainfall fitness, 2. Adaptation selection, 3. Reactivation selection, 4. Rainfall averaging, 5. 2Eqn model, 6. 3Eqn model,and 7. 4Eqn model
source("~/Documemts/Code/Model_Functions.r")

#Load model fitting likelihood functions:: 1.LSE/gaussina LLK, 2. Poisson LLK (tailored for each separate model)
source("~/Documemts/Code/Likelihood_Functions.r")

#########################################################################################################################################################
#########################################################################################################################################################
##### [1.] Fit data with the 2eqn model using LSE and  MCMC methds and save results
#########################################################################################################################################################

# Model initial conditions
Vint= c(100.0,40) 
#Vint= c(0.0,0.0) 

V_0 = Vint 

#Model parameters
#p_fix<-list(ui=0.4,p=0.07,Kc=17,cutoff=1.0e-34) # Fixed during fitting
#p_vary<-list(um=0.08,F=13,Kc=20.0) #To vary during fitting
#p_fix<-list(ui=0.4,p=0.07,Kc=17,R_min=1.0e-12,CoFF1=1.0e-7,CoFF2=1.0e-7)
p_fix<-list(ui=0.4,p=0.07,Kc=17,R_min=1.0e-12,CoFF1=1.0e-2,CoFF2=1.0e-1)

#p_vary<-list(um=0.07,F=5,T_val=80) #To vary during fitting
p_vary<-list(um=0.046176,F=2.62,T_val=80) #To vary during fitting
#Set simulation end point and iteration step size 
t_step<-length(T_R$Rain)*30.5 
h_step=0.05

#Carry model simulation and view population trajectories
# Select either the Null model
#[1.] Null model: Vector_2qn_seasonal_NULL

Basic_Model_Used<-Vector_2qn_seasonal_NULL
#Basic_Model_Used<-Vector_2qn_seasonal_NULL2

#Basic_Model_Used<-Vector_2qn_seasonal_Robust
Vector_proj<-Basic_Model_Used(p_fix,h_step,V_0,p_vary) #run simulations
tout<-seq(0,t_step,by=h_step) #Time
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2])) #Format output
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',lwd=2,col='blue') #plot output
points(AG$Time,AG$Control,pch=16,col='red')

FM<-Vector_popltn$Adult*p_vary$F
pI<-Vector_popltn$Aquatic*p_fix$p*0.5
plot(Vector_popltn$Time,log10(FM+1),pch=16,col='magenta',type='l')
lines(Vector_popltn$Time,log10(pM+1),pch=16,col='black',type='l')
plot(Vector_popltn$Time,pM,pch=16,col='black',type='l')

###############################################################################################################
# LSE-Model fitting using the no-aeastivation model:: H1 

#F,um,Kc upper and lower parameter bounce
#LB<-c(5.00,0.01,0,1.0e-35)
#UB<-c(20.0,0.200,100,1.0e-6) 
LB<-c(1.00,0.01,70)
UB<-c(20.0,0.20,100) 

#p_start<-c(F=20,um=0.05,T_val=20,cutoff=1.0e-30)              #Intial parameters guesses
p_start<-c(F=2.62,um=0.046176,T_val=80)              #Intial parameters guesses

Model_Used<-Vector_2qn_seasonal_NULL2

Fit_1 <- modFit(f=LSE_2qn_residuals,p =p_start,lower=LB, upper=UB) #LSE fitting

p_est=coef(Fit_1)

# Visualise fitting results
Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary=p_est)
tout<-seq(0,t_step,by=h_step)
rainfall<-Interpolate_data(T_R$Rain,h_step)
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2]))
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
par(new = TRUE)
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',ylim = c(0,60),lwd=2,col='red',xlab='Time',ylab='Adult Mosquitoes')
points(AG$Time,AG$Control,pch=16,col='red')


SF<-summary(Fit_1)
Var0 <- SF$modVariance
covIni <- SF$cov.scaled *2.4^2/length(p_est)

MCMC <- modMCMC(p=p_est,f=Posisson_Log_LKL_2Eqn,niter=110000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=50,jump=p_est*0.25)

q_par<-summary(as.mcmc(MCMC$pars))

#Extract the median
p_mcmc<-q_par$quantiles[,3]

# Extract 95% CI bounds
UBq<-q_par$quantiles[1:length(p_mcmc),5]
LBq<-q_par$quantiles[1:length(p_mcmc),1]

# Extract 50% CI bounds
UBq50<-q_par$quantiles[1:length(p_mcmc),4]
LBq50<-q_par$quantiles[1:length(p_mcmc),2]

h_step=0.05
rainfall<-Interpolate_data(T_R$Rain,h_step)
source("~/Documemts/Code/Plot_output.r")
# Save results:: Using the save plots function in Plot_outout.r script

path_name<-'~/Documemts/Code/Gambiea_'
Save_plots_results_basic(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H0_',path_name,y_lim=60)

path_name<-'~/Documemts/Results/MCMC_'
save_name<-'H0_'
mcmc.chain<-as.mcmc(MCMC$pars)
write.csv(mcmc.chain, paste(path_name,save_name,'_chain.csv',sep=''))

pdf(paste(path_name,save_name,'.pdf',sep=''))
plot(mcmc.chain,smooth=T)
dev.off()

par_posteror_sample<-as.mcmc(MCMC$pars)
lik_posterior_sample<-MCMC$SS
lik_func<-Posisson_Log_LKL_2Eqn
K_num<-length(p_mcmc)
D_num<-length(AG$Control)

DIC<-Model_DIC(par_posteror_sample,lik_posterior_sample,lik_func)
write.csv(DIC, paste(path_name,save_name,'_DIC.csv',sep=''))

