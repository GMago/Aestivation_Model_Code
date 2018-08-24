#### G Magombedze @imperial college 2016
#############################################################################################################################
################ Model fitting with LSE method and MCMC method
### 1. Fit data with the 2eqn model using LSE and  MCMC methds and save results
### 2. Fit data with the 3eqn model using LSE and  MCMC methds and save results

#Load the main library for model fitting
library('FME')

#Import/Load rainfall patterns/data
Rainfall_Boy_Rain <- read.csv("~/Documemts/Data/Boyila_Temp_Rain.csv")
T_R<-Rainfall_Boy_Rain 

#Load mosquito aestivation data
Mosq_data <- read.csv("~/Documemts/Data/Arabiensis_data.csv")
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
Vint= c(20.0,5) 
V_0 = Vint 

#Model parameters
p_fix<-list(ui=0.4,p=0.07,Kc=17,R_min=1.0e-12,CoFF1=3.0e-7,CoFF2=3.0e-7)
p_vary<-list(um=0.098,F=2.55,T_val=80) #To vary during fitting
#Set simulation end point and iteration step size 

t_step<-length(T_R$Rain)*30.5 
h_step=0.05

#Carry model simulation and view population trajectories
# Select either the Null model
#[1.] Null model: Vector_2qn_seasonal_NULL

Basic_Model_Used<-Vector_2qn_seasonal_NULL

Vector_proj<-Basic_Model_Used(p_fix,h_step,V_0,p_vary) #run simulations
tout<-seq(0,t_step,by=h_step) #Time
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2])) #Format output
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',lwd=2,col='blue') #plot output
points(AG$Time,AG$Control,pch=16,col='red')
#points(seq(30.5,(length(T_R$Time))*30.5,by=30.5),T_R$Rain*40/max(T_R$Rain),col='blue',type='o')

###############################################################################################################
# LSE-Model fitting using the no-aeastivation model:: H1 

#F,um,Kc upper and lower parameter bounce
LB<-c(0.01,0.01,70)
UB<-c(20.0,0.20,100) 

#p_start<-c(F=20,um=0.05,T_val=20,cutoff=1.0e-30)              #Intial parameters guesses
p_start<-c(F=2.55,um=0.098,T_val=80)              #Intial parameters guesses

Model_Used<-Vector_2qn_seasonal_NULL       #Fit data with the using the NULL model
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
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',ylim = c(0,8),lwd=2,col='red',xlab='Time',ylab='Adult Mosquitoes')
points(AG$Time,AG$Control,pch=16,col='red')


SF<-summary(Fit_1)
Var0 <- SF$modVariance
covIni <- SF$cov.scaled *2.4^2/length(p_est)

#MCMC <- modMCMC(p=p_est,f=LSE_2qn_residuals,niter=15000,lower=LB,upper=UB, var0=Var0,wvar0 = 1,updatecov=100,ntrydr=2,burninlength=5000,jump=p_est*0.05)
#MCMC variables:: ntrydr=delayed rejection, jump=proposal covariance matrix or proposal distribution, updatecov=steps for covariance update, var0=NULL (log likelihood is used), wvar0=0(prior is ignored, if 1 equal weights for the prior)
MCMC <- modMCMC(p=p_est,f=Posisson_Log_LKL_2Eqn,niter=5000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=100,ntrydr=2,burninlength=1000,jump=covIni,prior=NULL)

MCMC <- modMCMC(p=p_est,f=Posisson_Log_LKL_2Eqn,niter=5000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=50,ntrydr=5,jump=covIni)

MCMC <- modMCMC(p=p_mcmc,f=Posisson_Log_LKL_2Eqn,niter=11000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=30,jump=p_mcmc*0.25,burninlength = 1000)

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

path_name<-'~/Documemts/Results/Arabiensis_'
Save_plots_results_basic(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H0',path_name,y_lim=10)

path_name<-'~/Documemts/Results/MCMC_'
save_name<-'H0'
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

