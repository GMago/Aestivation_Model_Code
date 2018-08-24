#### G Magombedze @imperial college 2016
#############################################################################################################################

#Load the main library for model fitting
library('FME')

#Import/Load rainfall patterns/data
#Mosq_data <- read.csv("C:/Users/Gesham/Desktop/File Folders/Extracted data/Sudan_data/Adult_Female_vectors.csv")
Mosq_data <- read.csv("~/Documemts/Data/Adult_Female_Dry-Wet.csv")
Rain_Sudan <- read.csv("~/Documemts/Data/Mean-dairly-data-Khortum-district.csv")
T_R<-as.data.frame(cbind('Time'=(Rain_Sudan$Time),'Rain'=Rain_Sudan$Rainfall+25,'Humidity'=rowMeans(Rain_Sudan[,5:7])))

#DataF<-as.data.frame(cbind('Time'=(Mosq_data$Time),'Adult1'=Mosq_data$Dry,'Adult2'=Mosq_data$Wet))
DataF<-as.data.frame(cbind('Time'=(Mosq_data$Time),'Adult1'=Mosq_data$Fattasha,'Adult2'=Mosq_data$Nile))

AG<-as.data.frame(cbind('Time'=DataF$Time,'Control'=DataF$Adult2))

source("~/Documemts/Code/Model_Functions.r")

#Load model fitting likelihood functions:: 1.LSE/gaussina LLK, 2. Poisson LLK (tailored for each separate model)
source("~/Documemts/Code/Likelihood_Functions.r")

Vint= c(DataF$Adult2[1],DataF$Adult2[1],0,2) 
V_0 = Vint 

#Model parameters
p_fix<-list(k1=25,k2=0.5,ui=0.3,p=0.07,Fd=0,Kc=17,w=1) # Fixed during fitting
p_vary<-list(um=0.05,F=10,umd=0.01,d=0.08) #To vary during fitting

#Set simulation end point and iteration step size 
t_step<-length(T_R$Rain)*30.5 
h_step=0.05

#Carry model simulation and view population trajectories
#Select model Vector_3qn_seasonal_Eff
Model_Used<-Vector_3qn_seasonal_Eff

Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary) #run simulations
tout<-seq(0,t_step,by=h_step) #Time
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Active'=Vector_proj[,2],'Dormant'=Vector_proj[,3])) #Format output
plot(Vector_popltn$Time,Vector_popltn$Active,type='l',ylim=c(0,90),lwd=2,col='red',xlab='Time',ylab='Adult Mosquitoes') #plot output
lines(Vector_popltn$Time,Vector_popltn$Dormant,type='l',lwd=2,col='black')
points(DataF$Time,DataF$Adult2,pch=16,col='red')

###############################################################################################################
# LSE-Model fitting using the no-aeastivation model:: H1 

#F,um,Kc upper and lower parameter bounce

LB<-c(0.05,0.0100,01.00,0.0800)
UB<-c(0.10,0.0110,20.00,0.0810)


Aestivation_Model_Used<-Vector_3qn_seasonal_Eff
p_start<-c(um=0.098,umd=0.01,F=13,d=0.08) #Intial parameters guesses

Fit_1 <- modFit(f=LSE_3qn_Residuals_active_and_dormant_combined,p=p_start,lower=LB, upper=UB) #LSE fitting
#Fit_1 <- modFit(f=LSE_3qn_residuals_active,p=p_start,lower=LB, upper=UB) #LSE fitting
p_est=coef(Fit_1)

# Visualise fitting results
Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary=p_est)
tout<-seq(0,t_step,by=h_step) #Time
tout<-seq(0,t_step,by=h_step) #Time
rainfall<-Interpolate_data(T_R$Rain,h_step) #rainfall
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Active'=Vector_proj[,2],'Dormant'=Vector_proj[,3])) #Format output
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
par(new = TRUE)
plot(Vector_popltn$Time,Vector_popltn$Active,type='l',ylim=c(0,120),lwd=2,col='red',xlab='Time',ylab='Adult Mosquitoes') #plot output
lines(Vector_popltn$Time,Vector_popltn$Dormant,type='l',lwd=2,col='black')
points(DataF$Time,DataF$Adult2,pch=16,col='red')


SF<-summary(Fit_1)
Var0 <- SF$modVariance
covIni <- SF$cov.scaled *2.4^2/length(p_est)

#MCMC <- modMCMC(p=p_est,f=LSE_3qn_Residuals_active_and_dormant_combined,niter=2000,lower=LB,upper=UB, var0=Var0,wvar0 = 1,updatecov=100,ntrydr=5,burninlength=1000,jump=covIni)
#
MCMC <- modMCMC(p=p_est,f=Posisson_Log_LKL_3Eqn_Aestivation,niter=110000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=25,ntrydr=3,burninlength=10000,jump=p_est*0.1)

q_par<-summary(as.mcmc(MCMC$pars))

#Extract the median
p_mcmc<-q_par$quantiles[,3]

# Extract 95% CI bounds
UBq<-q_par$quantiles[1:length(p_mcmc),5]
LBq<-q_par$quantiles[1:length(p_mcmc),1]

# Extract 50% CI bounds
UBq50<-q_par$quantiles[1:length(p_mcmc),4]
LBq50<-q_par$quantiles[1:length(p_mcmc),2]

# Save results:: Using the save plots function in Plot_outout.r script
h_step=0.05
rainfall<-Interpolate_data(T_R$Rain,h_step)
source("~/Documemts/Code/Plot_output_Sudan.r")
path_name<-'~/Documemts/Results/Nile_'
#Save_plots_results_basic(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H2_Nile2p',path_name,y_lim=120)
Save_plots_results_aestivation(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H2_wet_All2_Nile',path_name,y_lim=120)

#Save chain and distributions
path_name<-'~/Documemts/Results/MCMC_'
save_name<-'H2_wet_All2_Nile'
mcmc.chain<-as.mcmc(MCMC$pars)
write.csv(mcmc.chain, paste(path_name,save_name,'_chain.csv',sep=''))

pdf(paste(path_name,save_name,'.pdf',sep=''))
plot(mcmc.chain,smooth=T)
dev.off()

par_posteror_sample<-as.mcmc(MCMC$pars)
lik_posterior_sample<-MCMC$SS
#lik_func<-Posisson_Log_LKL_3Eqn_active
lik_func<-Posisson_Log_LKL_3Eqn_Aestivation
K_num<-length(p_mcmc)
D_num<-length(DataF$Adult2)

DIC<-Model_DIC(par_posteror_sample,lik_posterior_sample,lik_func)
write.csv(DIC, paste(path_name,save_name,'_DIC.csv',sep=''))
