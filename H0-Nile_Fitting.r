#### G Magombedze @imperial college 2016
#############################################################################################################################

#Load the main library for model fitting
library('FME')

#Import/Load rainfall patterns/data
Mosq_data <- read.csv("~/Documemts/Data/Adult_Female_Dry-Wet.csv")

Rain_Sudan <- read.csv("~/Documemts/Data/Mean-dairly-data-Khortum-district.csv")
T_R<-as.data.frame(cbind('Time'=(Rain_Sudan$Time),'Rain'=Rain_Sudan$Rainfall+25,'Humidity'=rowMeans(Rain_Sudan[,5:7])))

#DataF<-as.data.frame(cbind('Time'=(Mosq_data$Time),'Adult1'=Mosq_data$Dry,'Adult2'=Mosq_data$Wet))
DataF<-as.data.frame(cbind('Time'=(Mosq_data$Time),'Adult1'=Mosq_data$Fattasha,'Adult2'=Mosq_data$Nile))

AG<-as.data.frame(cbind('Time'=DataF$Time,'Control'=DataF$Adult2))

source("~/Documemts/Code/Model_Functions.r")

#Load model fitting likelihood functions:: 1.LSE/gaussina LLK, 2. Poisson LLK (tailored for each separate model)
source("~/Documemts/Code/Likelihood_Functions.r")

Vint= c(0,DataF$Adult2[1]) 
V_0 = Vint 

#Model parameters
p_fix<-list(k1=25,k2=0.5,ui=0.4,p=0.07,Kc=17.0) # Fixed during fitting
#p_vary<-list(um=0.08,F=13,Kc=20.0) #To vary during fitting
p_vary<-list(um=0.1,F=10) #To vary during fitting
#Set simulation end point and iteration step size 
t_step<-length(T_R$Rain)*30.5 
h_step=0.05

#Carry model simulation and view population trajectories
Basic_Model_Used<-Vector_2qn_seasonal_NULL 

Vector_proj<-Basic_Model_Used(p_fix,h_step,V_0,p_vary) #run simulations
tout<-seq(0,t_step,by=h_step) #Time
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2])) #Format output
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',lwd=2,col='blue',ylim=c(0,120),xlab='Time',ylab='Adult Mosquitoes') #plot output
points(DataF$Time,DataF$Adult2,pch=16,col='red')

###############################################################################################################
LB<-c(1.00,0.01)
UB<-c(20.0,0.20) 

p_start<-c(F=10,um=0.1)              #Intial parameters guesses
Model_Used<-Vector_2qn_seasonal_NULL  

Fit_1 <- modFit(f=LSE_2qn_residuals,p=p_start,lower=LB, upper=UB) #LSE fitting
p_est=coef(Fit_1)

Fit_1 <- modFit(f=LSE_2qn_residuals,p=p_est,lower=LB, upper=UB) #LSE fitting
p_est=coef(Fit_1)

# Visualise fitting results
Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary=p_est)
tout<-seq(0,t_step,by=h_step) #Time
rainfall<-Interpolate_data(T_R$Rain,h_step) #rainfall
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2]))
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
par(new = TRUE)
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',ylim = c(0,120),lwd=2,col='blue',xlab='Time',ylab='Adult Mosquitoes')
points(DataF$Time,DataF$Adult2,pch=16,col='red')

SF<-summary(Fit_1)
Var0 <- SF$modVariance
covIni <- SF$cov.scaled *2.4^2/length(p_est)

#MCMC <- modMCMC(p=p_est,f=Posisson_Log_LKL_2Eqn,niter=5000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=50,jump=covIni)
MCMC <- modMCMC(p=p_est,f=Posisson_Log_LKL_2Eqn,niter=100000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=20,burninlength=10000,jump=p_est*0.25)

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
Save_plots_results_basic(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H0_wet_All_Nile',path_name,y_lim=120)

#Save chain and distributions
path_name<-'~/Documemts/Results/MCMC_'
save_name<-'H0_wet_All_Nile'
mcmc.chain<-as.mcmc(MCMC$pars)
write.csv(mcmc.chain, paste(path_name,save_name,'_chain.csv',sep=''))

pdf(paste(path_name,save_name,'.pdf',sep=''))
plot(mcmc.chain,smooth=T)
dev.off()

par_posteror_sample<-as.mcmc(MCMC$pars)
lik_posterior_sample<-MCMC$SS
lik_func<-Posisson_Log_LKL_2Eqn
K_num<-length(p_mcmc)
D_num<-length(DataF$Adult2)

DIC<-Model_DIC(par_posteror_sample,lik_posterior_sample,lik_func)
write.csv(DIC, paste(path_name,save_name,'_DIC.csv',sep=''))
