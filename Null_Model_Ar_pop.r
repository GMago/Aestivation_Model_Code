#### G Magombedze @imperial college 2016
#############################################################################################################################
################ Model fitting with LSE method and MCMC method
### Fit data with the 2eqn model using LSE and  MCMC methds and save results [H1]

#Load the main library for model fitting
library('FME')

#Import/Load rainfall patterns/data
Rainfall_Boy_Rain <- read.csv("~/Documemts/Data/Boyila_Temp_Rain.csv")
T_R<-Rainfall_Boy_Rain 

#load mosquito village data into a data frame
Mosq_dataV1 <- read.csv("~/Documemts/Data/ArabiensisV1_data.csv")
Mosq_dataV2 <- read.csv("~/Documemts/Data/ArabiensisV2_data.csv")

AG_Villages<-as.data.frame(cbind('Time'=round(Mosq_dataV1$Time),'ControlV1'=Mosq_dataV1$Control,'ControlV2'=Mosq_dataV2$Control))
Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))

#Load model functions:: 1.Rainfall fitness, 2. Adaptation selection, 3. Reactivation selection, 4. Rainfall averaging, 5. 2Eqn model, 6. 3Eqn model,and 7. 4Eqn model
source("~/Documemts/Code/Model_Functions.r")

#Load model fitting likelihood functions:: 1.LSE/gaussina LLK, 2. Poisson LLK (tailored for each separate model)
source("~/Documemts/Code/Likelihood_Functions.r")

#########################################################################################################################################################
#########################################################################################################################################################
##### [1.] Fit data with the 2eqn model using LSE and  MCMC methds and save results
#########################################################################################################################################################

Vint= c(20.0,5) 
V_0 = Vint 

p_fix<-list(k1=80,k2=0.5,ui=0.4,p=0.07,Kc=17.0) # Fixed during fitting
#p_vary<-list(um=0.08,F=13,Kc=20.0) #To vary during fitting
p_vary<-list(um=0.1,F=10) #To vary during fitting
#Set simulation end point and iteration step size 
t_step<-length(T_R$Rain)*30.5 
h_step=0.05

#Carry model simulation and view population trajectories
# Select either the H1 model or the Seasona Effects model
#[H1.] Seasonal Effects model: Vector_2qn_seasonal_Eff

Basic_Model_Used<-Vector_2qn_seasonal_NULL 

Vector_proj<-Basic_Model_Used(p_fix,h_step,V_0,p_vary) #run simulations
tout<-seq(0,t_step,by=h_step) #Time
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2])) #Format output
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',lwd=2,col='blue') #plot output
points(Data$time,Data$Adult1,pch=16,col='red')
points(Data$time,Data$Adult2,pch=16,col='green')

###############################################################################################################
# LSE-Model fitting using the no-aeastivation model:: H1 

#F,um,Kc upper and lower parameter bounce
LB<-c(5.00,0.01)
UB<-c(20.0,0.20) 

p_start<-c(F=10,um=0.1)              #Intial parameters guesses
Model_Used<-Vector_2qn_seasonal_NULL  

MCMC <- modMCMC(p=p_start,f=Log_LKL_poisson_2Eqn_Pop,niter=1100,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=100,ntrydr=2,burninlength=100,jump=p_start*0.1)

q_par<-summary(as.mcmc(MCMC$pars))

p_mcmc<-q_par$quantiles[,3]

Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
tout<-seq(0,t_step,by=h_step)
rainfall<-Interpolate_data(T_R$Rain,h_step)
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Aquatic'=Vector_proj[,1],'Adult'=Vector_proj[,2]))
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
par(new = TRUE)
plot(Vector_popltn$Time,Vector_popltn$Adult,type='l',ylim = c(0,10),lwd=2,col='blue',xlab='Time',ylab='Adult Mosquitoes')
points(Data$time,Data$Adult1,pch=16,col='red')
points(Data$time,Data$Adult2,pch=16,col='green')

MCMC <- modMCMC(p=p_start,f=Log_LKL_poisson_2Eqn_Pop,niter=5000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=100,ntrydr=5,burninlength=1000,jump=p_start*0.025)

q_par<-summary(as.mcmc(MCMC$pars))
#Extract the median
p_mcmc<-q_par$quantiles[,3]

# Extract 95% CI bounds
UBq<-q_par$quantiles[1:length(p_mcmc),5]
LBq<-q_par$quantiles[1:length(p_mcmc),1]

# Extract 50% CI bounds
UBq50<-q_par$quantiles[1:length(p_mcmc),4]
LBq50<-q_par$quantiles[1:length(p_mcmc),2]

#h_step=0.01
rainfall<-Interpolate_data(T_R$Rain,h_step)
source("~/Documemts/Code/Plot_output_pop.r")
# Save results:: Using the save plots function in Plot_outout.r script

path_name<-'~/Documemts/Results/Arabiensis_pop'
Save_plots_results_basic(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H0_Ar',path_name,y_lim=10)

#Save chain and distributions
path_name<-'~/Documemts/Results/MCMC_pop'
save_name<-'H0_Ar'
mcmc.chain<-as.mcmc(MCMC$pars)
write.csv(mcmc.chain, paste(path_name,save_name,'_chain.csv',sep=''))

pdf(paste(path_name,save_name,'.pdf',sep=''))
plot(mcmc.chain,smooth=T)
dev.off()

#Calculate model DIC

par_posteror_sample<-as.mcmc(MCMC$pars)
lik_posterior_sample<-MCMC$SS
lik_func<-Log_LKL_poisson_2Eqn_Pop
K_num<-length(p_mcmc)
D_num<-length(Data$Adult2l)

DIC<-Model_DIC_Gaus(par_posteror_sample,lik_posterior_sample,lik_func)
write.csv(DIC, paste(path_name,save_name,'_DIC.csv',sep=''))
