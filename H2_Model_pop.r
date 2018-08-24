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
Mosq_dataV1 <- read.csv("~/Documemts/Data/AgambiaeV1_Total.csv")
Mosq_dataV2 <- read.csv("~/Documemts/Data/AgambiaeV2_Total.csv")

AG_Villages<-as.data.frame(cbind('Time'=round(Mosq_dataV1$Time),'ControlV1'=Mosq_dataV1$ControlM,'ControlV2'=Mosq_dataV2$ControlM))
Data<-as.data.frame(cbind(time=(AG_Villages$Time),Adult1=AG_Villages$ControlV1,Adult2=AG_Villages$ControlV2))

#Load model functions:: 1.Rainfall fitness, 2. Adaptation selection, 3. Reactivation selection, 4. Rainfall averaging, 5. 2Eqn model, 6. 3Eqn model,and 7. 4Eqn model
source("~/Documemts/Code/Model_Functions.r")

#Load model fitting likelihood functions:: 1.LSE/gaussina LLK, 2. Poisson LLK (tailored for each separate model)
source("~/Documemts/Code/Likelihood_Functions.r")

#########################################################################################################################################################
#########################################################################################################################################################
##### [1.] Fit data with the 2eqn model using LSE and  MCMC methds and save results
#########################################################################################################################################################

# Model initial conditions
# Model initial conditions
Vint= c(100,40,0,15) 
V_0 = Vint 

#Model parameters
p_fix<-list(ui=0.4,p=0.07,Kc=17.0,R_min=1.0e-12,CoFF1=3.0e-7,CoFF2=3.0e-7) # Fixed during fitting
p_vary<-list(um=0.098,umd=0.01,F=13,d=0.02,w=1.0,T_val=80) #To vary during fitting

#Set simulation end point and iteration step size 
t_step<-length(T_R$Rain)*30.5 
h_step=0.05

#Carry model simulation and view population trajectories
#Select model Vector_3qn_seasonal_Eff

Model_Used<-Vector_3qn_seasonal_Eff

Vector_proj<-Model_Used(p_fix,h_step,V_0,p_vary) #run simulations
tout<-seq(0,t_step,by=h_step) #Time
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Active'=Vector_proj[,2],'Dormant'=Vector_proj[,3])) #Format output
plot(Vector_popltn$Time,Vector_popltn$Active,type='l',ylim=c(0,60),lwd=2,col='red',xlab='Time',ylab='Adult Mosquitoes') #plot output
lines(Vector_popltn$Time,Vector_popltn$Dormant,type='l',lwd=2,col='black')
points(Data$time,Data$Adult1,pch=16,col='red')
points(Data$time,Data$Adult2,pch=16,col='green')

###############################################################################################################
# LSE-Model fitting using the no-aeastivation model:: H1 

#F,um,Kc upper and lower parameter bounce
LB<-c(0.05,0.005,01.00,0.01,0.50,70)
UB<-c(0.10,0.025,20.00,0.10,1.00,100)

p_start<-c(um=0.098,umd=0.01,F=13,d=0.02,w=1,T_val=80) #Intial parameters guesses

Aestivation_Model_Used<-Vector_3qn_seasonal_Eff

MCMC <- modMCMC(p=p_start,f=Log_LKL_Active_3Eqn_pop,niter=110000,lower=LB,upper=UB,wvar0 = 1,var0=NULL,updatecov=100,ntrydr=5,burninlength=10000,jump=p_start*0.1)

q_par<-summary(as.mcmc(MCMC$pars))

p_mcmc<-q_par$quantiles[,3]

Vector_proj<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
tout<-seq(0,t_step,by=h_step) #Time
rainfall<-Interpolate_data(T_R$Rain,h_step) #rainfall
Vector_popltn<-as.data.frame(cbind('Time'=tout,'Active'=Vector_proj[,2],'Dormant'=Vector_proj[,3])) #Format output
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
par(new = TRUE)
plot(Vector_popltn$Time,Vector_popltn$Active,type='l',ylim=c(0,60),lwd=2,col='red',xlab='Time',ylab='Adult Mosquitoes') #plot output
lines(Vector_popltn$Time,Vector_popltn$Dormant,type='l',lwd=2,col='black')
#lines(Vector_popltn$Time,Vector_popltn$Dormant+Vector_popltn$Active,type='l',lwd=2,col='blue')
points(Data$time,Data$Adult1,pch=16,col='red')
points(Data$time,Data$Adult2,pch=16,col='green')


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
source("~/Documemts/Code/Plot_output_pop.r")
path_name<-'~/Documemts/Results/Gambiea_pop_'
Save_plots_results_aestivation(LBq,LBq50,UBq50,UBq,p_mcmc,file_name='H2x',path_name,y_lim=60)

#Save chain and distributions
path_name<-'~/Documemts/Results/MCMC_pop_'
save_name<-'H2'
mcmc.chain<-as.mcmc(MCMC$pars)
write.csv(mcmc.chain, paste(path_name,save_name,'_chain.csv',sep=''))

pdf(paste(path_name,save_name,'.pdf',sep=''))
plot(mcmc.chain,smooth=T)
dev.off()

par_posteror_sample<-as.mcmc(MCMC$pars)
lik_posterior_sample<-MCMC$SS
lik_func<-Log_LKL_Active_3Eqn_pop
K_num<-length(p_mcmc)
D_num<-length(Data$Adult2)

DIC<-Model_DIC(par_posteror_sample,lik_posterior_sample,lik_func)
write.csv(DIC, paste(path_name,save_name,'_DIC.csv',sep=''))
