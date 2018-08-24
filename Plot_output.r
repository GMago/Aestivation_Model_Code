#### G Magombedze @imperial college 2016
#############################################################################################################################
################ Plot and save estimated parameters
### 1. Save plots results H1:: Save results for 2Eqn model (Hypothesis 1)
### 2. Save plots results H2:: Save results for 2Eqn model (Hypothesis 2)
### 3. Save plots results H3:: Save results for 2Eqn model (Hypothesis 3)


Save_plots_results_H1<-function(LBq,LBq50,UBq50,UBq,p_mcmc,file_name){

tout<-seq(0,t_step,by=h_step)
Vector_pL<- Vector_2qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=LBq)
Vector_popL<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pL[,2]))

Vector_pL50<-Vector_2qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=LBq50)
Vector_popL50<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pL50[,2]))

Vector_pU50<-Vector_2qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=UBq50)
Vector_popU50<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pU50[,2]))

Vector_pU<-Vector_2qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=UBq)
Vector_popU<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pU[,2]))

Vector_pFit<-Vector_2qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pFit[,2]))

path_name<-'~/Documemts/Results/Gambiea-Poisson_Fit-95CIs_'
save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
plot(tout,Vector_popU$Non_Aquatic,type = 'n', xaxt='n', lwd=2,ylim=c(0,60),ylab='Adult Mosquito',xlab='Time in days',cex.lab=1.25,cex.axis=1.25)
axis(1, at=seq(1,14*30.5,by=(14*30.5/length(T_R$Time))), labels=T_R$Time,cex.axis=1.25)
polygon(c(rev(tout),(tout)), c(rev(Vector_popU$Non_Aquatic), (Vector_popL$Non_Aquatic)), col = "plum", border = NA)
polygon(c(rev(tout), (tout)), c(rev(Vector_popU50$Non_Aquatic), (Vector_popL50$Non_Aquatic)), col = "plum4", border = NA)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black')
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
plot(tout,Vector_popU$Non_Aquatic,type = 'n', xaxt='n', lwd=2,ylim=c(0,60),ylab='Adult Mosquito',xlab='Time in days',cex.lab=1.25,cex.axis=1.25)
axis(1, at=seq(1,14*30.5,by=(14*30.5/length(T_R$Time))), labels=T_R$Time,cex.axis=1.25)
polygon(c(rev(tout),(tout)), c(rev(Vector_popU$Non_Aquatic), (Vector_popL$Non_Aquatic)), col = "plum", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Non_Aquatic, rev(Vector_popL50$Non_Aquatic)), col = "plum4", border = NA)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black')
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
plot(tout,Vector_popU$Non_Aquatic,type = 'n', xaxt='n', lwd=2,ylim=c(0,60),ylab='Adult Mosquito',xlab='Time in days',cex.lab=1.25,cex.axis=1.25)
axis(1, at=seq(1,14*30.5,by=(14*30.5/length(T_R$Time))), labels=T_R$Time,cex.axis=1.25)
polygon(c(rev(tout),(tout)), c(rev(Vector_popU$Non_Aquatic), (Vector_popL$Non_Aquatic)), col = "plum", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Non_Aquatic, rev(Vector_popL50$Non_Aquatic)), col = "plum4", border = NA)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black')
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
dev.off()

#############################################################################
Estimated_parameters<-as.data.frame(cbind('Best Est'=p_mcmc,'UB95Ci'=UBq,'LB95Ci'=LBq,'UB50Ci'=UBq50,'LB50Ci'=LBq50))
write.csv(Estimated_parameters, paste(path_name,save_name,'_Est_pars.csv',sep=''))
}

########################################################################################################################################################
########## Save results for H2 and H3
#########################################################################################################################################################

Save_plots_results_H<-function(LBq,LBq50,UBq50,UBq,p_mcmc,file_name){

tout<-seq(0,t_step,by=h_step)
Vector_pL<- Vector_3qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=LBq)
Vector_popL<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pL[,2],'Dormant'=Vector_pL[,3]))

Vector_pL50<-Vector_3qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=LBq50)
Vector_popL50<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pL50[,2],'Dormant'=Vector_pL50[,3]))

Vector_pU50<-Vector_3qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=UBq50)
Vector_popU50<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pU50[,2],'Dormant'=Vector_pU50[,3]))

Vector_pU<-Vector_3qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=UBq)
Vector_popU<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pU[,2],'Dormant'=Vector_pU[,3]))

Vector_pFit<-Vector_3qn_seasonal_Eff(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pFit[,2],'Dormant'=Vector_pFit[,3]))

path_name<-'~/Documemts/Results/Gambiea-Poisson_Fit-95CIs_'
save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
plot(tout,Vector_popU$Active,type = 'n', xaxt='n', lwd=2,ylim=c(0,10),ylab='Adult Mosquito',xlab='Time in days',cex.lab=1.25,cex.axis=1.25)
axis(1, at=seq(1,14*30.5,by=(14*30.5/length(T_R$Time))), labels=T_R$Time,cex.axis=1.25)
polygon(c(rev(tout),tout),c(rev(Vector_popU$Dormant), Vector_popL$Dormant), col = "plum2", border = NA)
polygon(c(rev(tout),(tout)), c(rev(Vector_popU$Active), (Vector_popL$Active)), col = "lightcyan", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Active, rev(Vector_popL50$Active)), col = "skyblue", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Dormant, rev(Vector_popL50$Dormant)), col = "plum4", border = NA)
lines(tout,Vector_popFit$Active,lwd=3,col='black')
lines(tout,Vector_popFit$Dormant,lwd=3,col='red')
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
plot(tout,Vector_popU$Active,type = 'n', xaxt='n', lwd=2,ylim=c(0,10),ylab='Adult Mosquito',xlab='Time in days',cex.lab=1.25,cex.axis=1.25)
axis(1, at=seq(1,14*30.5,by=(14*30.5/length(T_R$Time))), labels=T_R$Time,cex.axis=1.25)
polygon(c(rev(tout),tout),c(rev(Vector_popU$Dormant), Vector_popL$Dormant), col = "plum2", border = NA)
polygon(c(rev(tout),(tout)), c(rev(Vector_popU$Active), (Vector_popL$Active)), col = "lightcyan", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Active, rev(Vector_popL50$Active)), col = "skyblue", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Dormant, rev(Vector_popL50$Dormant)), col = "plum4", border = NA)
lines(tout,Vector_popFit$Active,lwd=3,col='black')
lines(tout,Vector_popFit$Dormant,lwd=3,col='red')
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
plot(tout,Vector_popU$Active,type = 'n', xaxt='n', lwd=2,ylim=c(0,10),ylab='Adult Mosquito',xlab='Time in days',cex.lab=1.25,cex.axis=1.25)
axis(1, at=seq(1,14*30.5,by=(14*30.5/length(T_R$Time))), labels=T_R$Time,cex.axis=1.25)
polygon(c(rev(tout),tout),c(rev(Vector_popU$Dormant), Vector_popL$Dormant), col = "plum2", border = NA)
polygon(c(rev(tout),(tout)), c(rev(Vector_popU$Active), (Vector_popL$Active)), col = "lightcyan", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Active, rev(Vector_popL50$Active)), col = "skyblue", border = NA)
polygon(c(tout, rev(tout)), c(Vector_popU50$Dormant, rev(Vector_popL50$Dormant)), col = "plum4", border = NA)
lines(tout,Vector_popFit$Active,lwd=3,col='black')
lines(tout,Vector_popFit$Dormant,lwd=3,col='red')
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
dev.off()

#############################################################################
Estimated_parameters<-as.data.frame(cbind('Best Est'=p_mcmc,'UB95Ci'=UBq,'LB95Ci'=LBq,'UB50Ci'=UBq50,'LB50Ci'=LBq50))
write.csv(Estimated_parameters, paste(path_name,save_name,'_Est_pars.csv',sep=''))
}

#############################################################################################
# Generate CIs using FME SensRange function
Sense_cal_basic<-function(time,pars,state)
{
V_int<-c(state[1],state[2])
Output<-Model_Used(p_fix=0,h_step,V_0=Vint,p_vary=pars)
out<-as.data.frame(cbind('time'=time,'A'=Output[,2]))
return(out)
}


Sense_cal_aestivation<-function(time,pars,state)
{
V_int<-c(state[1],state[2],state[3],state[4])
Output<-Aestivation_Model_Used(p_fix=0,h_step,V_0=Vint,p_vary=pars)
out<-as.data.frame(cbind('time'=time,'A'=Output[,2],'D'=Output[,3]))
return(out)
}
###############################################################################################

Save_plots_results_basic<-function(LBq,LBq50,UBq50,UBq,p_mcmc,file_name,path_name,y_lim){
tout<-seq(0,t_step,by=h_step)
parRanges95 <- data.frame(min = LBq, max = UBq)
parRanges50 <- data.frame(min = LBq50, max = UBq50)
rownames(parRanges95)<- c('F','um','T_val')
rownames(parRanges50)<- c('F','um','T_val')
p_pass<-c(p_fix,p_mcmc)

Sense_95<- summary(sensRange(func=Sense_cal_basic,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c("A"),parRange=parRanges95,num=100))
Sense_50<- summary(sensRange(func=Sense_cal_basic,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c("A"),parRange=parRanges50,num=100))

Vector_pFit<-Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pFit[,2]))

#path_name<-'C:/Users/GMAGOMBE/Documents/Research Papers/Working progress/vector ecology/Aestivation_ODE_codes/Fitting_results/Gambiea-Poisson_Fit-95CIs_'
save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#############################################################################
Estimated_parameters<-as.data.frame(cbind('Best Est'=p_mcmc,'UB95Ci'=UBq,'LB95Ci'=LBq,'UB50Ci'=UBq50,'LB50Ci'=LBq50))
write.csv(Estimated_parameters, paste(path_name,save_name,'_Est_pars.csv',sep=''))
}

Save_plots_results_basic_CoFF<-function(LBq,LBq50,UBq50,UBq,p_mcmc,file_name,path_name,y_lim){
tout<-seq(0,t_step,by=h_step)
parRanges95 <- data.frame(min = LBq, max = UBq)
parRanges50 <- data.frame(min = LBq50, max = UBq50)
rownames(parRanges95)<- c('F','um','T_val','cutoff')
rownames(parRanges50)<- c('F','um','T_val','cutoff')
p_pass<-c(p_fix,p_mcmc)

Sense_95<- summary(sensRange(func=Sense_cal_basic,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c("A"),parRange=parRanges95,num=100))
Sense_50<- summary(sensRange(func=Sense_cal_basic,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c("A"),parRange=parRanges50,num=100))

Vector_pFit<-Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pFit[,2]))

#path_name<-'C:/Users/GMAGOMBE/Documents/Research Papers/Working progress/vector ecology/Aestivation_ODE_codes/Fitting_results/Gambiea-Poisson_Fit-95CIs_'
save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#############################################################################
Estimated_parameters<-as.data.frame(cbind('Best Est'=p_mcmc,'UB95Ci'=UBq,'LB95Ci'=LBq,'UB50Ci'=UBq50,'LB50Ci'=LBq50))
write.csv(Estimated_parameters, paste(path_name,save_name,'_Est_pars.csv',sep=''))
}
Save_plots_results_aestivation<-function(LBq,LBq50,UBq50,UBq,p_mcmc,file_name,path_name,y_lim)
{
tout<-seq(0,t_step,by=h_step)
parRanges95 <- data.frame(min = LBq, max = UBq)
parRanges50 <- data.frame(min = LBq50, max = UBq50)
rownames(parRanges95)<- c('um','umd','F','d','w','T_val')
rownames(parRanges50)<- c('um','umd','F','d','w','T_val')
p_pass<-c(p_fix,p_mcmc)

Sense_95<- summary(sensRange(func=Sense_cal_aestivation,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c('A','D'),parRange=parRanges95,num=100))
Sense_50<- summary(sensRange(func=Sense_cal_aestivation,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c('A','D'),parRange=parRanges50,num=100))

Vector_pFit<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pFit[,2],'Dormant'=Vector_pFit[,3]))

save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,which=c('A'),legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.7),ylim=c(0,y_lim),main='',cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_95,quant=F,which=c('D'),legpos=NULL,type='n',lwd=0,col=rgb(0.5,0.85,1,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
par(new=T)
plot(Sense_50,quant=F,which=c('A'),legpos=NULL,type='n',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
par(new=T)
plot(Sense_50,quant=F,which=c('D'),legpos=NULL,type='n',lwd=0,col=rgb(1.0,0.75,0.75,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
lines(tout,Vector_popFit$Active,lwd=3,col='red',ylim=c(0,y_lim))
lines(tout,Vector_popFit$Dormant,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,which=c('A'),legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.7),ylim=c(0,y_lim),main='',cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_95,quant=F,which=c('D'),legpos=NULL,type='n',lwd=0,col=rgb(0.5,0.85,1,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
par(new=T)
plot(Sense_50,quant=F,which=c('A'),legpos=NULL,type='n',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
par(new=T)
plot(Sense_50,quant=F,which=c('D'),legpos=NULL,type='n',lwd=0,col=rgb(1.0,0.75,0.75,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
lines(tout,Vector_popFit$Active,lwd=3,col='red',ylim=c(0,y_lim))
lines(tout,Vector_popFit$Dormant,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-10,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,which=c('A'),legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.7),ylim=c(0,y_lim),main='',cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_95,quant=F,which=c('D'),legpos=NULL,type='n',lwd=0,col=rgb(0.5,0.85,1,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
par(new=T)
plot(Sense_50,quant=F,which=c('A'),legpos=NULL,type='n',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
par(new=T)
plot(Sense_50,quant=F,which=c('D'),legpos=NULL,type='n',lwd=0,col=rgb(1.0,0.75,0.75,alpha=0.7),ylim=c(0,y_lim),ylab='',xlab='',main='',axes=FALSE)
lines(tout,Vector_popFit$Active,lwd=3,col='red',ylim=c(0,y_lim))
lines(tout,Vector_popFit$Dormant,lwd=3,col='black',ylim=c(0,y_lim))
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#############################################################################
Estimated_parameters<-as.data.frame(cbind('Best Est'=p_mcmc,'UB95Ci'=UBq,'LB95Ci'=LBq,'UB50Ci'=UBq50,'LB50Ci'=LBq50))
write.csv(Estimated_parameters, paste(path_name,save_name,'_Est_pars.csv',sep=''))
}

#################################################################################################
#Generate CIs using LHS package

Generate_CIs<-function(LBq,LBq50,UBq50,UBq)
{
CIs95_lhs<-cbind(min=LBq,max=UBq)
npar=500
set.seed(100)
CIs95_pars_ranges<-lhs(npar,CIs95_lhs)
CIs95_Ranges <- data.frame(CIs95_pars_ranges)
colnames(CIs95_Ranges)<- c('F','um','Kc')

CIs50_lhs<-cbind(min=LBq50,max=UBq50)
npar=500
set.seed(100)
CIs50_pars_ranges<-lhs(npar,CIs50_lhs)
CIs50_Ranges <- data.frame(CIs50_pars_ranges)
colnames(CIs50_Ranges)<- c('F','um','Kc')

p_fix<-list(k1=50,k2=0.025,ui=0.3125,p=0.07) # Fixed during fitting
V_0=c(550.0, 40.0)
t_step<-length(T_R$Rain)*30.5 
h_step=0.05
n_steps<-(t_step-0)/h_step

Sense_2qn_model<-Vector_2qn_seasonal_Eff

Sense_vector95 = matrix(0, n_steps+1,npar)
Sense_vector50 = matrix(0, n_steps+1,npar)
for (i in 1:npar)
{
Vector_proj95<-Sense_2qn_model(p_fix,h_step,V_0,p_vary=c(F=CIs95_Ranges[i,1],um=CIs95_Ranges[i,2],Kc=CIs95_Ranges[i,3]))
Vector_popltn95<-as.data.frame(cbind('Aquatic'=Vector_proj95[,1],'Adult'=Vector_proj95[,2]))
Sense_vector95[,i]<-Vector_popltn95$Adult

Vector_proj50<-Sense_2qn_model(p_fix,h_step,V_0,p_vary=c(F=CIs50_Ranges[i,1],um=CIs50_Ranges[i,2],Kc=CIs50_Ranges[i,3]))
Vector_popltn50<-as.data.frame(cbind('Aquatic'=Vector_proj50[,2],'Adult'=Vector_proj50[,2]))
Sense_vector50[,i]<-Vector_popltn50$Adult
}
col.rainbow <- rainbow(1)
matplot(tout,Sense_vector95,type='p',ylim=c(0,90),col=col.rainbow)
points((AG$Time),AG$Control,pch=16,col='darkred',cex=1.25)
}


