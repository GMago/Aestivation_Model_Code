#### G Magombedze @imperial college 2016
#############################################################################################################################
################ Plot and save estimated parameters
### 1. Save plots results H1:: Save results for 2Eqn model (Hypothesis 1)
### 2. Save plots results H2:: Save results for 2Eqn model (Hypothesis 2)
### 3. Save plots results H3:: Save results for 2Eqn model (Hypothesis 3)


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
rownames(parRanges95)<- c('F','um')
rownames(parRanges50)<- c('F','um')
p_pass<-c(p_fix,p_mcmc)

Sense_95<- summary(sensRange(func=Sense_cal_basic,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c("A"),parRange=parRanges95,num=100))
Sense_50<- summary(sensRange(func=Sense_cal_basic,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c("A"),parRange=parRanges50,num=100))

Vector_pFit<-Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Non_Aquatic'=Vector_pFit[,2]))

save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points(AG$Time,AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points(AG$Time,AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
axis(side=4, at = pretty(range(rainfall)),cex=1.25,cex.axis=1.25,cex.lab=1.5)
mtext("Rainfall (mm/month)",side=4,line=3,cex.lab=1.5,cex=1.5)
par(new = TRUE)
plot(Sense_95,quant=F,legpos=NULL,type='n',ylab='Adult Mosquito Population',xlab='Time in days',main='',lwd=0,col=rgb(0.5,0.5,0.5,alpha=0.70),ylim=c(0,y_lim),cex=1.25,cex.axis=1.25,cex.lab=1.5)
par(new=T)
plot(Sense_50,quant=F,legpos=NULL,type='n',ylab='',xlab='',main='',lwd=0,col=rgb(1.0,0.5,0.75,alpha=0.7),ylim=c(0,y_lim),axes=FALSE)
lines(tout,Vector_popFit$Non_Aquatic,lwd=3,col='black',ylim=c(0,y_lim))
points(AG$Time,AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
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
#rownames(parRanges95)<- c('um','F')
#rownames(parRanges50)<- c('um','F')
rownames(parRanges95)<- c('um','umd','F','d')
rownames(parRanges50)<- c('um','umd','F','d')
p_pass<-c(p_fix,p_mcmc)

Sense_95<- summary(sensRange(func=Sense_cal_aestivation,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c('A','D'),parRange=parRanges95,num=100))
Sense_50<- summary(sensRange(func=Sense_cal_aestivation,parms=c(p_pass),state=V_0,time=tout,dist="latin",sensvar=c('A','D'),parRange=parRanges50,num=100))

Vector_pFit<-Aestivation_Model_Used(p_fix,h_step,V_0,p_vary=p_mcmc)
Vector_popFit<-as.data.frame(cbind('Time'=tout,'Active'=Vector_pFit[,2],'Dormant'=Vector_pFit[,3]))

save_name<-file_name

#save png
png(paste(path_name,save_name,'.png',sep=''),width=7*300,height=5*300,res=300)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
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
points(AG$Time,AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save pdf
pdf(paste(path_name,save_name,'.pdf',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
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
points(AG$Time,AG$Control,pch=16,col='red',cex=1.25,ylim=c(0,y_lim))
dev.off()

#save svg
svg(paste(path_name,save_name,'.svg',sep=''),width=8,height=5,pointsize=15)
par(mar=c(5,4,4,5)+.1)
barplot(rainfall[seq(1,t_step/h_step,by=30.5/h_step)],col=rgb(0.0,0,1,alpha=0.6),lwd=3,ylim = c(-4,max(rainfall)+25),axes = FALSE,bty = "n", xlab = "", ylab = "") 
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
points(AG$Time,AG$Control,pch=16,col='darkred',cex=1.25,ylim=c(0,y_lim))
dev.off()

#############################################################################
Estimated_parameters<-as.data.frame(cbind('Best Est'=p_mcmc,'UB95Ci'=UBq,'LB95Ci'=LBq,'UB50Ci'=UBq50,'LB50Ci'=LBq50))
write.csv(Estimated_parameters, paste(path_name,save_name,'_Est_pars.csv',sep=''))
}

