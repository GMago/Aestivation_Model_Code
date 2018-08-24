#### G Magombedze @imperial college 2015
#############################################################################################################################
################ Functions for the mosquito aestivation and no-aestivation models ###########################################
### 1. Rainfall fitness function, based on rainfall history
### 2. Rainfall mosquito reactivation fitness
### 3. Temperature fitness function [not used]
### 4. Adaptation fitness function [mosquito adapatation to shortage of rainfall]
### 5. Two (aquatic and non-aquatic) equation model, ecological mosquito population model (no aestivation model), one phenotype population model
### 6. Four equation model (active and dormanat species, aquatic and non-aquatic). Two phenotype species model (aestivation model). 

############################################################################################################################################################
#Rainfall fitness function. Logistic hazard using 7 day moving average
Fitness_rain <-function(rainfall,moisture_retention,h_step,Rain_Threshold)
{
	di<-moisture_retention*(1/h_step)
	length_period<-length(rainfall)
	Fitness<-numeric(length(rainfall))
	rain<-numeric(length(rainfall))
	Mo=Rain_Threshold
		for (day in 1:length_period)
		{
			if(day<=di)
			{rain[day]<-mean(rainfall[1:day])}
			else
			{rain[day]<-mean(rainfall[((day-di)+1):day])}
		Fitness[day]<-1/(1+exp(Mo-rain[day]))
		}
return(Fitness)
}

# For calculating the carrying capacity
#Moving rainfall average
Average_rain <-function(rainfall,moisture_retention,h_step)
{
	di<-moisture_retention*(1/h_step)
	length_period<-length(rainfall)
	Average<-numeric(length(rainfall))
	rain<-numeric(length(rainfall))
		for (day in 1:length_period)
		{
			if(day<=di)
			{rain[day]<-mean(rainfall[1:day])}
			else
			{rain[day]<-mean(rainfall[((day-di)+1):day])}
		Average[day]<-rain[day]
		}
return(Average)
}
########################################################################################################################################################
# Dormant mosquiot reactivation fitness function
Reactivate_rain <-function(rainfall,moisture_retention,Rain_Threshold,h_step)
{
  di<-moisture_retention*(1/h_step)
  length_period<-length(rainfall)
  Fitness<-numeric(length(rainfall))
  rain<-numeric(length(rainfall))
  Mo=Rain_Threshold #saturation limit.
  for (day in 1:(length_period-1))
  {
    if((day<=di) && (mean(rainfall[1:day])<=rainfall[day+1]))
    {rain[day]<-mean(rainfall[1:day])
     Fitness[day]<-1/(1+1*exp(2*(Mo-rain[day])))}
    else if ((day>di) && (mean(rainfall[(day-di+1):day])<=rainfall[day+1]))
    {rain[day]<-mean(rainfall[((day-di)+1):day])
    Fitness[day]<-1/(1+1*exp(2*(Mo-rain[day])))}
    else Fitness[day]<-0
  }
  return(Fitness)
}

## Temperature fitness using the gumbel (extreme value) distribution
Fitness_temp<-function(Temp,Period,k1,k2)
{
	gumbel_fun<-exp(-0.001*(Temp-25)^2/20)*exp(-0.001*exp(0.01*(Temp-25)/20)) #fitness-generalised extreme value function   
	Temp_F<-gumbel_fun/(max(gumbel_fun)) # normalised fitness
	return(Temp_F)
}

#########################################################################################################################################################
# Mosquito dry season adaptation fitness function, using rainfall 7 day moving average
# Adaptation is based on normal distribution with means=optimal rainfall for adaptation
# Deviation, the window interval for adaptation
Fitness_rain_Adaptation<-function(Rainfall,Period,h_step,Std_rain,Rain_Threshold)
{
  mean_R=Rain_Threshold # mean optimal rainfall for adaptation
  dv=Std_rain # deviation
  no_days<-Period*(1/h_step)
  length_period<-length(Rainfall)
  Adapt<-numeric(length(Rainfall))
  rain<-numeric(length(Rainfall))
  for (day in 1:(length_period-1))
  {
        if((day<=no_days) && (mean(Rainfall[1:day])>=Rainfall[day+1]))
         {rain[day]<-mean(Rainfall[1:day])
          Adapt[day]<-exp(-(0.5*rain[day]-mean_R)^2/(dv^2))}
        else if ((day>no_days) && (mean(Rainfall[(day-no_days+1):day])>=Rainfall[day+1]))
         {rain[day]<-mean(Rainfall[((day-no_days)+1):day])
          Adapt[day]<-exp(-0.5*(rain[day]-mean_R)^2/(dv^2)) }
        else Adapt[day]<-0     
      
  }
  Adapt_F<-Adapt/(max(Adapt)) # normalised fitness
  return(Adapt_F)
}

######################################################################################################################################################
#Interpolate data 
Interpolate_data<-function(Data_var,h_step)
	{Data_obs<-Data_var
    time<-seq(0,length(Data_obs)*30.5, by=h_step)
		n=length(Data_obs)
		Tn<-seq(0,length(Data_obs)*30.5,length=n)
		Obs_simulated<-approx(Tn, Data_obs, time, method="linear", n=length(time), 1, length(Data_obs)*30.5, rule=1, f=0)
		Interpolated_data<-Obs_simulated$y
	return(Interpolated_data)	
	}

#############################################################################################################################################################	
#Vector population ecology model without aestivation (No-aestivation Model)
Vector_2qn_seasonal_Eff = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=7
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	Rain_Threshold=T_val
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  #print(MSQ_FT[t,1])
      Kmax<-rain_Kmax[t]+R_min
	  #Kmax=100000
      F_R<-Rain_F[t]
	  if((MSQ_FT[t,2]<CoFF1)&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   MSQ_FT[t,2]=0
	   }
	   
	   if(((MSQ_FT[t,1]<CoFF2))&&(t>1))
	   {
	   MSQ_FT[t,1]=0
	   #MSQ_FT[t,2]=0
	   }
      
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2] -p*MSQ_FT[t,1] -ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t,1] - um*MSQ_FT[t,2])  #Non aquatic vectors
        }
    return(MSQ_FT)
	})
}

Vector_2qn_seasonal_Euler_NULL = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=7
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  #R_min<-0.01*2*(um/p)*((ui*Kc)/((F*p/(2*um))-p-ui))-max(rain_Kmax)*0
	  #print(R_min)
	  #print(max(rain_Kmax))
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  #print(MSQ_FT[t,1])
      Kmax<-rain_Kmax[t]+R_min
	  #Kmax=100000
      F_R<-1
	  #MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2] -p*MSQ_FT[t,1] -ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      #MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t,1] - um*MSQ_FT[t,2])  #Non aquatic vectors

	  I = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2] -p*MSQ_FT[t,1] -ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      A = MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t,1] - um*MSQ_FT[t,2])  #Non aquatic vectors
     
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*A -p*I -ui*I*(1+Kc*(I/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*p*I - um*A)  #Non aquatic vectors
        }
    return(MSQ_FT)
	})
}

Vector_2qn_RK4_seasonal_NULL = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=7
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  #R_min<-0.01*2*(um/p)*((ui*Kc)/((F*p/(2*um))-p-ui))-max(rain_Kmax)*0
	  #print(R_min)
	  #print(max(rain_Kmax))
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  #print(MSQ_FT[t,1])
      Kmax<-rain_Kmax[t]+R_min
	  #Kmax=100000
      F_R<-1
	  #K1 and L1
	  K1 = h_step*(F_R*F*MSQ_FT[t,2] -p*MSQ_FT[t,1] -ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      L1 = h_step*(0.5*p*MSQ_FT[t,1] - um*MSQ_FT[t,2])  #Non aquatic vectors
      #K2 and L2
	  K2 = h_step*(F_R*F*(MSQ_FT[t,2]+0.5*L1) -p*(MSQ_FT[t,1]+0.5*K1) -ui*(MSQ_FT[t,1]+0.5*K1)*(1+Kc*((MSQ_FT[t,1]+0.5*K1)/Kmax)))  #Aquatic vectors
      L2 = h_step*(0.5*p*(MSQ_FT[t,1]+0.5*K1) - um*(MSQ_FT[t,2]+0.5*L1))  
      #K3 and L3
	  K3 = h_step*(F_R*F*(MSQ_FT[t,2]+0.5*L2) -p*(MSQ_FT[t,1]+0.5*K2) -ui*(MSQ_FT[t,1]+0.5*K2)*(1+Kc*((MSQ_FT[t,1]+0.5*K2)/Kmax)))  #Aquatic vectors
      L3 = h_step*(0.5*p*(MSQ_FT[t,1]+0.5*K2) - um*(MSQ_FT[t,2]+0.5*L2))  
      #K4 and L4	  
	  K4 = h_step*(F_R*F*(MSQ_FT[t,2]+L3) -p*(MSQ_FT[t,1]+K3) -ui*(MSQ_FT[t,1]+K3)*(1+Kc*((MSQ_FT[t,1]+K3)/Kmax)))  #Aquatic vectors
      L4 = h_step*(0.5*p*(MSQ_FT[t,1]+K3) - um*(MSQ_FT[t,2]+L3))  
      #Solution
      #MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+(K1+2*K2+2*K3+K4)/6  #Aquatic vectors
      #MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+(L1+2*L2+2*L3+L4)/6  #Non aquatic vectors
	  Aq = MSQ_FT[t,1]+(K1+2*K2+2*K3+K4)/6  #Aquatic vectors
      Ad = MSQ_FT[t,2]+(L1+2*L2+2*L3+L4)/6  #Non aquatic vectors
	  if((Aq<1.0e-5))
	  {
	  MSQ_FT[t + 1, 1] = 0  #Aquatic vectors
      MSQ_FT[t + 1, 2] = Ad  #Non aquatic vectors
	  }
	  else {
	  MSQ_FT[t + 1, 1] = Aq  #Aquatic vectors
      MSQ_FT[t + 1, 2] = Ad  #Non aquatic vectors
	  }
	  }
    return(MSQ_FT)
	})
}


Vector_2qn_seasonal_NULL = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=7
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  
   for (t in 1:n_steps) {
      Kmax<-rain_Kmax[t]+R_min
	  F_R<-1
	  
       if((MSQ_FT[t,2]<CoFF1)&&(t>1))
	   {
	   MSQ_FT[t,1]=0
	   MSQ_FT[t,2]=0
	   }
	   
	   if(((MSQ_FT[t,1]<CoFF2))&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   #MSQ_FT[t,2]=0
	   }
            
	  I = MSQ_FT[t,1]-(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))/(-p-ui*(1+2*Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      #print(I)
	  MSQ_FT[t + 1, 2] = (MSQ_FT[t,2]+h_step*(0.5*p*I))/(1+h_step*um)  #Non aquatic vectors
      MSQ_FT[t + 1, 1] = (MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t+1,2]))/(1+h_step*(p+ui*(1+Kc*(I/Kmax))))  #Aquatic vectors
	  	  
	  MSQ_FT[t + 1, 2] = (MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t+1,1]))/(1+h_step*um)  #Non aquatic vectors
      MSQ_FT[t + 1, 1] = (MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t+1,2]))/(1+h_step*(p+ui*(1+Kc*(MSQ_FT[t+1,1]/Kmax))))  #Aquatic vectors
     }
    return(MSQ_FT)
	})
}

Vector_2qn_seasonal_NULL2 = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=7
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  
   for (t in 1:n_steps) {
      Kmax<-rain_Kmax[t]+R_min
	  F_R<-1
	  
               
	  I = MSQ_FT[t,1]-(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))/(-p-ui*(1+2*Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      #print(I)
	  MSQ_FT[t + 1, 2] = (MSQ_FT[t,2]+h_step*(0.5*p*I))/(1+h_step*um)  #Non aquatic vectors
      MSQ_FT[t + 1, 1] = (MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t+1,2]))/(1+h_step*(p+ui*(1+Kc*(I/Kmax))))  #Aquatic vectors
	  	  
	  MSQ_FT[t + 1, 2] = (MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t+1,1]))/(1+h_step*um)  #Non aquatic vectors
      MSQ_FT[t + 1, 1] = (MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t+1,2]))/(1+h_step*(p+ui*(1+Kc*(MSQ_FT[t+1,1]/Kmax))))  #Aquatic vectors
	  
	  if((MSQ_FT[t+1,2]<CoFF1)&&(t>1))
	   {
	   #MSQ_FT[t+1,1]=0
	   MSQ_FT[t+1,2]=0
	   }
	   
	   if(((MSQ_FT[t+1,1]<CoFF2))&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t+1,1]=0
	   #MSQ_FT[t,2]=0
	   }
     }
    return(MSQ_FT)
	})
}



Vector_2qn_seasonal_NULL_p= function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=7
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  #R_min<-0.01*2*(um/p)*((ui*Kc)/((F*p/(2*um))-p-ui))-max(rain_Kmax)*0
	  #print(R_min)
	  #print(max(rain_Kmax))
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  #print(MSQ_FT[t,1])
      Kmax<-rain_Kmax[t]+0.05
	  #Kmax=100000
      F_R<-1
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2] -p*MSQ_FT[t,1] -ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(a*0.5*p*MSQ_FT[t,1] - um*MSQ_FT[t,2])  #Non aquatic vectors
        }
    return(MSQ_FT)
	})
}

##############################################################################################################################################
Vector_3qn_seasonal_Eff = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	Period_T=7
    #Rainfall_sd=20
    #Optimal_rain=80
		
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  Aestivation_F<-Fitness_rain_Adaptation(rainfall,Period=Period_T,h_step,Std_rain=T_val/2,Rain_Threshold=T_val)
      Reactivate_F<-Reactivate_rain(rainfall,soil_dry_time,Rain_Threshold=T_val,h_step)
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  F_ACT<-Reactivate_F[t]
      F_Adp<-Aestivation_F[t]
	  Kmax<-rain_Kmax[t]+R_min
	  F_R<-Rain_F[t]
	  
	  if((MSQ_FT[t,2]<CoFF1)&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   MSQ_FT[t,2]=0
	   }
	   
	   if(((MSQ_FT[t,1]<CoFF2))&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   #MSQ_FT[t,2]=0
	   }
      
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*((MSQ_FT[t,1])/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t,1]- um*MSQ_FT[t,2]+w*F_ACT*MSQ_FT[t,3]-d*F_Adp*MSQ_FT[t,2] )  #Non aquatic vectors
	  MSQ_FT[t + 1, 3] = MSQ_FT[t,3]+h_step*(-umd*MSQ_FT[t,3]-w*F_ACT*MSQ_FT[t,3] +d*F_Adp*MSQ_FT[t,2])  #Dormant adult vectors
	  MSQ_FT[t + 1, 4] = (1+Kc*((MSQ_FT[t,1])/Kmax))
        }
    return(MSQ_FT)
	})
}

Vector_3qn_seasonal_NULL= function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	Period_T=7
    #Rainfall_sd=20
    #Optimal_rain=80
		
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  Aestivation_F<-Fitness_rain_Adaptation(rainfall,Period=Period_T,h_step,Std_rain=T_val/2,Rain_Threshold=T_val)
      Reactivate_F<-Reactivate_rain(rainfall,soil_dry_time,Rain_Threshold=T_val,h_step)
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  F_ACT<-Reactivate_F[t]
      F_Adp<-Aestivation_F[t]
	  Kmax<-rain_Kmax[t]+R_min
	  F_R<-1
	  
	  if((MSQ_FT[t,2]<CoFF1)&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   MSQ_FT[t,2]=0
	   }
	   
	   if(((MSQ_FT[t,1]<CoFF2))&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   #MSQ_FT[t,2]=0
	   }
      
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*((MSQ_FT[t,1])/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t,1]- um*MSQ_FT[t,2]+w*F_ACT*MSQ_FT[t,3]-d*F_Adp*MSQ_FT[t,2] )  #Non aquatic vectors
	  MSQ_FT[t + 1, 3] = MSQ_FT[t,3]+h_step*(-umd*MSQ_FT[t,3]-w*F_ACT*MSQ_FT[t,3] +d*F_Adp*MSQ_FT[t,2])  #Dormant adult vectors
	  MSQ_FT[t + 1, 4] = (1+Kc*((MSQ_FT[t,1])/Kmax))
        }
    return(MSQ_FT)
	})
}

Vector_3qn_seasonal_Robust = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	Period_T=7
    Rainfall_sd=20
    Optimal_rain=80
	Emin<-0.05
	robust=1.0
	 
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  Aestivation_F<-Fitness_rain_Adaptation(rainfall,Period=Period_T,h_step,Std_rain=T_val/4,Rain_Threshold=T_val)
      Reactivate_F<-Reactivate_rain(rainfall,soil_dry_time,Rain_Threshold=T_val,h_step)
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  F_ACT<-Reactivate_F[t]
      F_Adp<-Aestivation_F[t]
	  Kmax<-rain_Kmax[t]+Emin
	  if((t>75/h_step)&&(t<275/h_step))
	  {robust=0.0}
	  else {robust=1.0}
	  
	  F_R<-Rain_F[t]
      Aquatic = MSQ_FT[t,1]+h_step*(robust*F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*((MSQ_FT[t,1])/Kmax)))  #Aquatic vectors
      Adult_A = MSQ_FT[t,2]+h_step*(0.5*robust*p*MSQ_FT[t,1]- um*MSQ_FT[t,2]+w*F_ACT*MSQ_FT[t,3]-d*F_Adp*MSQ_FT[t,2] )  #Non aquatic vectors
	  Adult_D = MSQ_FT[t,3]+h_step*(-umd*MSQ_FT[t,3]-w*F_ACT*MSQ_FT[t,3] +d*F_Adp*MSQ_FT[t,2])  #Dormant adult vectors
	  const = (1+Kc*((MSQ_FT[t,1])/Kmax))
	  
	  if(Adult_A<0.00001)
	  {
	  MSQ_FT[t + 1, 1] =0 
      MSQ_FT[t + 1, 2] =0 
	  MSQ_FT[t + 1, 3] =Adult_D 
	  MSQ_FT[t + 1, 4] =const
	  }
	  else{
	  MSQ_FT[t + 1, 1] =Aquatic 
      MSQ_FT[t + 1, 2] =Adult_A 
	  MSQ_FT[t + 1, 3] =Adult_D 
	  MSQ_FT[t + 1, 4] = const}
        }
    return(MSQ_FT)
	})
}


Vector_3qn_seasonal_Eff_p = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	Period_T=7
    Rainfall_sd=20
    Optimal_rain=80
	Emin<-0.05
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  Aestivation_F<-Fitness_rain_Adaptation(rainfall,Period=Period_T,h_step,Std_rain=T_val/4,Rain_Threshold=T_val)
      Reactivate_F<-Reactivate_rain(rainfall,soil_dry_time,Rain_Threshold=T_val,h_step)
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  F_ACT<-Reactivate_F[t]
      F_Adp<-Aestivation_F[t]
	  Kmax<-rain_Kmax[t]+R_min
	  F_R<-Rain_F[t]
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*((MSQ_FT[t,1])/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(a*0.5*p*MSQ_FT[t,1]- um*MSQ_FT[t,2]+w*F_ACT*MSQ_FT[t,3]-d*F_Adp*MSQ_FT[t,2] )  #Non aquatic vectors
	  MSQ_FT[t + 1, 3] = MSQ_FT[t,3]+h_step*(-umd*MSQ_FT[t,3]-w*F_ACT*MSQ_FT[t,3] +d*F_Adp*MSQ_FT[t,2])  #Dormant adult vectors
	  MSQ_FT[t + 1, 4] = (1+Kc*((MSQ_FT[t,1])/Kmax))
        }
    return(MSQ_FT)
	})
}

Vector_3qn_seasonal_Eff_Sense = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	#Period_T=7
    #Rainfall_sd=20
    #Optimal_rain=80
	#Emin<-0.05
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  Aestivation_F<-Fitness_rain_Adaptation(rainfall,Period=Period_T,h_step,Std_rain=T_val/4,Rain_Threshold=T_val)
      Reactivate_F<-Reactivate_rain(rainfall,soil_dry_time,Rain_Threshold=T_val,h_step)
   for (t in 1:n_steps) {
      #x_step=start+(t-1)*h_step
	  F_ACT<-Reactivate_F[t]
      F_Adp<-Aestivation_F[t]
	  Kmax<-rain_Kmax[t]+R_min
	  F_R<-Rain_F[t]
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*((MSQ_FT[t,1])/Kmax)))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*p*MSQ_FT[t,1]-um*MSQ_FT[t,2]+w*F_ACT*MSQ_FT[t,3]-d*F_Adp*MSQ_FT[t,2] )  #Non aquatic vectors
	  MSQ_FT[t + 1, 3] = MSQ_FT[t,3]+h_step*(-umd*MSQ_FT[t,3]-w*F_ACT*MSQ_FT[t,3] +d*F_Adp*MSQ_FT[t,2])  #Dormant adult vectors
	  MSQ_FT[t + 1, 4] = (1+Kc*((MSQ_FT[t,1])/Kmax))
        }
    return(MSQ_FT)
	})
}

####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
########################### Interventions model #####################################################################################################################################################################

Vector_4qn_seasonal_Eff_Control = function(p_fix,h_step,V_0,p_vary)
{
    start=0
    soil_dry_time=7
	t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	Period_T=7
    #Rainfall_sd=20
    #Optimal_rain=80
	#Emin<-0.05

	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
   {  
	 Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
	  Aestivation_F<-Fitness_rain_Adaptation(rainfall,Period=Period_T,h_step,Std_rain=T_val/4,Rain_Threshold=T_val)
      Reactivate_F<-Reactivate_rain(rainfall,soil_dry_time,Rain_Threshold=T_val,h_step)
	
	month=31/h_step
	Start_W=9*month
	Start_D=2*month
	Interval=6*month
	Yr=month*12
	sw0=sw1=sw2=dw0=dw1=dw2=0.0 #sw=switch for wet season, dw = switch for dry season, 0=larvicides, 1=IRS, 2=LLINs
                              # sw=0.0 or dw=0.0, means intervention is off, if (=1.0) intervention is on 
	count=0
	#print(count)
	irs=intervene
	lcid=intervene
	llin=intervene
	#print(intervene)

	for (t in 1:n_steps){
		#x_step=start+(t-1)*h_step
		######################################################################################################
		##### Interventions during the wet season
		######################################################################################################
		if((season==0)&&(opt==1)) #Season = 0 means intervene in the wet season, opt==1 means IRS spraying 
		{
		#print('here')
		if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw1=1 # IRS intervention 
			#print(sw1)
			sw0=sw2=0}# other interventions are switched off
	   else {
			sw0=sw1=sw2=0 # All interventions are switched off
			#print('this')
			}
		}
		
		else if ((season==0)&&(opt==2)) #Season = 0 means intervene in the wet season
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw2=1 # LLINs treated bed nets interventions 
			sw0=sw1=0}# other interventions are switched off
			else if(count>=3){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1
			sw0=sw1=0}# treatment off 
		}
		
		else if ((season==0)&&(opt==3)) #Season = 0 means wet season and opt=3, means larvicidal interventions
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=1
			sw1=sw2=0}# means treatment on
			else {
			sw0=sw1=sw2=0}# treatment off
		}
		
		else if ((season==0)&&(opt==4)) #Season = 0 means wet season and opt=4, means both IRS and LLIN
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=0
			sw1=sw2=1}# means treatment on
			else if(count>=3){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1 
			sw0=sw1=sw2=0}# treatment off
		}
		
		else if ((season==0)&&(opt==5)) #Season = 0 means wet season and opt=5, means IRS (or LLINs) and larvicides
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=sw1=1
			sw2=0}# means treatment on
			else if(count>=3){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1
			sw0=sw1=sw2=0}# treatment off
		}
	
		else if ((season==0)&&(opt==6)) #Season = 0 means wet season and opt=5, means IRS,LLIN and larvicides
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=sw1=sw2=1}# means treatment on
			else if(count>=3){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1
			sw0=sw1=sw2=0}# treatment off
		}
	
	######################################################################################################
	##### Interventions during the dry season
    ######################################################################################################
		else if((season==1)&&(opt==1)) #Season = 1 means intervene in the dry season, opt==1 means IRS spraying 
		{
		#print('here')
			if ((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3)))
			{dw1=1 # IRS intervention 
			#print(sw1)
			sw0=dw0=sw2=dw2=0}# other interventions are switched off
			else {
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 # All interventions are switched off
			#print('this')
			}
		}
		
		else if ((season==1)&&(opt==2)) #Season 1 means dry season and opt=2 means LLINs
			{
			if ((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3)))
			{dw2=1 # bet nets treated interventions 
			sw0=sw1=0}# other interventions are switched off
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=0.0 }# treatment off
			else {
			dw2=1
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off 
		}
	
		else if ((season==1)&&(opt==3)) #Season 1 means dry season and opt=3, means larvicidal interventions
		{
			if ((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3)))
			{dw0=1
			sw0=sw1=sw2=dw1=dw2=0.0 }# means treatment on
			else {
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
		}
		
		else if ((season==1)&&(opt==4)) #Season = 1 means dry season and opt=4, means both IRS and LLIN
		{
			if ((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3)))
			{dw0=0
			dw1=dw2=1}# means treatment on
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
			else {
			dw2=1 
			sw0=sw1=sw2=dw0=dw1=0.0 }# treatment off
		}
	
		else if ((season==1)&&(opt==5)) #Season = 1 means dry season and opt=5, means IRS,LLIN and larvicidals
			{
			if ((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3)))
			{dw0=dw1=dw2=1}# means treatment on
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
			else {
			dw2=1
			sw0=sw1=sw2=dw0=dw1=0.0}# treatment off
		}	 
	
	######################################################################################################
	##### Interventions for both active and dormant vectors
    ######################################################################################################
		else if((season==2)&&(opt==1)) #Season = 2 means intervene both active and dormant vectors dry/wet season, opt==1 means IRS spraying 
		{
		#print('here')
			if (((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval)&&(count<3)))||((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3))))
			{sw1=dw1=1 # IRS intervention 
			#print(sw1)
			sw0=dw0=sw2=dw2=0}# other interventions are switched off
			else {
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 # All interventions are switched off
			#print('this')
				}
		}
    
		else if ((season==2)&&(opt==2)) #Season = 2 means intervene both, opt=LLINs
			{
			if (((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval)&&(count<3)))||((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3))))			{sw2=dw2=1 # bet nets treated interventions 
			sw0=sw1=0}# other interventions are switched off
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=0.0 }# treatment off
			else {
			sw2=dw2=1
			sw0=sw1=dw0=dw1=0.0 }# treatment off 
		}
	
		else if ((season==2)&&(opt==3)) #Season = 2 means intervene both and opt=3, means larvicidal interventions
		{
			if (((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval)&&(count<3)))||((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3))))
			{dw0=sw0=1
			sw1=sw2=dw1=dw2=0.0 }# means treatment on
			else {
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
		}
	
		else if ((season==2)&&(opt==4)) #Season = 2 means intervene both and opt=4, means both IRS and LLIN
		{
			if (((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval)&&(count<3)))||((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3))))
			{dw0=0
			sw1=sw2=dw1=dw2=1}# means treatment on
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
			else {
			sw2=dw2=1 
			sw0=sw1=sw2=dw0=dw1=0.0 }# treatment off
		}
		
		else if ((season==2)&&(opt==5)) #Season = 2 means intervene both and opt=5, means both IRS and Larvicides
		{
			if (((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval)&&(count<3)))||((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3))))
			{sw2=dw2=0
			sw1=sw0=dw1=dw0=1}# means treatment on
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
			else {
			sw2=dw2=1 
			sw0=sw1=sw2=dw0=dw1=0.0 }# treatment off
		}
		
		else if ((season==2)&&(opt==6)) #Season = 2 means intervene both and opt=6 means IRS,LLIN and larvicides
		{
			if (((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval)&&(count<3)))||((t>(Start_D+Yr*count))&&(t<(Start_D+Yr*count+Interval)&&(count<3))))
			{sw0=sw1=sw2=dw0=dw1=dw2=1}# means treatment on
			else if(count>=3){
			dw2=0
			sw0=sw1=sw2=dw0=dw1=dw2=0.0 }# treatment off
			else {
			dw2=1
			sw0=sw1=sw2=dw0=dw1=0.0}# treatment off
		} 
	##############################################################################################
	#### End of interventions
	##############################################################################################
	
	####################### end of treatment options#######################################
		F_ACT<-Reactivate_F[t]
		F_Adp<-Aestivation_F[t]
		Kmax<-rain_Kmax[t]+R_min
		F_R<-Rain_F[t]
		if(((t%%(Yr))==0)&&(count<=3))
			{count=count+1}
		else {count=count}
		
		if((MSQ_FT[t,2]<CoFF1)&&(t>1)&&(Kmax<T_val))
	   {
	   #MSQ_FT[t,1]=0
	   MSQ_FT[t,2]=0
	   }
	   
	   if(((MSQ_FT[t,1]<CoFF2))&&(t>1)&&(Kmax<T_val))
	   {
	   MSQ_FT[t,1]=0
	   #MSQ_FT[t,2]=0
	   }
      
		MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2]-p*MSQ_FT[t,1]-ui*MSQ_FT[t,1]*(1+Kc*((MSQ_FT[t,1])/(Kmax*(1-lcid*sw0)))))  #Aquatic vectors
		MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*(1-lcid*sw0)*p*MSQ_FT[t,1]-um*MSQ_FT[t,2]+w*F_ACT*MSQ_FT[t,3]-d*F_Adp*MSQ_FT[t,2]-(irs*dw1+llin*sw2+irs*sw1)*MSQ_FT[t,2] )  #Non aquatic vectors
		MSQ_FT[t + 1, 3] = MSQ_FT[t,3]+h_step*(-umd*MSQ_FT[t,3]-w*F_ACT*MSQ_FT[t,3] +d*F_Adp*MSQ_FT[t,2]-(llin*dw2+irs*dw1)*MSQ_FT[t,3])  #Dormant adult vectors
		MSQ_FT[t + 1, 4] = (1+Kc*((MSQ_FT[t,1])/Kmax))
	    }
        return(MSQ_FT)}
	)
}

# intervention using the no-aestivation model
Vector_2qn_seasonal_Control = function(p_fix,h_step,V_0,p_vary)
{
    
    start=0
    soil_dry_time=7
	Period_T=4
    t_step<-length(T_R$Rain)*30.5 
	n_steps<-(t_step-start)/h_step
	#time<-seq(start,t_step,by=h_step)
    rainfall<-Interpolate_data(T_R$Rain,h_step)
	
	Nstg = length(V_0)  #Number of stages/classes
    MSQ_FT = matrix(0, (n_steps+1), Nstg)  #Initialize output matrix with time steps as rows, age classes as columns
    MSQ_FT[1, ] = V_0  #put initial values into first row of output matrix
    with(as.list(c(p_fix,p_vary)),
	{  
	  Rain_F<-Fitness_rain(rainfall,soil_dry_time,h_step,Rain_Threshold=T_val)
	  rain_Kmax<-Average_rain(rainfall,soil_dry_time,h_step)
   	  month=31/h_step
      Start_W=9*month
      Start_D=4*month
      Interval=6*month
      Yr=month*12
      sw0=sw1=sw2=dw0=dw1=dw2=0.0 #sw=switch for wet season, d = switch for dry season, 0=larvicides, 1=IRS, 2=LLINs
                              # sw=0.0 or dw=0.0, means intervention is off, if (=1.0) intervention is on 
     count=0
     irs=intervene
     lcid=intervene
     llin=intervene
	  for (t in 1:n_steps) {
	    if((season==0)&&(opt==1)) #Season = 0 means intervene in the wet season, opt==1 means IRS spraying 
		{
			#print('here')
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw1=1 # IRS intervention 
			#print(sw1)
			sw0=sw2=0}# other interventions are switched off
	   else {
			sw0=sw1=sw2=0 # All interventions are switched off
			#print('this')
			}
		}
		else if ((season==0)&&(opt==2)) #Season = 0 means intervene in the wet season
			{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw2=1 # LLINs treated bed nets interventions 
			sw0=sw1=0}# other interventions are switched off
			else if((count>=3)|(count<1)){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1
			sw0=sw1=0}# treatment off 
		}
		else if ((season==0)&&(opt==3)) #Season = 0 means wet season and opt=3, means larvicidal interventions
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=1
			sw1=sw2=0}# means treatment on
			else {
			sw0=sw1=sw2=0}# treatment off
		}
		else if ((season==0)&&(opt==4)) #Season = 0 means wet season and opt=4, means both IRS and LLIN
		{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=0
			sw1=sw2=1}# means treatment on
			else if((count>=3)|(count<1)){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1 
			sw0=sw1=0}# treatment off
		}
		else if ((season==0)&&(opt==5)) #Season = 0 means wet season and opt=5, means IRS (or LLINs) and larvicides
			{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=sw1=1
			sw2=0}# means treatment on
			else if((count>=3)|(count<1)){
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1
			sw0=sw1=0}# treatment off
		}
		else if ((season==0)&&(opt==6)) #Season = 0 means wet season and opt=5, means IRS,LLIN and larvicides
			{
			if ((t>(Start_W+Yr*count))&&(t<(Start_W+Yr*count+Interval))&&(count<3))
			{sw0=sw1=sw2=1}# means treatment on
			else if((count>=3)|(count<1))
			{
			sw2=0
			sw0=sw1=0}# treatment off
			else {
			sw2=1
			sw0=sw1=0}# treatment off
		}
	
	   if(((t%%(Yr))==0)&&(count<=3))
	    {count=count+1}
	   else {count=count}
	
	################################ Treatment options#######
      Kmax<-rain_Kmax[t]+R_min
	  F_R<-Rain_F[t]
	  if((MSQ_FT[t,2]<CoFF1)&&(t>1))
	   {
	   #MSQ_FT[t,1]=0
	   MSQ_FT[t,2]=0
	   }
	   
	   if(((MSQ_FT[t,1]<CoFF2))&&(t>1))
	   {
	   MSQ_FT[t,1]=0
	   #MSQ_FT[t,2]=0
	   }
      
      MSQ_FT[t + 1, 1] = MSQ_FT[t,1]+h_step*(F_R*F*MSQ_FT[t,2] -p*MSQ_FT[t,1] -ui*MSQ_FT[t,1]*(1+Kc*(MSQ_FT[t,1]/(Kmax*(1-lcid*sw0)))))  #Aquatic vectors
      MSQ_FT[t + 1, 2] = MSQ_FT[t,2]+h_step*(0.5*(1-lcid*sw0)*p*MSQ_FT[t,1] - um*MSQ_FT[t,2] -(llin*sw2+irs*sw1)*MSQ_FT[t,2])  #Non aquatic vectors
         
      }
    return(MSQ_FT)
	})
}

