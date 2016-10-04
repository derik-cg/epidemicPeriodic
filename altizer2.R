#This file replicates the simulations in box 2 of the Altizer paper
#Ecology letters (2006) 9: 467-484

library(deSolve)
rm(list=ls())

#define initial conditions
state<-c(S=(1E6-15), E=0, I=15)

#define the vector of output times
times<-seq(from=0,to=50,by=0.1)

#define parameter values
pars<-list(mu=0, sigma=45.625, gama=73, lambda0=0.02, 
           lambda1=1, beta0=1250, beta1=0.21, N=1E6)

#define the differential equations
epi<-function(t,state,pars)
{
  with(as.list(c(t,state,pars)),
       {
         #compute the lambda and beta for each time step
         lambda<-lambda0+(1+lambda1*cos(2*pi*t))
         #the next line is for a square lambda wave. 
         #lambda<-lambda0*(1+lambda1*)
         beta<-beta0*(1+beta1*if((t%%1 > 1/4)&(t%%1 < 3/4))-1 else 1)
         #here the ODEs
         dS<- lambda*(N-S)+mu*I-beta*S*(I/N)
         dE<- beta*S*(I/N)-lambda*E-sigma*E
         dI<- sigma*E-(gama+mu+lambda)*I
         return(list(c(dS,dE,dI)))
       })#this parenthesis closes the with
}

out<-ode(y=state,times=times,func=epi,parms = pars)
plot(out,ylim=c(0,1200))
out2<-out
#there is a transient period that ends about time 40
plot(out2[,1],out2[,4],xlim=c(0,5),ylim=c(400,35000),type="l")
