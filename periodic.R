#This file simulates an epidemiological model with 
#periodic infection and recovery rates. 
#Research with Geiser

library(deSolve)
rm(list=ls())

#define parameter values
pars<-list(iota=0.01, mu=0.001, beta=20)

#define initial conditions
state<-c(susceptibles=0.4, infected=0.1, recovered=0)

#define the vector of output times
times<-seq(from=0,to=500,by=0.5)

#define the differential equations
epi<-function(t,state,pars)
{

  with(as.list(c(t,state,pars)),
  {
  periodic<-iota*(1+0.16*cos(2*pi*t+0.15))
  #here the ODEs
    dx<- mu - mu*susceptibles - beta*susceptibles*infected
    dy<- beta*susceptibles*infected - periodic*infected - mu*infected
    dz<- periodic*infected - mu*recovered
    return(list(c(dx,dy,dz)))
  })#this parenthesis closes the with
}

out<-ode(y=state,times=times,func=epi,parms = pars)
plot(out)