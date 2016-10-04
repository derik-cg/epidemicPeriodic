#This file simulates an epidemiological model with 
#periodic recovery rates. 
#Research with Geiser

library(deSolve)
#rm(list=ls())



#define initial conditions
state<-c(susceptibles=0.4, infected=0.1, recovered=0)

#define the vector of output times
times<-seq(from=0,to=1000,by=0.5)

#define parameter values
pars<-list(iota=0.5, mu=0.05, beta=5)#mu=0.0041

#define the differential equations
epi<-function(t,state,pars)
{
  with(as.list(c(t,state,pars)),
       {
         periodic<-iota*(1+0.1*cos((0.032)*pi*t+0.15))
         #here the ODEs
         dsusceptibles<- mu - mu*susceptibles - beta*susceptibles*infected
         dinfected<- beta*susceptibles*infected - periodic*infected - mu*infected
         drecovered<- periodic*infected - mu*recovered
         return(list(c(dsusceptibles,dinfected,drecovered)))
       })#this parenthesis closes the with
}
out<-ode(y=state,times=times,func=epi,parms = pars)
plot(out)
#plot(out,ylim=c(0.034,0.05),xlim=c(100,1000))
#plot(out,ylim=c(0.034,0.045),xlim=c(60,100))
#measure the amplitude of oscillation in infected
sz<-dim(out)
ampl<-c(0.1,max(out[100:1000,2])-min(out[100:1000,2]))
