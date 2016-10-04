#simulation of equatios in Atlizer et al 2006 
#ecology letters 2006 9 467-484}
rm(list=ls())
#define parameters
pars<-list(mu=0, beta0=1250, beta1=0.1, sigma=45.625, gama=73, 
           lambda0=0.02, lambda1=1)

#definetimes
times<-seq(from=0,to=50,by=0.5)

#define state variables

state<-c(S=(5E6)-1, E=0, I=1, R=0, N=5E6)

#define equations
epi<-function(t,state,pars)
{
  with(as.list(c(t,state,pars)),
       {
         beta<-function(t)
         {
           beta0*(1+beta1*cos(2*pi*t))
         }
         
         lambda<-function(t)
         {
           lambda0*(1+lambda1*cos(2*pi*t))
         }
         
         dS<-lambda(t)*(N-S)+mu*I-beta(t)*S*I/N
         dE<-beta(t)*S*I/N-lambda(t)*E-sigma*E
         dI<-sigma*E-(gama+mu+lambda(t))*I
         R<-N-(S+E+I)
       })#close the with
  return(list(c(dS,dE,dI)))
}
out<-ode(y=state,times=times,func=epi,parms = pars)