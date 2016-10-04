#This file simulates an epidemiological model with 
#periodic recovery rates. 
#Research with Geiser

library(deSolve)
rm(list=ls())



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
         dsusceptibles<- mu - mu*susceptibles - beta*susceptibles*infected
         dinfected<- beta*susceptibles*infected - periodic*infected - mu*infected
         drecovered<- periodic*infected - mu*recovered
         return(list(c(dsusceptibles,dinfected,drecovered)))
       })#this parenthesis closes the with
}
#do a for loop to explore parameter space
#define sequences for three parameters, iota, mu and beta
seqiota<-seq(from=0.0001,to=2,length.out = 100)
seqmu<-seq(from=0.00001, to=0.2, length.out = 100)
seqbeta<-seq(from=1, to=100, length.out = 100)
parspace<-matrix(NA,nrow=1000000,ncol=4)
colnames(parspace)<-c("iota","mu","beta","ampl")
row<-0
for (i in seqiota)
{
  for (m in seqmu)
  {
    for (b in seqbeta)
    {
      #define parameter values
      pars<-list(iota=i, mu=m, beta=b)
      out<-ode(y=state,times=times,func=epi,parms = pars)
      row<-row+1
      #plot(out,ylim=c(0.0044,0.0045),xlim=c(60,100))
      
      #get the amplitude in the infected
      inf<-out[,2]
      #shorten it to the last 50
      long<-length(inf)
      inf<-inf[(long-50):long]
      ampl<-max(inf)-min(inf)
      parspace[row,]<-c(unlist(pars),ampl=ampl)
    }
  }
}
