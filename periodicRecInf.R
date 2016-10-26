#This file simulates an epidemiological model with 
#periodic recovery and infection rates.
#the periodicity of both is the important part in this
#file. How do they interact?
#Research with Geiser

setwd("~/Investigacion/Epidemiology")
rm(list=ls())

library(deSolve)
source('~/Investigacion/Epidemiology/f.peak.R')
source('~/Investigacion/Epidemiology/fft.ts.R')
source('~/Investigacion/Epidemiology/fft.osc.R')

#define initial conditions
state<-c(susceptibles=0.4, infected=0.4, recovered=0.4)

#define the vector of output times
endt<-365*20
lag<-7
times<-seq(from=0,to=endt,by=lag)
ntimes<-length(times)
nserie<-ntimes-round(ntimes/2)

#define parameter values
pars<-list(mu=0.0041, beta=(1/30), gamma=(1/60), amp1=0.1, amp2=0.001, p1=(pi*50), p2=365)

#define the differential equations
epi<-function(t,state,pars)
{
  with(as.list(c(t,state,pars)),
       {
         recovery<-gamma+(amp1/2)*cos((2*pi)/p1*t)
         infection<-beta+(amp2/2)*cos((2*pi)/p2*t)
         #here the ODEs
         dsusceptibles<- mu - mu*susceptibles - infection*susceptibles*infected
         dinfected<- infection*susceptibles*infected - recovery*infected - mu*infected
         drecovered<- recovery*infected - mu*recovered
         return(list(c(dsusceptibles,dinfected,drecovered)))
       })#this parenthesis closes the with
}
#the event function
eventfunc<-function(t,state,pars)
  
{
  state[1]<-0
  state[2]<-0
  state[3]<-0
  return(state)
}
#the root function
rootfunc<-function(t,state,parms)
{
  state[1]
  state[2]
  state[3]
  return(state)
}
#wrap the ode function in a tryCatch block
robustode<-function(state,times,epi,pars,rootfunc,eventfunc)
{
  tryCatch(ode(y=state,times=times,func=epi,parms = pars,
               rootfunc=rootfunc,events=list(func=eventfunc,root=T,
                                             terminalroot=T)),
           error=function(e){;NA})
}
#here a hypercube for the most important parameters
#amp1 amp2 p1 p2
lim.amp1<-seq(from=0.01,to=0.5,length.out=10)
lim.amp2<-seq(from=0.01,to=0.5,length.out=10)
lim.p1<-c(0.01,7,15,30,60,90)
lim.p2<-c(0.01,7,15,30,60,90)
#store the output in a table
#the size of the table is
n<-10*10*6*6
osc.out<-array(NA,dim=c(n,25))
colnames(osc.out)<-c("amp1","amp2","p1","p2", #1-4
                     "min.sus","max.sus",     #5-6
                     "min.inf","max.inf",     #7-8
                     "min.rec","max.rec",     #9-10
                     "per.acf.sus","per.acf.inf","per.acf.rec", #11-13
                     "per.fft.1.sus","peak.fft.1.sus", #14-15
                     "per.fft.2.sus","peak.fft.2.sus", #16-17
                     "per.fft.1.inf","peak.fft.1.inf", #18-19
                     "per.fft.2.inf","peak.fft.2.inf", #20-21
                     "per.fft.1.rec","peak.fft.1.rec", #22-23
                     "per.fft.2.rec","peak.fft.2.rec") #24-25
rw<-0
for (a in lim.amp1)
{
  for (b in lim.amp2)
  {
    for (c in lim.p1)
    {
      for (d in lim.p2)
      { 
        rw<-rw+1
        #if(b<0.065) next
        pars$amp1<-a
        pars$amp2<-b
        pars$p1<-c
        pars$p2<-d
        ####################################
        #### check with na how does this behaves
        ####################################
        #do the simulation
        out<-robustode(state,times,epi,pars,
                 rootfunc,eventfunc)
        if (is.null(out)) out<-NA
        else if (length(out[,1])<600)
        {
          #the simulation is not large enough to generate
          #oscillations. return NA
          out<-NA
        }
        else if(length(out[,1])<ntimes)
        {# There was a convergence error, so the simulation
          #does not end at ntimes, but before.
          #find time step 600 days before end
          tpend<-length(out[,1])#end of time period
          #using a square to find the given time as a 
          #minimum
          timemin<-((out[,1]-((out[length(out[,1]),1])-600))^2)
          tpstart<-which(timemin==min(timemin))
        }
        else
        {
          tpend<-ntimes #the end time point is the end of simulation
          tpstart<-ntimes/2
        }
        print(rw/n)
        if (!is.na(out[1])) plot(out)
        if (!is.na(out[1]))
        {
          osc.out[rw,1:4]<-c(a,b,c,d) #write pars
          ampl<-c(max(out[(tpstart):tpend,2]),min(out[(tpstart):tpend,2]),
                  max(out[(tpstart):tpend,3]),min(out[(tpstart):tpend,3]),
                  max(out[(tpstart):tpend,4]),min(out[(tpstart):tpend,4]))
          osc.out[rw,5:10]<-c(ampl)
          #find the period of the oscillations using a correlogram
          pers.acf<-c(f.peak(acf(out[(tpstart):tpend,2],lag.max=100,plot=F)),
                      f.peak(acf(out[(tpstart):tpend,3],lag.max=100,plot=F)),
                      f.peak(acf(out[(tpstart):tpend,4],lag.max=100,plot=F)))
          osc.out[rw,11:13]<-pers.acf
          #find the harmonics and magnitudes of the fft
          hm.fft<-c(fft.osc(out[(tpstart):tpend,2],nserie,lag),
                    fft.osc(out[(tpstart):tpend,3],nserie,lag),
                    fft.osc(out[(tpstart):tpend,2],nserie,lag))
          osc.out[rw,14:25]<-hm.fft
        }
        else
        {
          osc.out[rw,]<-NA
        }
      }
    }
    save.image("~/Investigacion/Epidemiology/oscilaciones.RData")
  }
}
save.image("~/Investigacion/Epidemiology/oscilaciones.RData")
osc.out$raw.ampl<-osc.out$max.sus-osc.out$min.sus
plot(osc.out$amp1,osc.out$amp2)
# dev.off()
# #plot(out,ylim=c(0.034,0.05),xlim=c(100,1000))
# #plot(out,ylim=c(0.034,0.045),xlim=c(60,100))
# #measure the amplitude of oscillation in infected
# #sz<-dim(out)
# #ampl<-c(0.1,max(out[100:200,2])-min(out[100:200,2]))
# #span of the stationary series
# strts<-floor(endt/lag/3)
# endts<-floor(endt/lag)
# source("fftauto.R")
# #span of the stationary series
#strts<-floor(endt/lag/3)
#endts<-floor(endt/lag)
#plot(out[strts:endts,1],out[strts:endts,4],type="l")
#spec<-fft.ts(out[strts:endts,3],endt,lag)
# 

