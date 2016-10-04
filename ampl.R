ampl<-function(st)
{
  #this function takes as the argument a time series, 
  #and computes the amplitude, assuming the serie has 
  #is stationary and has a strong deterministic 
  #component. Amplitude is measured as the difference of
  #the average maximum and minimum 
  
  #walk the time series and find maxima
  n<-length(st) #find the size
  maxes<-0 #array where maxes will be stored
  mins<-0 #array where mins will be stored
  
  for (i in 2:(n-1)) #for each step in the st
  { #find maximums and store in maxes
    if((st[i-1] < st[i]) & (st[i] > st[i+1]))
    {
      maxes<-c(maxes,st[i])
    }
    #find minima and store them in mins
    if((st[i-1] > st[i]) & (st[i] < st[i+1]))
    {
      mins<-c(mins,st[i])
    }
  }
  #compute the averages
  avgmax<-mean(maxes)
  avgmin<-mean(mins)
  
  return(c(avgmax,avgmin))
}