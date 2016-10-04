fft.ts<-function(serie,n,delta)
{
  #the first part is to obtain the fft on the crude data
  fft.out<-fft(serie)
  #now compute the nyquist number
  nyquist<-floor(length(fft.out)/2)
  #now compute the magnitude for each harmonic
  magn<-sqrt(Re(fft.out[1:nyquist])^2+Im(fft.out[1:nyquist])^2)*2/nyquist
  #now compute the harmonics
  harm<-(n*delta)/(1:nyquist)
  #eliminate the first element in harmonics and magnitude
  harm<-harm[2:nyquist]
  magn<-magn[2:nyquist]
  #plot harmonics versus magnitude
  plot(harm,magn,type="l",xlab="periodo",ylab="magnitud")
  #find the place with maximum magnitude
  indices<-0
  for (i in 2:(nyquist-2))
  {
    if ((magn[i-1] < magn[i]) & (magn[i] > magn[i+1]))
    {
    indices<-c(indices,i)
    }
  }
  #clip out the first element
  indices<-indices[2:length(indices)]
  #form a vector of maxima in magnitude and period
  peaks<-magn[indices]
  periods<-harm[indices]
  peak.data<-data.frame(periods,peaks,indices)
  peak.data<-peak.data[order(peak.data$peaks,decreasing=T),]
  ratio<-peak.data[1,"peaks"]/peak.data[2,"peaks"]
  return(list(harmonics=harm,magnitude=magn,mat=peak.data,ratio=ratio))
}