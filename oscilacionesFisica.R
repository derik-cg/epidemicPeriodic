#normal time cosine plot
time<-seq(from=0,to=30,by=0.1)
cosine<-cos(time)
plot(time,cosine)
cos.ts<-ts(cosine,start=0,frequency = 365)
fft.cos<-fft(cosine)
nyquist<-floor(length(cosine)/2)
magn<-sqrt(Re(fft.cos[1:nyquist])^2+Im(fft.cos[1:nyquist])^2)*2/nyquist
armonicos

#add the amplitude
R<-2 #amplitude
cos.amp<-R*cos(time)
plot(spec.cos$freq,spec.cos$spec,xlim=c(0,0.05))

#add the angular frequency
omega<-0.5 #angular frequency
cos.om<-cos(omega*time)
plot(time,cos.om)
# angular frequency of 0.5 doubles the time needed to
# reach the period. Angular frequency is how many time units
# it takes to reach the period. In terms of the unit circle
# it measures how many time units it take to do whole turn.

# Add the pase
phi<-1 #phase
cos.ph<-cos(time+phi)
plot(time,cos.ph)
#the phase moves the wave to the right or left, without
#altering the period or amplitude. Just stars later than
#by default.
