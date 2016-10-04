f.peak<-function(v)
{
  #This function takes a vector and finds the first peak
  #and returns the position of the peak.
  #Specially designed for finding periods using correlograms
  
  #find the lenght of the vector
  pos<-NA
  n<-length(v[[1]])
  for (i in 2:(n-1))
  {
    if ((v[[1]][i-1] < v[[1]][i]) & (v[[1]][i] > v[[1]][i+1])) 
    {
      pos<-i
      break
    }
  }
  return(pos)
}