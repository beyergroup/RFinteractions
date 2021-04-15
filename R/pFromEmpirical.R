######################################################################
### get P values from a 2D empirical null distribution
### x=observed x values
### y=observed y values
### xn=x values under the null
### yn=y values under the null
### ...=additional arguments to be passed to kde2d
######################################################################
pFromEmpirical = function(x,y,xn,yn,...){
  require(MASS)
  rgy = range(c(y,yn))
  rgx = range(c(x,xn))
  kd = kde2d(xn,yn,lims=c(rgx,rgy),...)
  z = numeric(length(x))
  names(z) = names(x)
  for(i in 1:length(z)){
    z[i] = kd$z[which.min(abs(x[i]-kd$x)),which.min(abs(y[i]-kd$y))]
  }
  p = numeric(length(x))
  names(p) = names(x)
  for(i in 1:length(p)) p[i] = sum(kd$z[kd$z<z[i]])/sum(kd$z)
  return(p)
}