library(deSolve)

# This function returns the derivative of each species abundance. r is a vector and A is a matrix.
# The function %*% is R's shortcut for matrix multiplication: A%*%n is a vector whose i-th component is Sum_j A_ij * n_j
TFM=function(time,abundances,parameters){
 f=parameters$f
 h=parameters$h
 T=parameters$T
 B=f*t(T)%*%(h/(T%*%(f*abundances)))
 dndt=abundances*(B-1)
 dndt[abundances<0]=0
 return(list(dndt))
}

#now use the ODE solver in package deSolve.  
output_TFM=function(N0,time_range=1:5000,func=TFM,parms=list(f=f,h=h,T=T)) ode(N0,time_range,func=func,parms=parms)

