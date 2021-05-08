library(deSolve)

# This function returns the derivative of each species abundance. r is a vector and A is a matrix.
# The function %*% is R's shortcut for matrix multiplication: A%*%n is a vector whose i-th component is Sum_j A_ij * n_j
HCCM=function(time,abundances,parameters){
	f=parameters$f		## fecundity: number of seeds per individual per unit time
	s=parameters$s		## 0 < s < infty: steepness parameter in competitive hierarchy. Flat when s = 0, step function when s = infty
	mu=min(f)			## mortality (common to all species). Its value defines the lowest viable fecundity in the community.
	K=parameters$K		## community-wide carrying capacity: the total number of sites available for colonization
	A=outer(f,f,function(f1,f2) (f1+f2)*1/2*(1-tanh(s*(f2-f1))))
	n=abundances
	n[which(n<1e-16)]=0
	dndt=n*((f-mu)-A%*%n/K)
	
	return(list(dndt))
}

#now use the ODE solver in package deSolve.  
output_HCCM=function(N0,time_range=1:5000,func=HCCM,parms=list(f=f,s=s,K=K)) ode(N0,time_range,func=func,parms=parms)

