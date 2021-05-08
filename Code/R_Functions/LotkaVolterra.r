library(deSolve)

# This function returns the derivative of each species abundance. r is a vector and A is a matrix.
# The function %*% is R's shortcut for matrix multiplication: A%*%n is a vector whose i-th component is Sum_j A_ij * n_j
LV=function(time,abundances,parameters){
 r=parameters$r
 A=parameters$A
 mu=parameters$mu
 abundances[abundances<1e-8]=0
 dndt=abundances*(r-A%*%abundances)+mu
 # dndt[abundances<1e-16]=0
 return(list(dndt))
}

#now use the ODE solver in package deSolve.  
output_LV=function(N0,time_range=1:5000,func=LV,parms=list(r=r,A=A,mu=mu)) ode(N0,time_range,func=func,parms=parms)

if(FALSE){
	#number of species
	S=100

	#set the list of time values we'll output the numerical solution of the growth curve at
	time_range=1:5000		

	#trait values
	set.seed(0)						## sets the seed of the random number generator. Useful for repeatable results.
	trait=seq(0,1,length=S)		## generates S trait values between 0 and 1, in regular increasing order. 


	#initial abundances
	N0=rep(1,S)	## vector with S identical elements 1. To make them random, set them to runif(S)	

	#immigration rates
	mu=rep(0,S)					## Alternatives: mu=rep(0,S); mu=rep(.001,S); mu<-.001*runif(S)
	mu<-.01*runif(S)

	#intrinsic growth rates
	r=rep(1,S)						## vector with S identical elements 1. Alternatives: r=rep(1,S); r=trait*(1-trait)
	r=trait*(1-trait)
	
	#trait differences
	d=as.matrix(dist(trait,diag=TRUE,upper=TRUE))	## matrix with the absolute trait differences |x_i - x_j| between all pairs of species 

	#competition coefficients
	w=0.3								## defines width of competition curve. Alternatives: w=0.15; w=seq(.01,.2,length=S)	
	A=exp(-(d/w)^4)+0*rnorm(S^2,0,.2)	## competition matrix - monotonic function of trait differences;
										## the second term adds a random perturbation to each coefficient
	A[A<0]=0; A[A>1]=1
	#diag(A)=1.1
	A=A/sum(A)*1e1


	#extract the relevant quantities
	sim=output_LV(N0=N0,parms=list(r=r,A=A,mu=mu))
	times=c(50,100,500,5000)			## times at which R will plot snapshots of the community
	# times=c(5,10,50,100)
	abunds1=sim[times[1],-1]		## abundances of all species at a specified time
	abunds2=sim[times[2],-1]
	abunds3=sim[times[3],-1]	
	abunds4=sim[times[4],-1]

	#set plot parameters
	if(names(dev.cur())!='null device') dev.off(); dev.new(width=12,height=12)	## tells R to kill the current plot device and open a new one
	par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,0))								## sets the margins of the new plot device

	#plot final abundance by trait
	plot(trait,abunds1,type="h",pch=20,las=1,main=paste("Time = ",times[1]))
	plot(trait,abunds2,type="h",pch=20,las=1,main=paste("Time = ",times[2]))
	plot(trait,abunds3,type="h",pch=20,las=1,main=paste("Time = ",times[3]))
	plot(trait,abunds4,type="h",pch=20,las=1,main=paste("Time = ",times[4]))	
	title(outer=TRUE,xlab="Trait",ylab="Abundance",main="Lotka-Volterra multiple species",line=.9,cex.lab=1.3)
}
if(FALSE){
  N0=N0_keep=abunds4
  N0[1:(S/2)]=0
  mu[1:(S/2)]=0
  sim=output_LV(N0=N0,parms=list(r=r,A=A,mu=mu))
  abunds1=sim[times[1],-1]		## abundances of all species at a specified time
  abunds2=sim[times[2],-1]
  abunds3=sim[times[3],-1]	
  abunds4=sim[times[4],-1]
  
  #set plot parameters
  par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,0))								## sets the margins of the new plot device
  
  #plot final abundance by trait
  plot(trait,abunds1,type="h",pch=20,las=1,main=paste("Time = ",times[1]))
  plot(trait,abunds2,type="h",pch=20,las=1,main=paste("Time = ",times[2]))
  plot(trait,abunds3,type="h",pch=20,las=1,main=paste("Time = ",times[3]))
  plot(trait,abunds4,type="h",pch=20,las=1,main=paste("Time = ",times[4]))	
  title(outer=TRUE,xlab="Trait",ylab="Abundance",main="Lotka-Volterra multiple species",line=.9,cex.lab=1.3)
}