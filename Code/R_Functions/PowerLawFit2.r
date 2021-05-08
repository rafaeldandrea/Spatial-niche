## Require library VGAM for the zeta function
library(VGAM)

## Maximum likelihood power law fit for the tail of the progeny, based on the Kolmogorov-Smirnov (KS) method. 
## The input is a plyr::count data frame, the vector of candidate thresholds for the power law, and the vector of candidate exponents.
## The output is a list with the estimated exponent and threshold, the minimum KS distance, and the log-likelihood.
pl.fit=function(dtf,xmins,xmax=max(dtf$x),cap=max(dtf$x),alphas=seq(1.45,1.60,by=.01)){	
	
	## Trim to values below potential exponential cutoff.
	w=subset(dtf,x<=xmax)
	
	## Hurwitz zeta function, calculated as (sum_{k=1,infty} - sum_{k=1,x0-1}) k^-alpha.
	## The Hurwitz is the normalization constant in the discrete pl pmf. Here I'm using vectorized notation to save computation effort: 
	## sum(f(seq(x0-1))) = sum(f(seq(x0))) - f(x0); cumsum(f(seq(x0-1))) = cumsum(f(seq(x0))) - f(seq(x0)).
	H=t(sapply(alphas,function(a) VGAM::zeta(a)-(cumsum(seq(max(xmins))^-a)-seq(max(xmins))^-a)))
	if(length(xmins)>1) H=H[,xmins]	else H=matrix(H,length(alphas),1)	
	
	## Log-likelihood.
	l1=outer(alphas,xmins,Vectorize(function(a,x0) with(subset(w,x>=x0),-a*sum(freq*log(x)))))
	l2=-t(t(log(H))*sapply(xmins,function(x0) sum(subset(w,x>=x0)$freq)))
	logL=l1+l2																		
	
	## Fit alpha by maximizing logL for each candidate xmin.
	alphafit=alphas[apply(logL,2,which.max)]										
	
	## Kolmogorov-Smirnov measure.
	KS=sapply(xmins,function(x0){													
		tmp=subset(w,x>=x0 & x<=cap)
		a=alphafit[match(x0,xmins)]
		
		## power law CDF
		theo=cumsum(as.integer(min(tmp$x):max(tmp$x))^-a)/H[match(a,alphas),match(x0,xmins)]				
		theo=theo[as.integer(min(tmp$x):max(tmp$x))%in%tmp$x]
		
		## Empirical CDF.
		emp=with(tmp,cumsum(freq/sum(freq)))												
		
		return(max(abs(emp-theo)))													
	})
	
	## Estimates are those that minimze the KS measure.
	D=min(KS); alphahat=alphafit[which.min(KS)]; xminhat=xmins[which.min(KS)]; 	
	log.likelihood=logL[match(alphahat,alphas),match(xminhat,xmins)]
	
	## Standard error of the exponent
	stder=(alphahat-1)/sqrt(sum(subset(w,x>=xminhat & x<=cap)$freq)) 	
	
	return(list(alphahat=alphahat,xminhat=xminhat,stder=stder,D=D,log.likelihood=log.likelihood))
}

## Fits a power law using pl.fit and calculates the empirical distribution function.
## The input is a plyr::count data frame and the vector of candidate thresholds for the power law.
## The output is the estimated parms of the power law and the fitted empirical distribution function.
cdfpl=function(w,xmins){
	
	plfit=pl.fit(w,xmins,cap=1e6)
	alpha=plfit$alphahat; xmin=plfit$xminhat; stder=plfit$stder
	
	## trim data
	x=subset(w,x>=xmin)$x; freq=subset(w,x>=xmin)$freq
	
	## normalization factor
	norm=if(xmin==1) 0 else sum((1:(xmin-1))^-alpha)
	
	## empirical distribution function
	vec=(zeta(alpha)-cumsum(seq(max(x))^-alpha))/(zeta(alpha)-norm)		
	
	return(list(xmin=xmin,alpha=alpha,stder=stder,log.likelihood=plfit$log.likelihood,y=sum(freq)*vec[x]))
}