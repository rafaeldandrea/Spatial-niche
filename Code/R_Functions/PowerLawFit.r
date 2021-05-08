library(VGAM)   # zeta function

## MLE power law fit for the tail of the progeny, based on the KS method. 
## The input is a plyr::count dataframe, the vector of candidate thresholds for the power law, and the vector of candidate exponents.
## The output is a list with the estimated exponent and threshold and the minimum KS.
pl.fit=function(dtf,xmins,xmax=max(dtf$x),cap=max(dtf$x),alphas=seq(1.45,1.60,by=.01)){	
	w=subset(dtf,x<=xmax)															## Trim to values below potential exponential cutoff.
	H=t(sapply(alphas,function(a) VGAM::zeta(a)-(cumsum(seq(max(xmins))^-a)-seq(max(xmins))^-a)))	## Hurwitz zeta function, calculated as 
	if(length(xmins)>1) H=H[,xmins]	else H=matrix(H,length(alphas),1)						## (sum_{k=1,infty} - sum_{k=1,x0-1}) k^-alpha.
																							## The Hurwitz is the normalization constant in the discrete pl pmf.					
																							## Here I'm using vectorized notation to save computation effort: 
																							## sum(f(seq(x0-1))) = sum(f(seq(x0))) - f(x0);
																							## cumsum(f(seq(x0-1))) = cumsum(f(seq(x0))) - f(seq(x0)).
	l1=outer(alphas,xmins,Vectorize(function(a,x0) with(subset(w,x>=x0),-a*sum(freq*log(x)))))
	l2=-t(t(log(H))*sapply(xmins,function(x0) sum(subset(w,x>=x0)$freq)))
	logL=l1+l2																		## Log-likelihood.			
	alphafit=alphas[apply(logL,2,which.max)]										## Fit alpha by maximizing logL for each candidate xmin.
	KS=sapply(xmins,function(x0){													## Kolmogorov-Smirnov measure.
		tmp=subset(w,x>=x0 & x<=cap)
		a=alphafit[match(x0,xmins)]
		theo=cumsum(as.integer(min(tmp$x):max(tmp$x))^-a)/H[match(a,alphas),match(x0,xmins)]				## power law CDF
		theo=theo[as.integer(min(tmp$x):max(tmp$x))%in%tmp$x]
		emp=with(tmp,cumsum(freq/sum(freq)))												## Empirical CDF.
		return(max(abs(emp-theo)))													
	})
	D=min(KS); alphahat=alphafit[which.min(KS)]; xminhat=xmins[which.min(KS)]; 	## Estimates are those that minimze the KS measure.
	stder=(alphahat-1)/sqrt(sum(subset(w,x>=xminhat & x<=cap)$freq)) 	
	log.likelihood=logL[match(alphahat,alphas),match(xminhat,xmins)]
	return(list(alphahat=alphahat,stder=round(stder,3),xminhat=xminhat,D=D,log.likelihood=log.likelihood))
}

pl.fit.sample=function(input,xmins,xmax=max(dat$x),alphas=seq(1.01,2,by=.01),numnulls=10,samplesize=sum(dat$freq)/2){
	input=subset(input,x<=xmax)
	tmp=t(sapply(seq(numnulls),function(k){
		jor=with(input,sample(x,size=samplesize,replace=1,prob=freq))
		mod=pl.fit(plyr::count(jor),xmins=xmins,xmax=xmax,alphas=alphas)
		return(c(a=mod$alphahat,s=mod$xminhat))
	}))
	return(list(a=mean(tmp[,1]),sa=sd(tmp[,1]),x0=mean(tmp[,2]),sx0=sd(tmp[,2])))
}


