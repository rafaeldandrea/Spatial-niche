## Functions for fitting James O'Dwyer's species abundance probability distribution 

library(dplyr)
library(pracma)
library(dgof)
library(zipfR)
library(dgof)
library(hypergeo)
library(orthopolynom)
## -------------------- James HPFQ series ---------------------------
JO_X=function(nbar,S,cd,co,rho,eta){
  x=S*co/(cd-co)
  uniroot(f=function(X){ 
      X-(x*genhypergeo(U=c(1,1,1+nbar*x+1/X),L=c(2,2+nbar*x),z=rho/eta*X,check_mod=FALSE))/
      ((1+nbar*x)*genhypergeo(U=c(1,1,1+nbar*x+1/X),L=c(2,1+nbar*x),z=rho/eta*X,check_mod=FALSE))
     },
     lower=1e-16,upper=1/10)$root
}

logPoch=function(z,n){
  lgamma(z+n)-lgamma(z)
}

PochRatio=function(X,z,n){
  exp(logPoch(1/X+z,n)-logPoch(z,n))
}

HypergeoSeriesTerm=function(k,X,nbar=100,x=S*co/(cd-co)){
  exp(
    2*logPoch(z=1,n=k)-
      logPoch(z=2,n=k)+
      k*log(nbar*X)-
      lgamma(1+k)
  )*PochRatio(X=X,z=1+x*nbar,n=k)
}

HypergeoSeriesApprox=function(k,X,nbar,x) sum(HypergeoSeriesTerm(0:k,X,nbar,x))

JO_pmf=function(nvec,X=NULL,nbar,S,cd,co,rho,eta,mode='exact'){ 
  x=S*co/(cd-co)
  if(is.null(X)) X=x/(1+x)/nbar
  if(mode=='exact'){ 
      normalization=1/(rho/eta*X*genhypergeo(U=c(1,1,1+nbar*x+1/X),L=c(2,1+nbar*x),z=rho/eta*X,check_mod=FALSE))
  }
  if(mode=='approx'){ 
      normalization=1/(rho/eta*X*HypergeoSeriesApprox(k=5e4,X=X,nbar=nbar,x=x))
  }
  f=1/nvec*(X*rho/eta)^nvec*PochRatio(X=X,z=1+x*nbar,n=nvec-1)
  
  return(f*normalization)
}

JO_pmf_ratio=function(X,n1,n2,x,nbar,rho,eta){
  if(n2>=n1) return(NA)
  n2/n1*PochRatio(X,1+x*nbar+n2,n1-n2)*(rho*X/eta)^(n1-n2)
}
  
JO_pmf_ratio_exact=function(X,x,nbar,rho,eta){
  mod=JO_pmf(nvec=c(1,nbar),X=X,nbar=nbar,S=S,cd=cd,co=co,rho=rho,eta=eta)
  return(mod[2]/mod[1])
}

Likelihood=function(X,data,nbar,S,cd,co,rho,eta,mode) sum(log(JO_pmf(data,X,nbar,S,cd,co,rho,eta,mode)))
  
## Cumulative distribution function
JO_cdf=function(X,nbar,S,cd,co,rho,eta,nmax=5e4,mode){
  n=seq(nmax)
  pmf=JO_pmf(n,X,nbar,S,cd,co,rho,eta,mode)
  cdf=cumsum(pmf)
  cum=sum(pmf)-cdf+pmf
  return(data.frame(x=n,pmf=pmf,cdf=cdf,cum=cum))
}

## Fit LS, get p-value from discrete Cramer-von Mises test 
## (implemented via library cvm.test)
JO_fit_cvm=function(n,nbar,S,cd,co,rho,eta,mode){ 
  
  ## Fit X parameter
  X=optimize(
    f=Likelihood,
    data=n,
    nbar=nbar,
    S=S,
    cd=cd,
    co=co,
    rho=rho,
    eta=eta,
    mode=mode,
    lower=0.8*1/nbar*S*co/(cd-co)/(1+S*co/(cd-co)),
    upper=1.2*1/nbar*S*co/(cd-co)/(1+S*co/(cd-co)),
    maximum=TRUE,
    tol=1e-5
  )$maximum
  
  # X=mean(1/(n/(S*co/(cd-co))+mean(n)))
  
  ## transform data into count data frame
  w=if(is.vector(n)) plyr::count(n) else n
  
  ## calculate ecdf given the data
  dtf=JO_cdf(X,nbar,S,cd,co,rho,eta,nmax=10*max(w$x),mode=mode)
  myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
  
  ## call function cvm.test from package dgof, 
  ## which performs the CVM test on discrete data/null distribution
  mod=dgof::cvm.test(with(w,rep(x,freq)),myfun)
  
  ## return the fitted parameter, the test statistic (the W2 distance) 
  ## and the p-value of the test
  return(list(dtf=dtf,parameter=X,statistic=mod$statistic,p.value=mod$p.value))
}

## -------------------- Log-series (LS), P(n) = -1/log(1-p)*p^n/n --------------------------------------------

## Probability mass function
LS_pmf=function(p,n) -1/log(1-p)*p^n/n

## Cumulative distribution function
LS_cdf=function(p,nmax=5e4){ 
	n=seq(nmax)
	pmf=LS_pmf(p,n)
	cdf=1+zipfR::Ibeta(p,n+1,1e-16)/log(1-p)
	cum=1-cdf+pmf
	return(data.frame(x=n,cdf=cdf,cum=cum))
}

## Maximum likelihood estimation of parameter p 
LS_p=function(mu) uniroot(f=function(p) -1/log(1-p)*p/(1-p)-mu,lower=1e-16,upper=1-1e-16)$root

## Fit LS, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
LS_fit=function(data){ 
	## transform data into count data frame
	w=if(is.vector(data)) plyr::count(data) else data
	
	## fit parameter p
	p=LS_p(weighted.mean(w$x,w$freq))
	
	## calculate ecdf given the data
	dtf=LS_cdf(p=p,nmax=10*max(w$x))
	myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
	
	## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
	mod=dgof::ks.test(with(w,rep(x,freq)),myfun)
	
	## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
	return(list(dtf=dtf,parameter=p,statistic=mod$statistic,p.value=mod$p.value))
}


## Fit LS, get p-value from discrete Cramer-von Mises test 
## (implemented via library cvm.test)
LS_fit_cvm=function(data){ 
  ## transform data into count data frame
  w=if(is.vector(data)) plyr::count(data) else data
  
  ## fit parameter p
  p=LS_p(weighted.mean(w$x,w$freq))
  
  ## calculate ecdf given the data
  dtf=LS_cdf(p=p,nmax=10*max(w$x))
  myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
  
  ## call function cvm.test from package dgof, 
  ## which performs the CVM test on discrete data/null distribution
  mod=dgof::cvm.test(with(w,rep(x,freq)),myfun)
  
  ## return the fitted parameter, the test statistic (the W2 distance) 
  ## and the p-value of the test
  return(list(dtf=dtf,parameter=p,statistic=mod$statistic,p.value=mod$p.value))
}



## ----------------- Inverse quadratic distribution (IQ), P(n) ~ 1/n/(1+a*n) -----------------------------------------

## Probability mass function 
IQ_pmf=function(a,n){ 
	k=1/(digamma(1+1/a)-digamma(1))
	k*1/(n*(1+a*n))
}

## Cumulative distribution function
IQ_cdf=function(a,nmax=5e4){
	n=seq(nmax)
	pmf=IQ_pmf(a,n)
	cdf=cumsum(pmf)
	cum=sum(pmf)-cdf+pmf
	return(data.frame(x=n,cdf=cdf,cum=cum))
}

## Maximum likelihood estimation of parameter a
IQ_a=function(data,amin,amax){
	w=if(is.vector(data)) plyr::count(data) else data
	x=w$x; y=w$freq
	a=uniroot(f=function(a) 1/a^2*pracma::psi(1,1+1/a)/(digamma(1+1/a)-digamma(1))*sum(y)-sum(x*y/(1+a*x)),lower=amin,upper=amax)$root
	return(c(a=a))
}

## Fit IQ, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
IQ_fit=function(data,amin=1e-16,amax=1-1e-16){ 
	## fit parameter a
	a=IQ_a(data,amin,amax)
	
	## transform data into count data frame
	w=if(is.vector(data)) plyr::count(data) else data
	
	## calculate ecdf given the data
	dtf=IQ_cdf(a=a,nmax=10*max(w$x))
	myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
	
	## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
	mod=dgof::ks.test(with(w,rep(x,freq)),myfun)
	
	## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
	return(list(dtf=dtf,parameter=a,statistic=mod$statistic,p.value=mod$p.value))
}


IQ_KS=function(data,amin=1e-16,amax=1-1e-16,numnulls=100){ 
	## fit parameter a
	a=IQ_a(data,amin,amax)
	
	## transform data into count data frame
	w=if(is.vector(data)) plyr::count(data) else data
	
	## calculate ecdf given the data
	dtf=IQ_cdf(a=a,nmax=10*max(w$x))
	myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
	
	## nulls
	n=seq(1e5)
	nulls=sapply(seq(numnulls),function(null){
		set.seed(null)
		sample(n,size=sum(w$freq),replace=TRUE,prob=IQ_pmf(a,n))
	})
	
	
	## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
	x=dgof::ks.test(rep(w$x,w$freq),myfun)$statistic
	z=apply(nulls,2,function(null) dgof::ks.test(null,myfun)$statistic)
	
	return(list(obs.D=x,exp.D=z))
}


## ------------- Poisson distribution (Pois), P(n) = exp(-lambda)*lambda^n/n! ---------------------------------------

## Probability mass function
Pois_pmf=function(n) dpois(n,mean(n)) 

## Cumulative distribution function
Pois_cdf=function(lambda,nmax=5e4){
	n=seq(nmax)
	pmf=dpois(n,lambda)
	cdf=zipfR::Rgamma(lambda,n+1)
	cum=sum(pmf)-cdf+pmf
	return(data.frame(x=n,cdf=cdf,cum=cum))
}

## Maximum likelihood estimation of parameter lambda
Pois_lambda=function(data){
	w=if(is.vector(data)) plyr::count(data) else data
	x=w$x; y=w$freq
	return(c(lambda=weighted.mean(x,y)))
}

## Fit Poisson, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
Pois_fit=function(data){ 
	## fit parameter lambda
	lambda=Pois_lambda(data)
	
	## transform data into count data frame
	w=if(is.vector(data)) plyr::count(data) else data
	
	## calculate ecdf given the data
	dtf=Pois_cdf(lambda=lambda,nmax=max(w$x))
	myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
	
	## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
	mod=dgof::ks.test(with(w,rep(x,freq)),myfun)
	
	## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
	return(list(dtf=dtf,parameter=lambda,statistic=mod$statistic,p.value=mod$p.value))
}


## ----------------- Kessler-Shnerb distribution (KeSh), P(n) ~ 1/n/(1+a*n)^(1+b) -----------------------------------------

## Normalization
Norm=function(E,sigma2,mu,nmax=5e4) 1/sum(1/seq(nmax)/(1+2*E/sigma2*seq(nmax))^(1+mu/E))

## Probability mass function 
KeSh_pmf=function(n,E,b,d,nmax=5e4){ 	## d=eta, b=epsilon*rho*NR/NS
	sigma2=b+d; mu=b-d
	p=1/n/(1+2*E/sigma2*n)^(1+mu/E)
	k=Norm(E,sigma2,mu,nmax)
	k*p
}

## Cumulative distribution function
KeSh_cdf=function(E,b,d,nmax=5e4){
	n=seq(nmax)
	pmf=Kesh_pmf(n,E,b,d,nmax)
	cdf=cumsum(pmf)
	cum=sum(pmf)-cdf+pmf
	return(data.frame(x=n,cdf=cdf,cum=cum))
}

## Maximum likelihood estimation of parameter E
KeSh_E=function(data,Emin,Emax){
	w=if(is.vector(data)) plyr::count(data) else data
	x=w$x; y=w$freq
	a=uniroot(f=function(a) 1/a^2*pracma::psi(1,1+1/a)/(digamma(1+1/a)-digamma(1))*sum(y)-sum(x*y/(1+a*x)),lower=amin,upper=amax)$root
	return(c(a=a))
}

## Fit IQ, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
KeSh_fit=function(data,amin=1e-16,amax=1-1e-16){ 
	## fit parameter a
	a=IQ_a(data,amin,amax)
	
	## transform data into count data frame
	w=if(is.vector(data)) plyr::count(data) else data
	
	## calculate ecdf given the data
	dtf=IQ_cdf(a=a,nmax=10*max(w$x))
	myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))
	
	## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
	mod=dgof::ks.test(with(w,rep(x,freq)),myfun)
	
	## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
	return(list(dtf=dtf,parameter=a,statistic=mod$statistic,p.value=mod$p.value))
}

## ----------------- Negative Binomial distribution (NB), P(n) ~ Γ(x+size)/(Γ(size) x!) prob^size (1-prob)^x -------------------------------------
## mu = size * prob/(1 - prob)
NB_fit=function(data,lambda,eta){
	## transform data into count data frame
	w=if(is.vector(data)) plyr::count(data) else data
	
	## fit the negative binomial
	nbfit=sads::fitnbinom(data,trunc=NULL,start.value=c(prob=lambda/eta,mu=mean(data)))
	mu=as.numeric(nbfit@details$par['mu'])
	size=as.numeric(nbfit@details$par['size'])

	## calculate ecdf given the data
	dtf=data.frame(x=0:(2*max(data)),cdf=pnbinom(0:(2*max(data)),size=size,mu=mu)) %>% mutate(cum=1-cdf)
	myfun=stepfun(w$x,c(0,subset(dtf,x%in%w$x)$cdf))

	## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
	mod=dgof::ks.test(with(w,rep(x,freq)),myfun)
	
	## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
	return(
	  list(
	    dtf=dtf,
	    parameter=c(size=size,mu=mu),
	    statistic=mod$statistic,
	    p.value=mod$p.value
    )
  )
}



## ---------------- Wrapper ----------------
fitSAD=function(data,dbn,method='KS',...){
	stopifnot(dbn%in%c('LS','IQ','Pois','NB'))
  if(method=='CVM' & dbn=='LS'){ 
    return(LS_fit_cvm(data))
	}else if(dbn=='LS') return(LS_fit(data))
	if(dbn=='IQ') return(IQ_fit(data,...))
	if(dbn=='Pois') return(Pois_fit(data))	
	if(dbn=='NB') return(NB_fit(data,...))
}


## ===============================================================================================
## ================================== DISCARDED CODE =============================================
## ===============================================================================================

## Formats data for KS test	
# wcc=function(consumers){ 
	# wc=plyr::count(consumers) %>% 
		# merge(.,data.frame(x=seq(max(.$x))),by='x',all.y=TRUE) %>% 
		# mutate(freq=ifelse(is.na(freq),0,freq)) %>% 
		# mutate(cum=sum(freq)-cumsum(freq)+freq,cdf=cumsum(freq)/sum(freq))
	# return(data.frame(x=wc$x,cum=wc$cum,cdf=wc$cdf))
# }

## Chi-squared and G index for the log-series fit
# lschisq=function(consumers){ 
	# wc=plyr::count(consumers) %>% 
		# merge(.,data.frame(x=seq(max(.$x))),by='x',all.y=TRUE) %>% 
		# mutate(freq=ifelse(is.na(freq),0,freq)) %>% 
		# mutate(cum=sum(freq)-cumsum(freq)+freq) %>%
		# mutate(pred=length(consumers)*cdfls(seq(max(consumers)),p=p.pred))		
	# return(data.frame(chisq=with(wc,sum((pred-cum)^2/pred)),G=with(wc,2*sum(cum*log(cum/pred)))))
# }


# IQ_fit.a=function(data,nmax){
	# mod=optimize(IQ_ks.stat,c(1e-5,1),data=data,nmax=nmax)
	# a=mod$minimum
	# KS=mod$objective
	# return(data.frame(a=a,KS=KS))
# }


# IQ_ks=function(a,data,nmax=5e4){
	# if(is.data.frame(data)) o=data
	# if(is.vector(data)) o=wcc(data)
	# e=IQ_cdf(a,nmax=max(o$x))
	# return(rename(merge(o,e,by='x'),o=cdf.x,e=cdf.y))
# }


# IQ_ks.stat=function(a,data,nmax=5e4){ 
	# K=IQ_ks(a,data,nmax)
	# max(abs(K$o-K$e))
# }


# IQ_plot.fit=function(a,data,main=FALSE,legend=FALSE,nmax=5e4){
	# wc=mutate(plyr::count(data),cum=sum(freq)-cumsum(freq)+freq,cdf=cumsum(freq))
	# plot(0,t='n',xlim=c(min(wc$x),max(wc$x)),ylim=c(1,length(data)),xlab='Abundance',ylab='Cumulative No. Species',log='x')
	# if(main) title(main='Consumers')
	# w=mutate(merge(wc,IQ_cdf(a,nmax),by='x',all.x=TRUE),cum.y=length(data)*cum.y)	
	# with(w,{
		# lines(x,cum.x,lwd=2)
		# lines(x,cum.y,lwd=2,col=red)
		# lines(x,length(data)*cdfls(x,p=fit.logser(data)$summary['X']),col=blue,lwd=2)
		# if(legend) legend('topright',bty='n',col=c(red,blue),lwd=2,legend=c('P(n) ~ (1/n)*(a*n+1)^(-1)','Log-series'),cex=1.5)
	# })
# }


# IQ_pval=function(data,numnulls,nmax=5e4){
	# mod=optimize(IQ_ks.stat,c(1e-5,1),data=data,nmax=nmax)
	# a=mod$minimum
	# KS=mod$objective
	# set.seed(0)
	# nullKS=sapply(seq(numnulls),function(k){ 
		# ndata=sample(nmax,size=if(is.data.frame(data)) data$cum[1] else length(data),replace=TRUE,prob=1/seq(nmax)/(1+a*seq(nmax))) %>%
			# plyr::count(.) %>%
			# mutate(cdf=cumsum(freq)/sum(freq),cum=sum(freq)-cumsum(freq)+freq)
		# return(IQ_ks.stat(a,ndata,nmax))
	# })
	# return(c(a=a,KS=KS,p.value=sum(nullKS>=KS)/numnulls))
# }

## Fits Poisson to data
# PoisKS=function(n){
	# if(is.vector(n)){
		# wc=plyr::count(n) %>% 
			# merge(data.frame(x=seq(min(n),max(n))),.,all.x=TRUE) %>% 
			# mutate(freq=ifelse(is.na(freq),0,freq)) %>%
			# mutate(cdf=cumsum(freq)/sum(freq))
		# mu=mean(n)
	# }
	# if(is.data.frame(n)){
		# wc=n
		# mu=with(wc,weighted.mean(x,freq))
	# }
	# x=min(wc$x):max(wc$x)
	# pred=ncol(wc)*(zipfR::Rgamma(1+x,mu)+dpois(x,mu))
	# obs=wc$cdf
	# KS=max(abs(obs-pred))
	# return(KS)	
# }

## Find p-value for Poisson fit based on Kolmogorov-Smirnov distance
# KS_Pvalue=function(n){
	# KS=PoisKS(n)
	# mu=if(is.vector(n)) mean(n) else if(is.data.frame(n)) with(n,weighted.mean(x,freq))
	# S=if(is.vector(n)) length(n) else if(is.data.frame(n)) nrow(n)
	# nulls=sapply(seq(numnulls),function(.) sample(1e4,size=S,replace=TRUE,prob=dpois(seq(1e4),mu)))
	# nKS=apply(nulls,2,function(.) PoisKS(.))
	# p.value=sum(nKS>=KS)/numnulls
	# return(c(KS=KS,p.value=p.value))
# }

