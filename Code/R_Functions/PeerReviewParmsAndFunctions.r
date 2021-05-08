library(plyr)

editorial.methods=c('none','tiebreak','overrule','blacklist.50','blacklist.90','match.quality','revision') 
method.names=c('None','Tiebreak','Bypass refs','Blacklist (p = 0.5)','Blacklist (p = 0.9)','Match refs','Revision')

N=1000																				## Total number of scientists
sigma.author=10																		## std dev of distribution of proficiency among authors
sigma.quality=5																		## std dev of dbn of quality of submitted ms by a single author
sigma.lognormal=.5																	## std dev of quality dbn of authors & submissions in lognormal scenario
lambda=.1																			## Controls ``stickiness'' of the moving standard
steps=500																			## Total number of review cycles

## --------------------- Functions ---------------------------------------------------
Decision=function(ref,q.sub,q.min,ref.type,Q.author,Qmin.submit,Q.min.fixed.standards,dbn){
	if(ref.type[ref]=='moving.standards') return(1*(q.sub>=q.min))
	if(ref.type[ref]=='indifferent') return(1*(q.sub>=Qmin.submit & q.sub<=Q.author[ref]))
	if(ref.type[ref]=='random') return(round(runif(1)))
	if(ref.type[ref]=='fixed.standards') return(1*(q.sub>=Q.min.fixed.standards))
	if(ref.type[ref]=='conscientious') return(1*(q.sub>=Q.min.fixed.standards & q.sub<=Q.author[ref]))
	if(ref.type[ref]=='narcissist' & dbn=='normal') return(1*(abs(q.sub-Q.author[ref])<=1.96*sigma.quality))
	if(ref.type[ref]=='narcissist' & dbn=='lognormal') return(1*(abs(log(q.sub)-log(Q.author[ref]))<=1.96*sigma.lognormal))
}

PeerReview=function(fr=0,fc=0,frnd=0,ff=0,fd=0,fn=0,run=0,method='none',num.refs=2,dbn='normal'){
	stopifnot(method%in%editorial.methods)
	if(fc+fr+frnd+ff+fd+fn!=1) stop('Ref types do not add to 1')
	set.seed(run)	
	
	if(dbn=='normal'){ 
		Q.author=rnorm(N,mean=100,sd=sigma.author) 											## Proficiency is normally distributed
		Qmin.submit=90																		## Min quality accepted by indifferent selfish refs,
																							## corresponds to the 18th percentile of submissions.
		Q.min.fixed.standards=qnorm(.9,mean=100,sd=sqrt(sigma.author^2+sigma.quality^2))	## Min quality accepted by fixed-standard honest refs
	}
	if(dbn=='lognormal'){
		Q.author=rlnorm(N,mean=log(100),sd=sigma.lognormal) 										## Proficiency is lognormal-distributed
		Qmin.submit=qlnorm(pnorm(90,mean=100,sd=sqrt(sigma.author^2+sigma.quality^2)),mean=log(100),sd=sqrt(2)*sigma.lognormal)	
																					## Min quality accepted by indifferent selfish refs,
																					## corresponds to the 18th percentile of submissions.
		Q.min.fixed.standards=qlnorm(.9,mean=log(100),sd=sqrt(2)*sigma.lognormal)	## Min quality accepted by fixed-standard honest refs
	}
	
	ref.type=sample(c(	rep('moving.standards',round(fc*N)),
						rep('indifferent',round(fr*N)),
						rep('random',round(frnd*N)),
						rep('fixed.standards',round(ff*N)),
						rep('conscientious',round(fd*N)),
						rep('narcissist',round(fn*N))
	))

	simtime=0
	Q.accept.avg=Q.submit.avg=Rejected=A.list=R.list=NULL
	Q.min=M=Q.accept=0
	while(simtime<steps){
		simtime=simtime+1
		if(simtime==1) authors=sample(N,size=N/2) else authors=setdiff(seq(N),authors)
		if(dbn=='normal') Q.submitted=rnorm(length(authors),mean=sample(Q.author[authors]),sd=sigma.quality)
		if(dbn=='lognormal') Q.submitted=rlnorm(length(authors),meanlog=log(sample(Q.author[authors])),sdlog=sigma.lognormal)
		if(method=='match.quality'){ 
			refs=sapply(seq_along(authors),function(j) sample(seq(N)[-authors[j]],size=3,prob=dnorm(Q.author[-authors[j]]-Q.submitted[j],mean=5,sd=5)))
		}else refs=sapply(authors,function(author) sample(seq(N)[-author],size=3))
		accepted=apply(rbind(refs,Q.submitted),2,function(v){
			ref=as.integer(v[1:3]); q.sub=v[4]
			recs=sapply(ref,function(r) Decision(r,q.sub,Q.min,ref.type,Q.author,Qmin.submit,Q.min.fixed.standards,dbn))
			prob.accept=if(num.refs==2 & method!='tiebreak') mean(recs[1:2]) else median(recs)
			if(method=='revision' & prob.accept<1){
				q.sub=q.sub+abs(rnorm(1,mean=0,sd=sigma.quality))
				recs=sapply(ref,function(r) Decision(r,q.sub,Q.min,ref.type,Q.author,Qmin.submit,Q.min.fixed.standards,dbn))
				prob.accept=if(num.refs==2) mean(recs[1:2]) else median(recs)
			}
			return(sample(c(1,0),1,prob=c(prob.accept,1-prob.accept)))
		})		
		if(method=='overrule'){ 
			accepted[Q.submitted>quantile(Q.submitted,.9)]=1
			# accepted[Q.submitted<quantile(Q.submitted,.5)]=0
		}
		M=lambda*M+(1-lambda)*mean(Q.accept)
		Q.min=M
		if(any(accepted==1)){ 
			Q.accept=Q.submitted[accepted==1]
			if(simtime>=.2*steps) A.list=c(A.list,Q.accept)			## Lists the quality of every accepted paper
		}
		if(any(accepted==0)){ 
			Q.reject=Q.submitted[accepted==0]
			if(simtime>=.2*steps) R.list=c(R.list,Q.reject)			## Lists the quality of every rejected paper
		}
		Q.accept.avg=c(Q.accept.avg,mean(Q.accept))					## Average quality of accepted manuscripts in the round
		Q.submit.avg=c(Q.submit.avg,mean(Q.submitted))				## Average quality of submitted manuscripts in the round
		Rejected=c(Rejected,sum(accepted==0)/length(authors))		## Percentage of total submitted manuscripts that were rejected in the round
		
	}
	return(list(
		data=data.frame(Q.accept.avg=Q.accept.avg,Q.submit.avg=Q.submit.avg,Rejected=Rejected),
		A.list=A.list,
		R.list=R.list
	))
}

PeerReviewBlacklist=function(fr,fc,p0,career,run=0,Qa=NULL,num.refs=2,dbn='normal'){
	time.factor=2				## Converts time from years to review cycles. Eg time.factor = 2 means 2 review cycles per year
	career=career*time.factor	## This converts the time in years which is fed to the function to the time in review cycles.
	
	stopifnot(fc+fr==1)
	
	if(is.null(Qa)) Qa=c(129,114,111,109,107,105,104,103,101,100)[match(fr,c(.001,1:9/10))]	## Those values were numerically obtained by running PeerReview(fr) 
	
	set.seed(run)	
	ref.type=sample(c(rep('moving.standards',round(fc*N)),rep('indifferent',round(fr*N))))

	if(dbn=='normal'){ 
		Q.author=rnorm(N,mean=100,sd=sigma.author) 											## Proficiency is normally distributed
		Qmin.submit=90																		## Min quality accepted by indifferent selfish refs,
																							## corresponds to the 18th percentile of submissions.
		Q.min.fixed.standards=qnorm(.9,mean=100,sd=sqrt(sigma.author^2+sigma.quality^2))	## Min quality accepted by fixed-standard honest refs
	}
	if(dbn=='lognormal'){
		Q.author=rlnorm(N,mean=log(100),sd=sigma.lognormal) 										## Proficiency is lognormal-distributed
		Qmin.submit=qlnorm(pnorm(90,mean=100,sd=sqrt(sigma.author^2+sigma.quality^2)),mean=log(100),sd=sqrt(2)*sigma.lognormal)	
																					## Min quality accepted by indifferent selfish refs,
																					## corresponds to the 18th percentile of submissions.
		Q.min.fixed.standards=qlnorm(.9,mean=log(100),sd=sqrt(2)*sigma.lognormal)	## Min quality accepted by fixed-standard honest refs
	}
	
	simtime=0; age=sample(rep(seq(career),each=N/career)); age=c(age,sample(career,N-length(age)))
	Q.accept.avg=Q.submit.avg=Rejected=A.list=R.list=NULL
	Q.min=M=Q.accept=0
	
	tally=data.frame(ref=seq(N),n=0,k=0)											## keeps score of the number of reviews (n) and disagreements (k)
	mat=Blacklist(fc=fc,fr=fr,Qa=Qa,Qmin.submit,Q.min.fixed.standards)				## lists the probs that a ref is indifferent based on n and k
	
	blacklist=integer(0)
	reftrack=refrecord=NULL
	while(simtime<steps){
		simtime=simtime+1
		
		## submission stage
		if(simtime==1) authors=sample(N,size=N/2) else authors=setdiff(seq(N),authors)
		if(dbn=='normal') Q.submitted=rnorm(length(authors),mean=sample(Q.author[authors]),sd=sigma.quality)
		if(dbn=='lognormal') Q.submitted=rlnorm(length(authors),meanlog=log(sample(Q.author[authors])),sdlog=sigma.lognormal)
		
		## choosing refs
		refs=sapply(authors,function(author) sample(seq(N)[-author],size=num.refs,prob=ifelse(seq(N)[-author]%in%blacklist,0,1)))
		reftrack=c(reftrack,as.numeric(refs))
		
		## reviewing process
		matreviews=sapply(seq(ncol(refs)),function(j){	## reviews for each manuscript
			q.sub=as.numeric(Q.submitted[j])
			msrefs=apply(refs,1,function(v) as.integer(v[j]))
			recommendations=sapply(msrefs,function(ref) Decision(ref,q.sub,Q.min,ref.type,Q.author,Qmin.submit,Q.min.fixed.standards,dbn))
		}) 
		if(num.refs==2) reviews=apply(matreviews,2,function(revs) mean(revs))
		if(num.refs==3) reviews=apply(matreviews,2,function(revs) median(revs))
		
		## accepting and rejecting submitted papers 
		accepted=sapply(reviews,function(p) sample(c(1,0),1,prob=c(p,1-p)))
		M=lambda*M+(1-lambda)*mean(Q.accept)
		Q.min=M
		
		if(any(accepted==1)){ 
			Q.accept=Q.submitted[accepted==1]
			if(simtime>=.2*steps) A.list=c(A.list,Q.accept)
		}
		if(any(accepted==0)){ 
			Q.reject=Q.submitted[accepted==0]
			if(simtime>=.2*steps) R.list=c(R.list,Q.reject)
		}		
		Q.accept.avg=c(Q.accept.avg,mean(Q.accept))
		Q.submit.avg=c(Q.submit.avg,mean(Q.submitted))
		Rejected=c(Rejected,sum(accepted==0)/length(authors))
		
		## blacklisting refs 
		if(num.refs==2){
			disagreements=sapply(reviews,function(p) 1*(p==.5))
			res=ddply(data.frame(ref=c(refs[1,],refs[2,]),n=1,k=rep(disagreements,2)),.(ref),function(v) c(n=sum(v$n),k=sum(v$k)))
			ind=tally$ref%in%res$ref; tally$n[ind]=tally$n[ind]+res$n; tally$k[ind]=tally$k[ind]+res$k
			tmp=subset(merge(tally,mat,all.x=TRUE),!is.na(pR))
			blacklist=with(tmp,ref[pR>=p0])
		}
		if(num.refs==3){
			mark=apply(matreviews,2,function(p) 1*(p!=median(p)))
			res=ddply(data.frame(ref=as.numeric(t(refs)),n=1,k=as.numeric(t(mark))),.(ref),function(v) c(n=sum(v$n),k=sum(v$k)))
			ind=tally$ref%in%res$ref; tally$n[ind]=tally$n[ind]+res$n; tally$k[ind]=tally$k[ind]+res$k
			tmp=subset(merge(tally,mat,all.x=TRUE),!is.na(pR))
			blacklist=with(tmp,ref[pR>=p0])
		}
		
		## scholar turnover 
		age=age+1
		Q.author[age>career]=rnorm(sum(age>career),mean=100,sd=sigma.author)
		blacklist=setdiff(blacklist,seq(N)[age>career])
		tally$n[age>career]=0; tally$k[age>career]=0
		ref.type[age>career]=sample(c('moving.standards','indifferent'),size=sum(age>career),replace=TRUE,prob=c(fc,fr))
		retiring=which(age>career); refrecord=c(refrecord,sapply(retiring,function(x) sum(reftrack==x)))
		reftrack=setdiff(reftrack,retiring)
		age[age>career]=1		
	}
	refrecord=c(refrecord,plyr::count(reftrack)$freq)
	
	return(list(
		data=data.frame(Q.accept.avg=Q.accept.avg,Q.submit.avg=Q.submit.avg,Rejected=Rejected),
		A.list=A.list,
		R.list=R.list,
		ref.record=refrecord
	))
}

Blacklist=function(fc,fr,Qa,Qmin.submit,Q.min.fixed.standards){
	PA=function(b) dnorm(b,mean=100,sd=sigma.author) 
	PP=function(x) dnorm(x,mean=100,sd=sqrt(sigma.author^2+sigma.quality^2))
	
	I0=Vectorize(function(b) PA(b)*integrate(PP,lower=Qa,upper=b)$value)
	qrc=1-integrate(PA,lower=-Inf,upper=ifelse(fr>0,Qmin.submit,Q.min.fixed.standards))$value-integrate(I0,lower=Qa,upper=Inf)$value
	
	I2=Vectorize(function(b2,b1) PA(b2)*integrate(PP,lower=b1,upper=b2)$value)
	I1=Vectorize(function(b1) PA(b1)*integrate(I2,lower=b1,upper=Inf,b1=b1)$value)
	qrr=integrate(I1,lower=ifelse(fr>0,Qmin.submit,Q.min.fixed.standards),upper=Inf)$value
	
	qc=fr*qrc; qr=fr*qrr+fc*qrc
	res=ddply(merge(data.frame(n=1:50),data.frame(k=0:50)),.(n,k),function(v){
		n=v$n; k=v$k
		if(n>=k){ 
			PkRn=qr^k*(1-qr)^(n-k)
			PkCn=qc^k*(1-qc)^(n-k)
			return(1/(1+PkCn/PkRn*fc/fr))
		} else return(NA)
	}); names(res)[3]='pR'
	return(res)
}
