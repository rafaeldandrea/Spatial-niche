ClusterFit=function(N,S=NULL,S0=NULL,print=FALSE,plot=TRUE,barcode=TRUE){
	P=function(S0,S,k) choose(S,k)/choose(S0,k)/(S/(S0-S+1))
	
	if(is.null(S) | is.null(S0)){ 
		y=sign(N)
		S0=length(y)
		S=sum(y)
	}else{ y=sample(c(rep(1,S),rep(0,S0)))}

	s=paste(y,collapse='')
	z=unlist(lapply(strsplit(s,split='0')[[1]],nchar))
	w=plyr::count(z[z>0]) 
	
	p=P(S0,S,w$x)
	O=w$freq
	E=sum(O)*p/sum(p)
	
	G=2*sum(O*log(O/E)); k=nrow(w)-1
	pval=1-pchisq(G,df=k)
	zscore=(G-k)/sqrt(2*k)

	if(nrow(w)<=4 & S>80 & max(w$x)>5){ zscore=Inf; pval=0}		## When the number of strings of different sizes is too small the G-test fails. 
																## But this is a highly nonrandom situation. The null distribution of this number for
																## S > 80 excludes nrow(w)<=4.
	
	if(barcode==TRUE){ dev.new(width=9,height=3); par(mar=c(2,2,2,2)); plot(rep(1,S0),t='h',col=2-y,lwd=4,ylim=c(0,1),xaxt='n',yaxt='n',ann=FALSE)}
	
	if(plot==TRUE){
		# dev.new()
		with(w,plot(x,O/sum(O),las=1,pch=20,ylim=c(0,1),xlab='Length of the unbroken chain',ylab='Probability'))
		lines(w$x,E/sum(E),lwd=2,col=2)
	}	
	
	if(print==TRUE) print(paste('G-statistic = ',round(G,1),'   p-value = ',round(pval,3)))

	return(list(lattice.size=S0,extant.species=S,chain.length=w$x,counts=O,expected=E,G=G,p.value=pval,z.score=zscore,df=k))
}

BinForClusterFit=function(dat,nbins){
	if(!'trait'%in%names(dat)) stop('Data does not have a "trait" column')
	dat=subset(dat,!is.na(trait) & is.finite(trait)); dat=dat[order(dat$trait),]
	intvl=findInterval(dat$trait,seq(min(dat$trait),max(dat$trait)+.01,l=nbins+1))
	dtf=data.frame(bin=intvl,N=dat$N)
	res=ddply(dtf,.(bin),function(v) sum(v$N)); names(res)=c('trait','N')
	return(res)
}

