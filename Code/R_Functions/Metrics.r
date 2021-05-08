## Runs the metrics on the data

require(plyr)

Metrics=function(dat,k){
	stopifnot(all(c('N','trait')%in%names(dat)))
	dat=subset(dat,N>0); dat=dat[order(dat$trait),]
	trait=dat$trait; N=dat$N
	
	z=kmeans(rep(trait,N),centers=seq(min(trait),max(trait),l=k),iter.max=100)
	centroids=sapply(z$centers,function(x) trait[which.min(abs(trait-x))])
	
	f=unique(data.frame(trait=rep(trait,N),clus=z$cluster))
	dtf=merge(dat,f,all.x=TRUE)
	
	peaks=sort(ddply(dtf,.(clus),function(v) v$trait[which.max(v$N)])[,2])
	others=setdiff(trait,peaks)
	flankers=sapply(peaks,function(p) others[which.min(abs(others-p))])
	
	##	Correlation
	Npeaks=N[trait%in%peaks]; Nflankers=sapply(flankers,function(f) N[match(f,trait)])
	metric.cor=ifelse(length(Npeaks)>2,cor(Npeaks,Nflankers),NA)
	
	## Depression
	vicinity=ddply(subset(dtf,trait%in%others),.(clus),function(v) mean(v$N)); names(vicinity)[2]='meanN'
	datflankers=data.frame(clus=sort(unique(dtf$clus)),flanker=flankers)
	metric.dep=with(merge(vicinity,datflankers,all.x=TRUE),mean(flanker/meanN))
	
	## Gap Rao
	interv=others[others>peaks[1] & others<peaks[k]]
	ind=findInterval(interv,peaks)
	newdtf=data.frame(trait=interv,N=N[match(interv,trait)],ind=ind)
	y=plyr::count(newdtf$ind); names(y)=c('ind','freq')
	newdtf=subset(merge(newdtf,y,all.x=TRUE),freq>=3,select=-freq)
	metric.rao=mean(ddply(newdtf,.(ind),function(v){
		p=v$N/sum(v$N); d=as.matrix(dist(v$trait))
		as.numeric(p%*%(d%*%p)/2)
	})[,2])
	
	return(data.frame(depression=metric.dep,correlation=metric.cor,gaprao=metric.rao))
}

MetricsPam=function(dat,k){
	require(cluster)
	stopifnot(all(c('N','trait')%in%names(dat)))
	dat=subset(dat,N>0); dat=dat[order(dat$trait),]
	trait=dat$trait; N=dat$N
	
	z=clara(rep(trait,N),k=k,samples=200)
	f=unique(data.frame(trait=rep(trait,N),clus=z$clustering))
	dtf=merge(dat,f,all.x=TRUE)
	
	peaks=sort(z$medoids)
	others=setdiff(trait,peaks)
	flankers=sapply(peaks,function(p) others[which.min(abs(others-p))])
	
	##	Correlation
	Npeaks=N[trait%in%peaks]; Nflankers=sapply(flankers,function(f) N[match(f,trait)])
	metric.cor=ifelse(length(Npeaks)>2,cor(Npeaks,Nflankers),NA)
	
	## Depression
	vicinity=ddply(subset(dtf,trait%in%others),.(clus),function(v) mean(v$N)); names(vicinity)[2]='meanN'
	datflankers=data.frame(clus=sort(unique(dtf$clus)),flanker=flankers)
	metric.dep=with(merge(vicinity,datflankers,all.x=TRUE),mean(flanker/meanN))
	
	## Gap Rao
	interv=others[others>peaks[1] & others<peaks[k]]
	ind=findInterval(interv,peaks)
	newdtf=data.frame(trait=interv,N=N[match(interv,trait)],ind=ind)
	y=plyr::count(newdtf$ind); names(y)=c('ind','freq')
	newdtf=subset(merge(newdtf,y,all.x=TRUE),freq>=3,select=-freq)
	metric.rao=mean(ddply(newdtf,.(ind),function(v){
		p=v$N/sum(v$N); d=as.matrix(dist(v$trait))
		as.numeric(p%*%(d%*%p)/2)
	})[,2])
	
	return(data.frame(depression=metric.dep,correlation=metric.cor,gaprao=metric.rao))
}

CVTrend=function(dat){
	stopifnot(all(c('trait','N')%in%names(dat))) 
	dat=subset(dat,N>0); trait=dat[order(dat$N,decreasing=TRUE),]$trait
	cv=sapply(nosp<-length(trait):5,function(n){
		d=diff(sort(trait[1:n]))
		sd(d)/mean(d)
	})	
	return(data.frame(nosp=nosp,cv=cv))
}

## Rao's Quadratic Entropy and Laliberte Functional Dispersion
LitMetrics=function(dat,onedim=TRUE){
	dat=dat[with(dat,!is.na(trait)&N>0),]
	if(onedim) tr=dat$trait else tr=with(dat,cbind(SLA,leaf.area,max.height,seed.mass))
	N=dat$N
	
	##Rao's quadratic entropy
	dm=as.matrix(dist(tr,upper=TRUE,diag=TRUE))
	p=N/sum(N)
	Q=p%*%(dm%*%p)/2
	
	##Functional Dispersion
	centroid=as.numeric((t(tr)%*%N)/sum(N))
	if(onedim) z=abs(tr-centroid) else z=sqrt(rowSums((tr-centroid)^2))
	FDis=(z%*%N)/sum(N)
	
	##CV
	d=diff(sort(dat$trait))
	CV=sd(d)/mean(d)
	
	return(data.frame(Rao=Q,FDis=FDis,CV=CV))
}
