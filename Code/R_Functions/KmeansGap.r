## Reference for the gap statistic method: Tibshirani, Walther, and Hastie. J. R. Statist. Soc. B (2001) 63, Part 2, pp. 411-423

## KmeansGap(): Given a set of species traits and abundances, determines whether species are statistically clustered by traits,
##				finds the number of clsuters. Needs packages 'plyr' and 'pracma'.
## ARGUMENTS
## dat 				: data.frame	Should have an abundance column named N and either a trait column named trait 
##										or multiple trait columns for higher-dimensional tests.
## nozeros 			: log 			If TRUE, discards species with zero abundance from the null communities.
## multiD 			: log 			If FALSE, runs analysis on a single trait axis, uses the column from dat named "trait" (returns error if no such column exists).
##										If TRUE, runs analysis using all columns of dat other than N (multitrait analysis).
## nullmodel 		: chr 			If 'shuffle', observed abundances are permutated across observed traits. 
##                    					If 'draw', null traits are drawn de novo from U(0,1), and observed abundances permutated across them. 
## numnulls 		: int 			Number of null communities to test against.
## mink 			: int 			Minimum number of clusters to search for.
## maxk 			: int 			Maximum number of clusters to search for.
## nstartingpoints 	: int			Fed to argument "centers" of function kmeans(). Number of different random starting points for the clusters.
## weighting 		: chr 			Type of weight applied to within-cluster dispersion. 
##										If weighting = 'tibshirani', uses dispersion as defined in Tibshirani et al 2001, J R Stat Soc B 63, part 2, pp 411-423.
##										If weighting = 'yan', uses correction proposed in Yan & Ye 2007, Biometrics 63, 1031-1037.
## shortcut 		: log 			If TRUE, initial centroids are evenly distributed. If FALSE, initial centroids are randomly distributed, 
##										and metric repeats analysis for a total of nstartingpoints times. 
##										TRUE runs faster, but increases the risk of missing the global optimum.
## internal.peaks 	: log			If TRUE, the metric ignores the option of a single cluster, looking instead for internal peaks in the gap curve.
##									Use when looking specifically for substructure inside a single community-wide cluster, 
##										as e.g. caused by environmental filters for intermediate traits. 
## plot 			: log 			If TRUE, plots results as a a gap curve.
## bands			: log			If TRUE, adds 95th quantile of the distribution of gap index across null communities for each number of clusters.
## plotquant90 		: log 			If TRUE, adds 90th quantile of the null communities.
## verbose 			: log 			If TRUE, prints dots on console indicating which number of clusters is currently being tested.


## VALUES
## data				: data.frame	Rows show values corresponding to each number of clusters tested. Columns are as follows
##				 						k			: num number of clusters tested (all integers between mink and maxk)
##										gap			: num gap index = difference in log dispersal between observed community and mean of null communities
##										Egap		: num mean gap index across null communities 
##				 						sdgap		: num standard deviation of the gap index across null communities 
##										nullquant95	: num 95th quantile of the gap index across null communities
##										nullquant90	: num 90th quantile of the gap index across null communities
##										logWk		: num natural logarithm of the within-cluster dispersal returned from kmeans()
##										ElogWk		: num mean of the values of logWk across the null communities
## nullmodel		: chr 			Either 'draw' or 'shuffle'
## khat				: int 			Number of clusters estimated for the observed community = number of clusters at which the gap index was maximal
## maxgap			: num 			Gap statistic = maximum value of the gap index across all number of clusters tested in the observed community
## maxnullquant95	: num 			95th quantile of the gap statistics of the null communities
## maxnullquant90	: num			90th quantile of the gap statistics of the null communities
## mink				: num 			Minimum number of clusters tested
## maxk				: num 			Maximum number of clusters tested
## z.score			: num 			Defined as (maxgap-mean(maxnullgap))/sd(maxnullgap)
## p.value			: num 			Fraction of null communities with a higher gap statistic than the observed community.
## weighting		: chr 			Prints the input value of weighting. See ARGUMENTS.
## dispersal		: chr 			Prints the corresponding dispersion weighting. Let D be the within-cluster sum of squared distances:
##										D/2n is used in Tibshirani et al. 2001. 
##										D/(2n(n-1)) is used in Yan & Ye 2007. 
KmeansGap=function(
					dat,
					nozeros=FALSE,
					multiD=FALSE,
					nullmodel='shuffle',
					numnulls=100,
					mink=1,
					maxk=NULL,
					nstartingpoints=100,
					weighting='tibshirani',
					shortcut=FALSE,
					internal.peaks=FALSE,
					plot=FALSE,
					bands=TRUE,
					plotquant90=TRUE,
					verbose=TRUE
				){
	library(pracma)
	set.seed(0)
	
	if(class(dat)!='data.frame') stop('Input to KmeansGap must be a data frame')
	if(!multiD & !'trait'%in%names(dat)) stop('Input data frame for the 1d case must have a "trait" column')
	if(!'N'%in%names(dat)) stop('Input data frame must have a "N" column')
	
	## prepare data for analysis
	dat=dat[complete.cases(dat),] 						# remove NAs	
	if(nozeros) dat=dat[dat$N>0,]						# remove species with zero abundance
	dat$N=dat$N/min(1,min(dat$N[dat$N>0]))				# normalize abundances such that the rarest species with positive abundance has abundance 1
	if(!multiD) dat=dat[order(dat$trait),]				# order species by trait
	
	rescale=function(x) (x-min(x))/(max(x)-min(x))
	traitdata=if(multiD) as.data.frame(apply(subset(dat,select=-N),2,rescale)) else data.frame(trait=rescale(dat$trait))
	
	## If maxk not provided, set to minimum of 20 or half the number of species
	if(is.null(maxk)) maxk=round(min(20,sum(dat$N>0)/2))
	
	## vector of number of clusters to test
	kvec=mink:maxk	
	
	## take trait data frame and abundance vector, create a data frame/vector where each row/element is an individual if multiD=TRUE/FALSE
	rep.row=function(r,a) colwise(function(x) rep(x,a))(r)
	Com=function(traits,abuns) if(ncol(traits)>1) rep.row(r=traits,a=abuns) else rep(traits$trait,abuns)	
	
	## kms() applies kmeans() to either a vector where each element is the trait of one individual in the community (multiD = FALSE), 
	## or a data frame where each row lists the traits of one individual in the community (multiD = TRUE).
	## If shortcut = TRUE, will first try kmeans algorithm with evenly spaced initial centroids; 
	## If algorithm fails to converge or shortcut = FALSE, will assign random initial locations and repeat the process nstartingpoints times.
	kms=function(com,k){
		if(k==1){ 
			z=kmeans(com,centers=1,iter.max=100)
		}else{
			nullcenters=if(!multiD) seq(min(com),max(com),l=k) else apply(unique(com),2,function(v) sample(seq(min(v),max(v),l=k)))
			if(shortcut){
				z=try(kmeans(com,centers=nullcenters,iter.max=100),silent=TRUE)
				if(class(z)=='try-error') z=kmeans(com,centers=k,nstart=nstartingpoints,iter.max=100)
			}else z=kmeans(com,centers=k,nstart=nstartingpoints,iter.max=100)
		}
		return(z)
	}
	
	## read within-cluster dispersion (sum of squared trait distances between individuals and centroids) and apply desired weighting.
	## Note: kmeans()$withinss is equivalent to Tibshirani's D_r/2n_r (the summand in Tibshirani's Eq 2, Yan & Ye's Eq 1), where
	## D_r is the summed pairwise squared distances within a cluster (summed among all pairs, with ij counting and ji counting as well).
	## The scaling by 2n_r, where n_r is the number of elements in cluster r, makes it identical to the summed squared distances to the centroid, 
	## i.e. k-means' withinss.
	Dispersion=function(dtf,weight) with(dtf,{
		if(weight=='tibshirani') return(tot.withinss)						## Tibshirani's W_k (Eq 2 in Tibshirani et al 2001, Yan & Ye 2007)
		if(weight=='yan') return(sum(withinss/ifelse(size>1,size-1,1)))		## Correction from Yan & Ye 2007 (Eq 5 in Yan & Ye 2007)
	})
	
	## perform the kmeans test on the observed community
	Wk=sapply(kvec,function(k){ 	## Wk = within-cluster sum of squares (ie squared distances); D in Tibshirani et al 2001
		com=Com(traitdata,dat$N)
		mod=kms(com,k)
		return(Dispersion(dtf=mod,weight=weighting))
	})
	
	## perform the kmeans test on the null communities.
	## If nullmodel = 'shuffle', will permutate observed species abundances across observed species traits. 
	## If nullmodel = 'draw', will generate de novo traits, and permutate observed species abundances.
	nullWk=sapply(kvec,function(k){
		if(verbose) cat(".",if(k%%10==0 | k==maxk) paste(k,"\n"))
		
		ndat=sapply(seq(numnulls),function(run) sample(dat$N))
		return(apply(ndat,2,function(N){ 
			if(nullmodel=='shuffle') nulltraitdata=traitdata
			if(nullmodel=='draw'){
				if(!multiD) nulltraitdata=data.frame(trait=scale(sort(runif(nrow(dat)))))
				if(multiD) nulltraitdata=apply(nulltrait,2,function(v) sort(runif(length(v),min=min(v),max=max(v))))
			}
			com=Com(nulltraitdata,N)
			mod=kms(com,k)
			return(Dispersion(dtf=mod,weight=weighting))
		}))
		
	})
	
	nullN=sapply(seq(numnulls),function(run) sample(dat$N))
	if(nullmodel=='shuffle') nulltrait=traitdata
	
	logWk=log(Wk)										## vector length(kvec) 
	lognullWk=log(nullWk)								## matrix numnulls x length(kvec)
	ElogWk=apply(lognullWk,2,mean,na.rm=TRUE)			## vector length(kvec)
	Gapk=ElogWk-logWk									## vector length(kvec)
	nullGapk=t(ElogWk-t(lognullWk))						## matrix numnulls x length(kvec)
	
	if(internal.peaks){
		peaks=pracma::findpeaks(Gapk)
		maxgap=ifelse(!is.null(peaks),max(peaks[,1]),NA)					## scalar
		khat=ifelse(!is.null(peaks),peaks[which.max(peaks[,1]),2],NA)		## scalar	
		maxnullgap=apply(nullGapk,1,function(.){							## vector numnulls
			peaks=pracma::findpeaks(.)
			maxgap=ifelse(!is.null(peaks),max(peaks[,1]),0)
		})
	}else{
		maxgap=max(Gapk)													## scalar
		khat=kvec[which.max(Gapk)]											## scalar
		maxnullgap=apply(nullGapk,1,max,na.rm=TRUE)							## vector numnulls
	}
	
	z.score=(maxgap-mean(maxnullgap))/sd(maxnullgap)	## scalar
	p.value=sum(maxnullgap>maxgap)/length(maxnullgap)	## scalar
	
	mngap=apply(nullGapk,2,mean,na.rm=TRUE); 			## mean  of the null gaps --- vector length(kvec)
	sdgap=apply(nullGapk,2,sd,na.rm=TRUE); 				## stdev of the null gaps --- vector length(kvec)
	nullquant95=apply(nullGapk,2,function(vec) as.numeric(quantile(vec,.95,na.rm=TRUE)))	## 95% quantile of the null gaps --- vector length(kvec)
	nullquant90=apply(nullGapk,2,function(vec) as.numeric(quantile(vec,.9,na.rm=TRUE)))	## 90% quantile of the null gaps --- vector length(kvec)
	maxnullquant95=quantile(maxnullgap,.95)				## 95% quantile of the max null gaps --- scalar
	maxnullquant90=quantile(maxnullgap,.9)				## 90% quantile of the max null gaps --- scalar
	
	if(plot){		
		plot(0,t='n',xlim=c(min(kvec),max(kvec)),ylim=c(0,max(Gapk,maxnullquant95)),las=1,xlab='No. clusters',ylab='Gap')
		if(bands) polygon(x=c(kvec,rev(kvec)),y=c(mngap,rev(nullquant95)),col='grey90',border=NA)
		lines(kvec,rep(maxnullquant95,length(kvec)),col=2)
		if(plotquant90){ 
			if(bands) polygon(x=c(kvec,rev(kvec)),y=c(mngap,rev(nullquant90)),col='grey80',border=NA)
			lines(kvec,rep(maxnullquant90,length(kvec)),col=2,lty=2)
		}
		points(kvec,Gapk,t='o',pch=20,lwd=2)
		polygon(x=c(kvec,rev(kvec)),y=c(mngap,rep(-maxnullquant95,length(kvec))),col='white',border=NA)
		box()
	}
	
	if(weighting=='tibshirani') disp='D/2n'; if(weighting=='yan') disp='D/(2n(n-1))'
	
	return(list(
		data=data.frame(
			k=kvec,
			gap=Gapk,
			Egap=mngap,
			sdgap=sdgap,
			nullquant95=nullquant95,
			nullquant90=nullquant90,
			logWk=logWk,
			ElogWk=ElogWk
		),
		nullmodel=nullmodel,
		khat=khat,
		maxgap=maxgap,
		maxnullquant95=maxnullquant95,
		maxnullquant90=maxnullquant90,
		mink=mink,
		maxk=maxk,
		z.score=z.score,
		p.value=p.value,
		weighting=weighting,
		dispersion=disp
	))
}


## Calculates gap statistic from already-run kmeans metric. 
## Receives as input the data file and the null file containing all null statistics.
GapC=function(data,nulls,weighting,internal.peaks=FALSE){
	library(plyr)
	if(weighting=='tibshirani'){ 
		data$logWk=data$log_Wk_tibshirani
		nulls$logWk=nulls$log_Wk_tibshirani
	}
	if(weighting=='yanye'){ 
		data$logWk=data$log_Wk_yanye
		nulls$logWk=nulls$log_Wk_yanye
	}
	ElogWk=ddply(nulls,.(k),function(.) c(ElogWk=mean(.$logWk)))			## data frame with columns k, ElogWk
	nulls=merge(nulls,ElogWk,by='k') %>% mutate(Gapk=ElogWk-logWk)
	data=merge(data,ElogWk,by='k') %>% mutate(Gapk=ElogWk-log_Wk_tibshirani)						## data frame with columns k, Gapk
	
	if(internal.peaks){
		peaks=pracma::findpeaks(data$Gapk)
		maxgap=ifelse(!is.null(peaks),max(peaks[,1]),NA)					## scalar
		khat=ifelse(!is.null(peaks),data$k[which.max(peaks[,1]),2],NA)		## scalar	
		maxnullgap=ddply(nulls,.(run),function(.){ 
			peaks=pracma::findpeaks(.$Gapk,na.rm=TRUE)
			return(c(maxnullgap=ifelse(!is.null(peaks),max(peaks[,1]),0)))
		})$maxnullgap	## vector numnulls
		
	}else{
		maxgap=max(data$Gapk)
		khat=with(data,k[which.max(Gapk)])
		maxnullgap=ddply(nulls,.(run),function(.) c(maxnullgap=max(.$Gapk,na.rm=TRUE)))$maxnullgap	## vector numnulls
	}
	
	z.score=(maxgap-mean(maxnullgap))/sd(maxnullgap)	## scalar
	p.value=sum(maxnullgap>maxgap)/length(maxnullgap)	## scalar

	nullstats=ddply(nulls,.(k),function(.) with(.,c(
		Egap=mean(Gapk,na.rm=TRUE),						## mean  of the null gaps --- vector length(kvec)
		sdgap=sd(Gapk,na.rm=TRUE),						## stdev of the null gaps --- vector length(kvec)
		nullquant95=quantile(Gapk,.95,na.rm=TRUE),		## 95% quantile of the null gaps --- vector length(kvec)
		nullquant90=quantile(Gapk,.90,na.rm=TRUE)		## 90% quantile of the null gaps --- vector length(kvec)
	)))
	maxnullquant95=as.numeric(quantile(maxnullgap,.95,na.rm=TRUE))	## 95% quantile of the max null gaps --- scalar
	maxnullquant90=as.numeric(quantile(maxnullgap,.90,na.rm=TRUE))	## 90% quantile of the max null gaps --- scalar
	
	return(list(
		data=cbind(
			data.frame(
				k=data$k,
				gap=data$Gapk,
				logWk=data$logWk,
				ElogWk=data$ElogWk),
			nullstats),
		nullmodel='permutation',
		khat=khat,
		maxgap=maxgap,
		maxnullquant95=maxnullquant95,
		maxnullquant90=maxnullquant90,
		mink=1,
		maxk=20,
		z.score=z.score,
		p.value=p.value,
		dispersion=weighting,
		numnulls=length(unique(nulls$run)),
		mink=min(data$k),
		maxk=max(data$k)
	))
}

## Plots the gap curve, i.e. gap index vs number of clusters.
## Adds red line indicating the 95th quantile of the gap statistic across null communities.
## If line90 = TRUE, adds red line indicating the 90th quantile of the gap statistic across null communities.
## If bands = TRUE, plots the 95th and/or 90th quantile of the null distribution of gap index for each number of clusters.
## If whiten = TRUE, trims the plot at gap >= 0 on the y axis.
PlotKmeansGap=function(gap,ymin=NULL,ymax=NULL,linewidth=1,ylabb='Gap',xlabb='Number of clusters',line90=FALSE,whiten=FALSE,bands=FALSE,...){
	kvec=gap$data$k; Gapk=gap$data$gap; mngap=gap$data$Egap
	nullquant95=gap$data$nullquant95; nullquant90=gap$data$nullquant90
	maxnullquant95=gap$maxnullquant95; maxnullquant90=gap$maxnullquant90
	if(is.null(ymin)) ymin=0
	if(is.null(ymax)) ymax=max(Gapk,maxnullquant95)
	polymin=rep(ymin,length(kvec))
	plot(0,t='n',xlim=c(min(kvec),max(kvec)),ylim=c(ymin,ymax),las=1,xlab=xlabb,ylab=ylabb,...)
	if(bands){
		polygon(x=c(kvec,rev(kvec)),y=c(polymin,rev(nullquant95)),col='grey90',border=NA)
		polygon(x=c(kvec,rev(kvec)),y=c(polymin,rev(nullquant90)),col='grey80',border=NA)
	}
	lines(kvec,rep(maxnullquant95,length(kvec)),col=2)
	if(line90) lines(kvec,rep(maxnullquant90,length(kvec)),col=2,lty=2)
	points(kvec,Gapk,t='o',pch=20,lwd=linewidth)
	if(whiten) polygon(x=c(kvec,rev(kvec)),y=c(polymin,rep(-maxnullquant95,length(kvec))),col='white',border=NA)
	box()
	return(NULL)
}


## Plots community as abundance-by-trait stem plot, with species trait represented on the x-axis and abundance on the y-axis.
## Colors stems by cluster membership, as determined by kmeans algorithm, for the specified number of clusters K.
PlotClusters=function(trait,N,K,plot=TRUE,colors=TRUE,nstartingpoints=10000,xlab=NULL,...){
	N=N/min(1,min(N[N>0])) 
	k=kmeans(rep(trait,N),centers=K,nstart=nstartingpoints)
	dtf=unique(data.frame(trait=rep(trait,N),N=rep(N,N),cluster=k$cluster))
	dtf=dtf[order(dtf$trait),]
	dtf$cluster=match(dtf$cluster,unique(dtf$cluster))
	cols=if(colors) dtf$cluster else ifelse(dtf$cluster%%2==1,1,2)
	if(plot) with(dtf,plot(trait,N,t='h',col=cols,las=1,xlab=ifelse(is.null(xlab),'Trait',xlab),ylab='Abundance',...))
	return(dtf)	
}