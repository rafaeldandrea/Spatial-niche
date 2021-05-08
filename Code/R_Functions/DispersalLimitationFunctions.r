## "Melts" the borders of the island using heat diffusion. m0 is the landscape matrix specifying stress level of
## each site. order is the number of iterations of the diffusion procedure. k is the thermal conductivity, ie how 
## intense each melting step is. pbc = periodic boundary conditions.
Diffuse=function(m0,order,k,pbc=TRUE,plot=FALSE){
	m=m0
	L=nrow(m)
	if(order>0) for(tmp in seq(order)){
		f=matrix(0,L,L)
		for(i in 2:(L-1)) for(j in 2:(L-1)) f[i,j]=m[i+1,j]+m[i-1,j]+m[i,j+1]+m[i,j-1]-4*m[i,j]
		if(pbc){
			for(i in 2:(L-1)){ 
				f[i,1]=m[i+1,1]+m[i-1,1]+m[i,2]+m[i,L]-4*m[i,1]
				f[i,L]=m[i+1,L]+m[i-1,L]+m[i,1]+m[i,L-1]-4*m[i,L]
			}		
			for(j in 2:(L-1)){
				f[1,j]=m[2,j]+m[L,j]+m[1,j+1]+m[1,j-1]-4*m[1,j]
				f[L,j]=m[1,j]+m[L-1,j]+m[L,j+1]+m[L,j-1]-4*m[L,j]
			}
			f[1,1]=m[2,1]+m[L,1]+m[1,2]+m[1,L]-4*m[1,1]
			f[L,L]=m[1,L]+m[L-1,L]+m[L,1]+m[L,L-1]-4*m[L,L]
		}
		m=m+k*f
	}
	if(plot){
		image(seq(L),seq(L),m,las=1,xlab='Longitude',ylab='Latitude',col=terrain.colors(100))
		abline(h=seq(0,L,l=L+1)+L/(1+L)/2,v=seq(0,L,l=L+1)+L/(1+L)/2,col='grey')
	}
	return(m)
}

## Instead of melting island boundaries, creates an s-by-s patchwork of equal sized stress domains, where s^2 = S (no. of stress levels),
## and generates the stress landscape by either drawing one of S values for each pixel based on distances to domain centers (discrete=TRUE), 
## or by taking a weighted average between the S by distance to domain centers (discrete=FALSE). The latter option produces a smooth stress landscape.
## Parameter autocor controls the unbroken solidity of the patches. autocor = 0 leads to a random matrix, autocor = inf leads to a patchwork of perfect squares.
## Only works when S is a square number.
Autocorrelate=function(stlevels,autocor,discrete=TRUE,plot=FALSE){
	root=sqrt(stlevels)
	stopifnot(round(root)==root)
	tmp=merge(x=seq(1,2*root-1,by=2)*L/2/root,y=seq(1,2*root-1,by=2)*L/2/root); tmp=tmp[sample(nrow(tmp)),]
	xspots=tmp$x; yspots=tmp$y
	values=seq(stlevels)
	stress=matrix(,L,L)
	for(x in seq(L)) for(y in seq(L)){
		d=sqrt((x-xspots)^2+(y-yspots)^2)
		p=exp(-(autocor*d/(sqrt(2)*L))^2)
		if(discrete & any(p>1e-15)) stress[x,y]=sample(values,size=1,prob=p)
		if(!discrete) stress[x,y]=mean(values*p/sum(p))
		if(discrete & all(p<1e-15)){ 
			xi=unique(sort(tmp$x))[findInterval(x,seq(1,L,l=root+1),rightmost.closed=1)]
			yi=unique(sort(tmp$y))[findInterval(y,seq(1,L,l=root+1),rightmost.closed=1)]
			stress[x,y]=values[which(tmp$x==xi & tmp$y==yi)]
		}
	}
	if(plot){
		image(seq(L),seq(L),stress,las=1,xlab='Longitude',ylab='Latitude',col=terrain.colors(L^2),las=1)
		abline(h=seq(0,L,l=L+1)+L/(1+L)/2,v=seq(0,L,l=L+1)+L/(1+L)/2,col='grey')
		text(xspots,yspots,labels=values)
	}
	return(stress)
}

## Creates a stress landscape with (number) square stress patches of side (side). A proportion (noise) of the sites in the landscape
## have a random stress value, thus decreasing stress uniformity within patches.
Patches=function(side,number,noise,L=50,seed=0,plot=FALSE){
	set.seed(seed)
	l=side*sqrt(number); stopifnot(round(l)==l); stopifnot(L>=side*sqrt(number)); stopifnot(noise>=0 & noise<=1)
	tmp=merge(x=1:L,y=1:L) 
	s=with(tmp,1+(x-1)%/%side+10*(1+(y-1)%/%side))						## creates (number) square patches of side (side) with unique, uniform stress levels
	s[tmp$x>l | tmp$y>l]=NA												## apply NA to sites beyond the sublandscape of interest
	s=match(s,sort(unique(s)))											## normalizes stress levels to [1:(number)]
	s[!is.na(s)]=match(s[!is.na(s)],sample(setdiff(unique(s),NA)))		## randomizes stress levels across patches
	k=which(!is.na(s)); lk=length(k)
	s[sample(k,size=round(noise*lk))]=sample(s[k],size=round(noise*lk))	## randomizes stress within patches according to parm (noise)
	
	if(plot){
		image(seq(L),seq(L),matrix(s,L,L),las=1,xlab='Longitude',ylab='Latitude',col=terrain.colors(L^2),las=1)
		segments(x0=seq(0,l,l=l+1)+l/(1+l)/2,y0=rep(l/(1+l)/2,l+1),y1=l+rep(l/(1+l)/2,l+1),col='grey')
		segments(y0=seq(0,l,l=l+1)+l/(1+l)/2,x0=rep(l/(1+l)/2,l+1),x1=l+rep(l/(1+l)/2,l+1),col='grey')
	}
	return(s)
}

RandomField0=function(L=50,beta=10,psill=25,range=15,noise=0,varfrac=1,totvar=NULL,plot=FALSE,seed=0){
	# set seed of random number generator
	set.seed(seed)
	
	# unconditional simulations on a L x L grid using gstat
	library(gstat); library(sp)
	
	# create structure
	xy=expand.grid(1:L,1:L)
	names(xy)=c("x","y")
	
	# define the gstat object (spatial model)
	g.dummy=gstat(formula=z~1,locations=~x+y,dummy=TRUE,beta=beta,model=vgm(psill=psill,model="Gau",range=range+1e-16),nmax=20)
	
	# make simulation based on the stat object
	yy=predict(g.dummy,newdata=xy,nsim=1); gridded(yy)=~x+y
	tmp=yy@data[,1]; yy@data[,1]=(tmp-min(tmp))/(max(tmp)-min(tmp))
	stress=yy@data[,1]
	
	# add noise
	if(noise>0){ 
		ind=sample(L^2,size=noise*L^2)
		stress[ind]=stress[sample(ind)]
	}
	
	# scale the variance -- affects regional heterogeneity
	stress=stress*sqrt(varfrac)
	if(!is.null(totvar)) stress=stress*sqrt(totvar/var(stress))
	
	# plot
	yy@data[,1]=stress
	if(plot) print(spplot(yy,col.regions=terrain.colors(L^2)))
	
	return(stress)
}

RandomField=function(L=50,rangepar,sillpar,nuggetpar,seed,plot=FALSE){
	stopifnot(nuggetpar>=0 & nuggetpar<=1)
	library(RandomFields)
	RFoptions(seed=seed)
	stress=RFsimulate(RMgauss(scale=rangepar+1e-16,var=2*sillpar*(1-nuggetpar))+RMtrend(mean=0)+RMnugget(var=2*sillpar*nuggetpar),x=1:L,y=1:L)@data$variable1
	# if(plot) image(1:L,1:L,matrix(stress,L,L),col=terrain.colors(L^2),las=1)
	if(plot) plot(raster::raster(matrix(stress,L,L)),las=1,xlab='x-coordinate',ylab='y-coordinate')
	return(stress)
}

GaussVgm=function(range,h,sill,nugget) (sill-nugget)*(1-exp(-h^2/range^2))+nugget*(h>0)		## when h = range, the variogram is at 0.632 of the sill
sumGauss=function(range,h,gamma,sill,nugget) sum((gamma-GaussVgm(range,h,sill,nugget))^2)
FitRange=function(h,gamma,sill,nugget,L) optimize(f=sumGauss,interval=c(0,L),h,gamma,sill,nugget)$minimum

Variogram=function(stress,plot=FALSE,coloroption=0,cutoff){
	if(is.matrix(stress)) L=nrow(stress) else if(is.numeric(stress)) L=sqrt(length(stress)) else stop('Invalid input')
	
	library(sp)
	library(gstat)
	
	stress=as.numeric(stress)
	
	data=cbind(merge(data.frame(x=seq(L)-1),data.frame(y=seq(L)-1)),stress=stress)
	coordinates(data)=~x+y
	sample.vgm=variogram(stress~1,data=data,width=1,cutoff=cutoff)	

	fitted.vgm=vgm(	model="Gau",
					psill=var(stress),				
					nugget=sample.vgm$gamma[1],		
					range=FitRange(h=sample.vgm$dist,gamma=sample.vgm$gamma,sill=var(stress),nugget=sample.vgm$gamma[1],L=L)
	)
	
	sill=fitted.vgm[2,2]		## semivariance at the landscape level --> regional heterogeneity
	nugget=fitted.vgm[1,2]		## semivariance at distance = 1	--> local uniformity
	range=fitted.vgm[2,3]		## distance at which the semivariance reaches 63% of the sill

	res=data.frame(dist=c(0,sample.vgm$dist),gamma=c(0,sample.vgm$gamma),fitted=GaussVgm(range,h=c(0,sample.vgm$dist),sill,nugget))
	
	range=sqrt(3)*range		## redefining the range as the distance at which the semivariance reaches 95% of the sill
	
	if(plot){
		par(mfrow=c(1,2),mar=c(4.5,4.5,2,2),las=1,pch=20,cex.lab=1.5)
		if(coloroption==1){
			y=sort(findInterval(stress,vec=seq(0,1,l=10*L^2)))
			image(1:L,1:L,matrix(stress,L,L),xlab='Longitude',ylab='Latitude',col=terrain.colors(10*L^2)[y])
		} else image(1:L,1:L,matrix(stress,L,L),xlab='Longitude',ylab='Latitude',col=terrain.colors(10*L^2))
		
		with(res,{
			plot(dist,gamma,t='p',xlim=c(0,cutoff),ylim=c(0,1.2*sill),ylab=expression(gamma),xlab='Distance')
			lines(dist,fitted,col=blue,lwd=3)
			abline(h=c(sill,nugget),col=c(green,yellow)); abline(v=range,col=red)
			text(x=range-1.5,y=sill/2,labels=paste('range = ',round(range,1)),col=red,srt=90)
			text(x=5,y=sill*1.02,labels=paste('sill = ',round(sill,2)),col=green)
			text(x=cutoff-10,y=nugget+.02*sill,labels=paste('nugget = ',round(nugget,2)),col=yellow)
			box()
		})		
	}
	return(list(data=res,nugget=nugget,sill=sill,range=range))	
}

RandomFieldBin=function(x=NULL,levels,...){ 
	if(is.null(x)) stress=as.numeric(cut(RandomField(...),breaks=seq(min(x)-1e-16,max(x)+1e-16,l=levels)))
	if(!is.null(x)) stress=as.numeric(cut(x,breaks=seq(min(x)-1e-16,max(x)+1e-16,l=levels)))
	return(stress)
}

BrokenStick=function(n,seed){
	set.seed(seed)
	p=numeric(n)
	p[1]=runif(1)
	for(i in 2:(n-1)) p[i]=runif(1,min=p[i-1],max=1)
	p[n]=1
	p=c(p[1],diff(p))
	return(sample(p))
}