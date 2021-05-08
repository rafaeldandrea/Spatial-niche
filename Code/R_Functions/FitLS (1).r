## Predicte number of species for a given abundance (r)
pred.logser=function(r,alfa,N){
	X=N/(alfa+ N)
	alfa*X^r/r
}

## Logseries PDF
## Suppl Mat. of Alonso et al. 2008 Ecol.Lett. 11:93-105
dlr <- function(n,N,alfa){
	x<-N/(N+alfa)
	gama<- function(x) {1/log(1/(1-x))}
	gama(x)*(x^n)/n
}

cdf.logser=function(k,X) -1/log(1-X)*zipfR::Ibeta(X,k,1e-10)


## Fit logseries given an abundance vector or
## just the number of individuals and species
## if abundance vector is given, returns also the log-likelihood and AIC 
fit.logser=function(x=NULL,size=NULL,rich=NULL){
  if(!is.null(x)){
    S <- length(x)
    N <- sum(x)
    if(!is.null(size)|!is.null(rich)){
      warning(paste("Model fitted with size = ",N," and rich = ",S," \n calculated from supplied abundances"))
    }
  }
  if(is.null(x)&!is.null(size)&!is.null(rich)){
    S <- rich
    N <- size
  }

  if(is.null(x)&is.null(size)&is.null(rich)){
    stop("Please provide size and species number or a vector of species abundances")
  }
  f1=function(X) {S/N - ((1-X)/X)*(-log(1-X))}
  sol=uniroot(f1,interval=c(0.1,0.999999))
  X=sol$root
  alfa=N*(1-X)/X
  results <- list(summary=c(N=as.integer(N),S=as.integer(S),alfa=alfa,X=X,loglik=NULL,AIC=NULL),sol=sol)
  if(!is.null(x)){
    LL <- sum(log(dlr(x,N,alfa)))
    aic <- -2*LL+2
    results$summary <- c(results$summary,loglik=LL,AIC=aic)
    results
  }
  results
}

## Creates a spline function for the logseries
## given an abundance vector
## if Fisher alpha is not given, alpha is estimated by maximum likelihood 
# rad.logserf <- function(x,alfa=NULL){
  # if(is.null(alfa)){
    # a <- fit.logser(x)$summary["alfa"]
  # }
# else a <- alfa
  # N <- sum(x)
  # splinefun(x=cumsum(pred.logser(N:1,a,N)),y=N:1,"natural")
# }

## Testing the functions
## Package untb to generate log-serie communities
# library(untb)
## Package vegan for rank-ab plots & model fitting
## Generating community 10000 individuals, 100 species
# ex1 <- fisher.ecosystem(10000,100,10000)
## Rank-abund plot
# plot(radfit(ex1))
## Adding the logseries curve
# curve(rad.logserf(ex1)(x),add=T, col="darkorange",lwd=2, lty=2)