## Functions for calculating varD(t), where D = (n(t)-n0)/sqrt(n0),
## using the formula for P(n,n0,t) from Chisholm & O'Dwyer 2014.

Prob_n_n0_t=function(t,n,n0,nu){
  omega=exp(-nu*t)
  Omega=(1/omega-1)*(1-nu)/nu
  
  lprob=-n*log(1/omega/Omega+1) +
    -n0*log(Omega*(1+omega*Omega)) + 
    lchoose(n+n0-1,n) +
    log(Re(hypergeo::hypergeo(1-n0,-n0,1-n0-n,(1-Omega)*(1+omega*Omega))))
  
  prob=exp(lprob)
}

E_n0_t=function(t,n0,nu,maxn=1e3){
  n=1:maxn
  sum(n*Prob_n_n0_t(t,n,n0,nu))
}

Var_n0_t=function(t,n0,nu,maxn=1e3){
  n=1:maxn
  sum(n^2*Prob_n_n0_t(t,n,n0,nu))-E_n0_t(t,n0,nu,maxn)^2
}

VarD_n0_t=function(t,n0,nu,maxn=1e3){
  Var_n0_t(t,n0,nu,maxn)/n0
}
