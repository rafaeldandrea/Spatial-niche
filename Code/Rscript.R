## Code for Spatial Niches

## Libraries
library(dplyr)    ## for data wrangling
library(ggplot2)  ## for plotting
library(tidyr)    ## for function crossing()
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps

## Load functions
dummy <-
  list.files(path='~/R_Functions/',pattern = "[.][R]$",full.names=TRUE,ignore.case=TRUE,recursive=TRUE) %>%
  lapply(source)

## Set working directory
setwd('~/SpatialNiche/')

## Set/create folder to save data files in
datadir='~/SpatialNiche/Data/20191001/'
dir.create(datadir,showWarnings=FALSE)

## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc=as.numeric((nodename<-as.character(Sys.info()['nodename']))!='RAF-PC')

## set seed of random number generator
set.seed(0)

## parameter data frame
scenarios <-
  crossing(
    tibble(
      S=100,
      NST=20,
      A=5*S*NST,
      m=1,
      b=1,
      soil_type_distribution=c('homogeneous','normal')
    ),
    tibble(
      non_preferred_site_affinity=seq(0,1,l=5)
    )
  )

## Read scenario from either input on command line (if running on hpcc) 
## or chosen scenario listed in the line below (if running on PC)
ind=ifelse(hpcc,as.numeric(commandArgs(TRUE)[1]),6)

## Choose row of 'parameters'data frame to specify scenario
scen=scenarios[ind,]

## Assign each column of scen to a variable matching its name
list2env(as.list(scen),envir=globalenv())


## total number of sites per soil type
if(soil_type_distribution=='homogeneous') P=A/NST
if(soil_type_distribution=='normal'){
  P=0
  while(length(P)!=NST){
    foo=sample(NST,size=A,replace=TRUE,prob=dnorm(seq(NST),mean=0,sd=5))
    P=sort(plyr::count(foo)$freq,decreasing=TRUE)
  }
}

## affinity matrix
foobar=c(1,rep(non_preferred_site_affinity,NST-1))
C=matrix(non_preferred_site_affinity,S,NST); for(i in 1:S) C[i,]=guyrot(v=foobar,n=i-1)

## initial condition
if(soil_type_distribution=='homogeneous') N=matrix(A/S/NST,S,NST)
if(soil_type_distribution=='normal'){
  N=matrix(0,S,NST)
  for(j in 1:NST){
    chosen=sample(S,replace=TRUE,size=P[j])
    bar=as_tibble(table(chosen)) %>% mutate(chosen=as.numeric(chosen))
    N[bar$chosen,j]=bar$n
  }
}
R=P-colSums(N)

## dynamics
test=FALSE
iteration=biotime=save_index=0
DT=BT=test_statistic=NULL
t0=as.numeric(Sys.time()); Sys.sleep(1)
while(test==FALSE){
  if((as.numeric(Sys.time())-t0)>(3600*20)) test=TRUE
  iteration=iteration+1
  biotime=biotime+rexp(1,rate=sum(b*rowSums(N)*C%*%R)+sum(m*N))
  event=sample(c('birth','death'),size=1,prob=c(sum(b*rowSums(N)*C%*%R),sum(m*N)))
  if(event=='death'){
    position=sample(seq(S*NST),size=1,prob=as.numeric(m*N))
    species=1+(position-1)%%S
    if(rowSums(N)[species]>1){
      N[position]=N[position]-1
      site_type=1+(position-1)%/%S
      R[site_type]=R[site_type]+1
    }
  }
  if(event=='birth'){
    position=sample(seq(S*NST),size=1,prob=as.numeric(C*outer(b*rowSums(N),R)))
    N[position]=N[position]+1
    site_type=1+(position-1)%/%S
    R[site_type]=R[site_type]-1
  }
  
  ## stop the simulation when the following condition is met
  if(hpcc==0){
    DT=cbind(DT,m*rowSums(N))
    BT=cbind(BT,(rowSums(N)*C)%*%R)
  }
  if(hpcc==0 & iteration%%1000==0){
    lsfit=fitSAD(data=rowSums(N),dbn='LS',method='CVM')
    test_statistic=c(test_statistic,lsfit$statistic)
    #test_statistic=c(test_statistic,max(abs(1-rowMeans(DT)/rowMeans(BT))))
    print(c(length(test_statistic),round(test_statistic[length(test_statistic)],3)))
    DT=BT=NULL
    
    if(length(test_statistic)==100){
      mod=lm(test_statistic~seq_along(test_statistic))
      pvalue=summary(mod)$coef[2,4]
      plot(test_statistic)
      abline(mod)
      #if(pvalue>.05) test=TRUE
      test_statistic=NULL
    }
    
  }
  
  ## plot intermediate results
  if(hpcc==0 & iteration%%1e4==0){
    # plot(
    #   rowSums(N),
    #   t='h',
    #   ylim=c(0,max(rowSums(N))),
    #   xlab='Species',
    #   ylab='Abundance',
    #   las=1,
    #   main=paste(iteration,'iterations')
    # )
    
    fit=lsfit$dtf %>%
      as_tibble
    
    data=plyr::count(rowSums(N)) %>%
      mutate(
        cdf=cumsum(freq),
        cum=(sum(freq)-cdf+freq)/sum(freq)
      ) %>%
      as_tibble
    
    dtf <-
      data %>%
      inner_join(fit,by='x')
    
    myplot <- 
      dtf %>% 
      ggplot() + 
      geom_line(aes(x,cum.y),color=red) +
      geom_point(aes(x,cum.x)) + 
      scale_x_log10() + 
      scale_y_log10() +
      theme_bw() +
      theme(aspect.ratio = 1) +
      xlab('Abundance') +
      ylab('Cumulative probability')
    
    gridExtra::grid.arrange(myplot)
  }
  
  ## save current state
  if(hpcc==1 & round(as.numeric(Sys.time())-t0)%%1800==0){
    save_index=save_index+1
    if(!file.exists(paste0(datadir,'Parameters.RData'))){
      parameters=list(
        P=P,
        C=C,
        S=S,
        NST=NST,
        A=A,
        m=m,
        b=b
      )
      save(parameters,file=paste0(datadir,'Parameters.RData'))
    }
    
    thedata=list(
      N=N,
      biotime=biotime,
      iteration=iteration
    )
    save(thedata,file=paste0(datadir,formatC(ind,width=2,flag="0"),'_',formatC(save_index,width=2,flag="0"),'.RData'))
    
    Sys.sleep(1)
  }
}

