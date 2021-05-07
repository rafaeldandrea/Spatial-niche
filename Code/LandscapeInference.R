## generate landscapes for inference

## load libraries
library(tidyr)
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps

## Load functions
loadfuns <- 
  lapply(
    list.files(
      path='~/R_Functions/',
      pattern = "[.][R]$",
      full.names=TRUE,
      ignore.case=TRUE,
      recursive=TRUE
      ),
    source
  )


scenarios <-
  tibble(
      S=5, ## number of species
      NST=5, ## number of soil types
      A=1e4, ## number of sites
      NC=25, ## number of cells
      non_preferred_site_affinity=0.2, ## affinity for non-preferred sites (between 0 and 1)
      seed=0 ## seed for the random number generator
  )

## Assign each column of scen to a variable matching its name
list2env(as.list(scenarios),envir=globalenv())

## affinity matrix
foobar=c(1,rep(non_preferred_site_affinity,NST-1))
C=matrix(non_preferred_site_affinity,S,NST); for(i in 1:S) C[i,]=guyrot(v=foobar,n=i-1)

## landscape matrix [sqrt(A) x sqrt(A)]
foo=Patches(side=sqrt(A/NC),number=NC,noise=0,L=sqrt(A))
cell=match(foo,unique(foo))
set.seed(seed) ## change the parameter 'seed' to obtain different landscapes
Landscape=sample(NST,replace=TRUE,size=NC)[match(foo,unique(foo))]

## community matrix [sqrt(A) x sqrt(A)]
Community=sapply(seq(A),function(site) sample(S,size=1,prob=C[seq(S),Landscape[site]]))

## Tabulation matrix [S x NC]
Tabulation=tibble(cell=cell,soiltype=Landscape,species=Community) %>%
  group_by(cell) %>%
  count(species) %>%
  spread(key=species,value=n,sep='_')

## show landscape matrix
image(matrix(Landscape,sqrt(A),sqrt(A)),main='Landscape',las=1)

## show community matrix
image(matrix(Community,sqrt(A),sqrt(A)),main='Community matrix',las=1)
