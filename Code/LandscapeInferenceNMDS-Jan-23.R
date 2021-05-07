
## libraries ==========
library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps
library(fpc) ## for clustering analysis
library(numbers) ## for function divisors()
library(MASS) ## for generalized inverse (left inverse) of Cmatrix when S > NST
library(dbscan) ## for HDBSCAN clustering method
library(pracma) ## for function Mode()
library(plyr)

count <- dplyr::count
select <- dplyr::select
summarize <-dplyr::summarize
red <- "#DD5144"

## Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc <- !(user <- Sys.info()["user"]) %in% c("rafael", "wangdz")


## Read scenario from either input on command line (if running on hpcc)
## or chosen scenario listed in the line below (if running on PC)
simind <- ifelse(hpcc, as.numeric(commandArgs(TRUE)[1]), 1601)

## Load functions =========
loadfuns <- {
  lapply(
    list.files(
      path = switch(
        match(user, c("rafael", "wangdz", "rdandrea", "dw9")),
        "~/R_Functions",
        "/Users/wangdz/Documents/UIUC/Odywer/R_Functions",
        "/data/rdandrea-dw9/R_Functions",
        "/data/rdandrea-dw9/R_Functions"
      ),
      pattern = "[.][R]$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ),
    source
  )
}

## operation directives ==========
landscape.patches <- FALSE
landscape.randomfield <- TRUE

clustering.method <- "pam"

do.kmeansgap <- (clustering.method == "kmeans")
do.pamk <- (clustering.method == "pam")
do.hierarch <- (clustering.method == "hierarch")
do.NMDS <- (clustering.method == "nmds")
do.hdbscan <- (clustering.method == "hdbscan")

plot.landscape_community_C_matrices <- FALSE
plot.landscape_community_variogram <- FALSE
plot.clustering.statistics.against.cell.size <- FALSE

analyze.data <- FALSE
analyze.data2 <- FALSE

calculate.likelihood <- FALSE
simulate.likelihood <- FALSE
simulate.likelihood2 <- TRUE



## scenarios ===========
scenarios <- {
  crossing(
    tibble(
      S = c(10, 100), ## number of species
      NST = c(5, 10) ## number of soil types
    ),
    tibble(
      A = 120^2, ## number of sites
      NC = A / numbers::divisors(sqrt(A))^2, ## number of cells
      non_preferred_site_affinity = 0.1, ## affinity for non-preferred sites (between 0 and 1)
    ),
    tibble(
      C_sd = c(0, .01, .05, .1)
    ),
    tibble(
      landscape_autocor_parm = c(15, 30, 100)
    )
  ) %>%
    mutate(
      mode = "spec",
      C_mean = NA,
      C_cv = NA
    ) %>%
    rbind(
      crossing(
        tibble(
          S = c(10, 100), ## number of species
          NST = c(10, 100) ## number of soil types
        ),
        tibble(
          A = 120^2, ## number of sites
          NC = A / numbers::divisors(sqrt(A))^2, ## number of cells
          non_preferred_site_affinity = NA, ## affinity for non-preferred sites (between 0 and 1)
        ),
        tibble(
          mode = "gen",
          C_mean = 1,
          C_cv = c(.2, .4, .6),
          C_sd = NA
        ),
        tibble(
          landscape_autocor_parm = c(15, 30, 100)
        )
      )
    ) %>%
    filter(NC > NST & NC < A) %>% ## remove values of NC that do not make sense given A and NSTrd
    crossing(tibble(run = 1:100)) %>% ## seed for the random number generator
    arrange(desc(mode), desc(landscape_autocor_parm), C_sd, S, NC)
}

simind_vec <- 1 + (simind - 1) * 748 + 0:747
simind_vec <- simind_vec[simind_vec <= nrow(scenarios)]

## functions
{
  
  # ReorderMatrix=function(data){
  #   order_matrix=NST+1-t(apply(data,1,rank))
  #   preferred=t(apply(order_matrix,1,order))
  #   mylist=NULL
  #   for(i in 1:nrow(data)){ 
  #     for(j in 1:ncol(data)){ 
  #       if(!preferred[i,j] %in% mylist){ 
  #         mylist=c(mylist,preferred[i,j])
  #         break
  #       }
  #     }
  #   }
  #   data=data[,mylist]
  #   return(data)
  # }
  
  ReorderMatrix=function(data){
    order_matrix=ncol(data)+1-t(apply(data,1,rank))
    preferred=t(apply(order_matrix,1,order))
    mylist=NULL
    for(i in 1:nrow(data)){ 
      for(j in 1:ncol(data)){ 
        if(!preferred[i,j] %in% mylist){ 
          mylist=c(mylist,preferred[i,j])
          break
        }
      }
    }
    if(length(mylist)<ncol(data)) mylist=c(mylist,setdiff(seq(ncol(data)),mylist))
    data=data[,mylist]
    return(data)
  }
  
  ## Landscape: vector whose A elements are the soil types of each site ========
  Generate_Landscape <- function(seed,landscape_autocor_parm) {
    if (landscape.patches) {
      patches <- Patches(side = sqrt(A / NC), number = NC, noise = 0, L = sqrt(A))
      set.seed(seed) ## change the parameter 'seed' to obtain different landscapes
      Landscape <- sample(NST, replace = TRUE, size = NC)[match(patches, unique(patches))]
    }
    if (landscape.randomfield) {
      ## create a Gaussian random field
      field <- RandomField(L = sqrt(A), rangepar = landscape_autocor_parm, sillpar = 1, nuggetpar = 0, seed = seed)
  
      ## simplify the field by boxing sites into one of NST soil types
      Landscape <- 1 + findInterval(field, quantile(field, seq(NST) / NST), rightmost.closed = TRUE)
    }
    return(Landscape)
  }
  
  ## define cell grid based on NC. (Uses pos2coord from R_Functions) ============
  Generate_Grid <- function(A, NC) {
    cellgrid <-
      as_tibble(pos2coord(pos = seq(A), dim.mat = c(sqrt(A), sqrt(A)))) %>%
      mutate(
        cx = 1 + (V1 - 1) %/% (max(V1) / sqrt(NC)),
        cy = 1 + (V2 - 1) %/% (max(V2) / sqrt(NC)),
        cellname = paste0(
          formatC(cx, width = nchar(sqrt(NC)), flag = "0"),
          formatC(cy, width = nchar(sqrt(NC)), flag = "0")
        ),
        cell = match(cellname, unique(cellname))
      ) %>%
      pull(cell)
    return(cellgrid)
  }
  
  ## affinity matrix Cmatrix ========
  Generate_Cmatrix <- function(seed,mode='spec') {
    set.seed(seed)
    if (mode == "spec") {
      foobar <- c(1, rep(non_preferred_site_affinity, NST - 1))
      Cmatrix <- matrix(non_preferred_site_affinity, S, NST)
      for (i in 1:S) Cmatrix[i, ] <- guyrot(v = foobar, n = i - 1)
      Cmatrix <- ReorderMatrix(Cmatrix)
      Cmatrix[Cmatrix != 1] <- Cmatrix[Cmatrix != 1] + rnorm(length(Cmatrix[Cmatrix != 1]), mean = 0, sd = C_sd)
      Cmatrix[Cmatrix < 0] <- 0
    }
    if (mode == "gen") {
      Cmatrix <- matrix(runif(NST * S, min = C_mean * (1 - sqrt(3) / 2 * C_cv), max = C_mean * (1 + sqrt(3) / 2 * C_cv)), S, NST)
      Cmatrix[Cmatrix < 0] <- 0
      Cmatrix <- ReorderMatrix(Cmatrix)
    }
    return(Cmatrix)
  }
  
  ## community matrix [sqrt(A) x sqrt(A)] =======
  Generate_Community <- function(Landscape, Cmatrix, seed) {
    set.seed(seed)
    R <- sapply(seq(NST), function(soiltype) sum(Landscape == soiltype))
    Community <-
      sapply(
        seq(A),
        function(site) {
          # sample(S,size=1,prob=pmax(0,rowSums(ginv(Cmatrix)))*Cmatrix[seq(S),Landscape[site]])
          sample(S, size = 1, prob = R[Landscape[site]] * Cmatrix[seq(S), Landscape[site]])
        }
      )
    return(Community)
  }
  
  ## Tabulation matrix [S x NC] ======
  Generate_Tabulation <- function(cellgrid, Community) {
    Tabulation <-
      tibble(
        cell = cellgrid,
        species = Community
      ) %>%
      group_by(cell) %>%
      count(species) %>%
      spread(key = species, value = n, sep = "_") %>%
      ungroup()
    Tabulation[is.na(Tabulation)] <- 0
    return(Tabulation)
  }
  
  ## estimate Cmatrix matrix =========
  Generate_estC <- function(est_Nij, R) {
    est_C <- t(matrix(as.numeric(unlist(est_Nij[, -1])), nrow(est_Nij), ncol(est_Nij) - 1))
    est_C <- est_C / rowSums(est_C) / R
    # est_C=round(est_C/max(est_C),1)
    #est_C <- ReorderMatrix(est_C)
    return(est_C)
  }
  
  ## find the point in a curve with the maximum |2nd derivative| ====
  Max2Deriv <- function(x, y) {
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    indices <- which.max(abs(d2))
    return(indices)
  }
  
}

# for(simind in simind_vec){
#   try({


if (simulate.likelihood2) {
  clustering.mode <- 'hierarchical'
  
  run <- 2
  A <- 120^2
  NST <- 5
  S <- 5
  landscape_autocor_parm <- 50
  C_sd <- 0.5
  C_mean <- 1
  C_cv <- 1
  Landscape <- Generate_Landscape(seed = run, landscape_autocor_parm = landscape_autocor_parm)
  R <- as.numeric(table(Landscape))
  non_preferred_site_affinity <- 0.025
  Cmatrix <- Generate_Cmatrix(seed = run, mode = 'spec')
  # Cmatrix=diag(S)
  Community <- Generate_Community(Landscape, Cmatrix, seed = run)
  
  ## knowing both C and the Landscape
  N <- as.numeric(table(Community))
  P0ij <- t(t(N*Cmatrix)/colSums(N*Cmatrix))
  foo <- tibble(i=Community,j=Landscape)
  bar <- apply(foo,1,function(v) P0ij[v[1],v[2]])
  loglik0 <- sum(log(bar))
  
  res <- 
    ddply(
      crossing(tibble(NC = unique(scenarios$NC)),tibble(NST=3:30)) %>% 
        filter(NC > S & NC <= A/9 &NC > 25), .(NST,NC), function(v) {
          NST <- v$NST
          NC <- v$NC
          print(paste(NST,NC))
          cellgrid <- Generate_Grid(A, NC)
          Tabulation <- Generate_Tabulation(cellgrid, Community)
          abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
          tN <- t(abunds_by_cells)
          tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
          
          ## estimate matrix Nij based on the identification of cells to soil types
          
          ## Using k-means
          if(clustering.mode=='kmeans'){
            clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
            cluster_index=clustermodel$cluster
            clustering.statistic=clustermodel$tot.withinss
            numclusters_cells <- NST
          }
          
          ## Using hierarchical clustering
          if(clustering.mode=='hierarchical'){
            dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
            clus_sites <- hclust(dis_sites)
            cluster_index <- cutree(clus_sites, NST)
            clustering.statistic <- NA
            numclusters_cells <- NST
          }
          
          ## Using pam clustering
          if(clustering.mode=='pam'){
            dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
            clustermodel=cluster::pam(x=dis_sites, k=NST)
            cluster_index <- clustermodel$clustering
            clustering.statistic <- as.numeric(clustermodel$objective['swap'])
            numclusters_cells <- NST
          }
            
          est_Nij <-
            Tabulation %>%
            mutate(soiltype = cluster_index) %>%
            select(-cell) %>%
            group_by(soiltype) %>%
            summarize_all(sum) %>%
            filter(soiltype > 0)
          
          est_C <- Generate_estC(est_Nij, R)
          Pij <- t(t(N*est_C)/colSums(N*est_C))
          
          
          logprob=log(Pij[, cluster_index])
          loglik <- sum(tN*logprob*is.finite(logprob),na.rm=TRUE)
          
          tr <- function(m) sum(diag(m))
          squarify <- function(m) m %*% t(m)
          center <- function(m) t(t(m) - colMeans(m))
          CC <- squarify(center(Cmatrix))
          eCeC <- squarify(center(est_C))
          RV <- tr(CC %*% eCeC) / sqrt(tr(CC %*% CC) * tr(eCeC %*% eCeC))
          
          return(
            data.frame(
              NC = NC,
              numclusters = numclusters_cells,
              loglik = loglik,
              statistic = clustering.statistic,
              RV = RV
            )
          )
        }
    ) %>%
    mutate(cellsize = sqrt(A / NC)) %>%
    as_tibble
  
  {
    df_maxloglik <-
      res %>% 
      group_by(NST) %>%
      summarize(
        bestcellsize=cellsize[which.max(loglik)],
        maxloglik=loglik[which.max(loglik)],
        maxRV=RV[which.max(loglik)]
      )
    nst_continuum <- with(df_maxloglik,seq(min(NST),max(NST),l=1e3))
    smoothed_maxloglik <- with(df_maxloglik, predict(loess(maxloglik~NST),nst_continuum))
    index <- Max2Deriv(nst_continuum, smoothed_maxloglik)
    bestNST <- unique(round(nst_continuum[index]))
    bestcellsize <- df_maxloglik %>% filter(NST==bestNST) %>% pull(bestcellsize)
    
    nst=df_maxloglik$NST
    
    NC=A/bestcellsize^2
    NST <- bestNST
    cellgrid <- Generate_Grid(A, NC)
    Tabulation <- Generate_Tabulation(cellgrid, Community)
    abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
    tN <- t(abunds_by_cells)
    tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
    
    ## Using k-means
    if(clustering.mode=='kmeans'){
      clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
      cluster_index=clustermodel$cluster
      clustering.statistic=clustermodel$tot.withinss
      numclusters_cells <- NST
    }
    
    ## Using hierarchical clustering
    if(clustering.mode=='hierarchical'){
      dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
      clus_sites <- hclust(dis_sites)
      cluster_index <- cutree(clus_sites, NST)
      clustering.statistic <- NA
      numclusters_cells <- S
    }
    
    ## Using pam clustering
    if(clustering.mode=='pam'){
      dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
      clustermodel=cluster::pam(x=dis_sites, k=NST)
      cluster_index <- clustermodel$clustering
      clustering.statistic <- as.numeric(clustermodel$objective['swap'])
      numclusters_cells <- NST
    }
    
    est_Nij <-
      Tabulation %>%
      mutate(soiltype = cluster_index) %>%
      select(-cell) %>%
      group_by(soiltype) %>%
      summarize_all(sum) %>%
      filter(soiltype > 0)
    
    est_C <- ReorderMatrix(Generate_estC(est_Nij, as.numeric(table(cluster_index))))
    
    }
  
  cellsize=unique(res$cellsize)
  plot0 <-
    res %>%
    ggplot() +
    theme_bw() +
    #scale_x_log10() +
    ylab('Log Likelihood') +
    xlab('Cell Linear Size') +
    scale_x_continuous(breaks=cellsize, labels=cellsize, minor_breaks = NULL)
  
  plot_Landscape <- 
    tibble(landscape=Landscape,community=Community) %>% 
    cbind(expand.grid(x=1:120,y=1:120)) %>%
    ggplot(aes(x,y,fill=landscape)) + 
    theme_bw() +
    geom_tile() + 
    theme(legend.position = 'none') +
    ggtitle('Landscape')  +
    scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
    scale_y_continuous(breaks=NULL,minor_breaks = NULL)
  
  plot_Community <- 
    tibble(landscape=Landscape,community=Community) %>% 
    cbind(expand.grid(x=1:120,y=1:120)) %>%
    ggplot(aes(x,y,fill=community)) + 
    theme_bw() +
    geom_tile()  + 
    theme(legend.position = 'none') +
    ggtitle('Community') +
    scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
    scale_y_continuous(breaks=NULL,minor_breaks = NULL)
  
  plot_Cmatrix <- 
    tibble(C=as.numeric(Cmatrix)) %>% 
    cbind(expand.grid(species=seq(nrow(Cmatrix)),soiltype=seq(ncol(Cmatrix)))) %>%
    ggplot(aes(x=species,y=soiltype,fill=C)) + 
    theme_bw() +
    geom_tile()  + 
    theme(legend.position = 'none') +
    ggtitle('C matrix') 
  
  plot_Cest <- 
    tibble(C=as.numeric(est_C)) %>% 
    cbind(expand.grid(species=seq(S),soiltype=seq(NST))) %>%
    ggplot(aes(x=species,y=soiltype,fill=C)) + 
    theme_bw() +
    geom_tile()  + 
    theme(legend.position = 'none') +
    ggtitle('estimated C')
  
  plot_cellsize_by_loglik <- 
    plot0 + 
    geom_line(aes(cellsize, loglik/1e3,group=as.factor(NST),color=as.factor(NST))) + 
    geom_point(aes(cellsize, loglik/1e3,group=as.factor(NST),color=as.factor(NST))) +
    ggtitle('Estimating both C and the Landscape') +
    geom_hline(yintercept=loglik0/1e3,color=red) +
    theme(legend.position='none')
  
  plot_cellsize_by_RV <- 
    plot0 + 
    geom_line(aes(cellsize, RV,group=as.factor(NST),color=as.factor(NST))) + 
    geom_point(aes(cellsize, RV,group=as.factor(NST),color=as.factor(NST))) +
    ggtitle('Match between true C and estimated C') +
    ylab('RV') +
    ylim(c(0,1)) +
    theme(legend.position='none')
  
  plot_NST_by_maxloglik <-
    df_maxloglik %>%
    ggplot(aes(NST,maxloglik/1e3)) +
    geom_vline(aes(xintercept = bestNST),color=red) +
    geom_line(aes(x, y/1e3),data=tibble(x=nst_continuum,y=smoothed_maxloglik)) +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=nst, labels=nst, minor_breaks = NULL) +
    xlab('Number of soil types') +
    ylab('Max log likelihood (/ 1e3)')
  
  plot_NST_by_maxRV <-
    df_maxloglik %>%
    ggplot(aes(NST,maxRV)) +
    geom_vline(aes(xintercept = bestNST),color=red) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=nst, labels=nst, minor_breaks = NULL) +
    xlab('Number of soil types') +
    ylab('Max RV') +
    ylim(c(0,1))
  
  plot_NST_by_bestcellsize <-
    df_maxloglik %>%
    ggplot(aes(NST,bestcellsize)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=nst, labels=nst, minor_breaks = NULL) +
    scale_y_continuous(breaks=cellsize, labels=cellsize, minor_breaks = NULL) +
    xlab('Number of soil types') +
    ylab('Best cell size') +
    ylim(range(cellsize))
  
  plotstat <- 
    plot0 + 
    geom_line(aes(cellsize, statistic)) + 
    geom_point(aes(cellsize, statistic))
  
  plotRV <- 
    plot0 + 
    geom_line(aes(cellsize, RV)) + 
    geom_point(aes(cellsize, RV)) +
    ggtitle('Match between true C and estimated C') +
    ylab('RV') +
    ylim(c(0,1))
  
  gridExtra::grid.arrange(
    plot_Landscape, 
    plot_Cmatrix, 
    plot_Community,
    plot_cellsize_by_loglik,
    plot_cellsize_by_RV,
    plot_NST_by_maxloglik,
    plot_NST_by_maxRV,
    plot_NST_by_bestcellsize,
    plot_Cest,
    nrow=3
  )
  
  # gridExtra::grid.arrange(plotL, plot1, plot2, plotC, plot3, plotRV, nrow = 2)
  # gridExtra::grid.arrange(plotL, plotCmatrix, plot1, plotC, plotCest, plotRV, nrow = 2)
}

SAD = TRUE
if(SAD==TRUE){
abundance = colSums(Tabulation)[-1]
#infer landscape
linearsize  = sqrt(size(Landscape)[2]/NC)
fooo = rep(cluster_index,each = linearsize ) # repeat every element by cell linear size
# split into size= 120, repeat every chunk by 3 and combine 
newlandscape= unlist(rep(split(fooo, ceiling(seq_along(fooo)/120)),each=linearsize),use.names = FALSE)
newCommunity= Generate_Community(newlandscape,est_C,seed = run)
newabundance = table(newCommunity)



plot_newLandscape <- 
  tibble(landscape=newlandscape) %>% 
  cbind(expand.grid(x=1:sqrt(size(newlandscape)[2]),y=1:sqrt(size(newlandscape)[2]))) %>%
  ggplot(aes(x,y,fill=landscape)) + 
  theme_bw() +
  geom_tile() + 
  theme(legend.position = 'none') +
  ggtitle('NewLandscape')  +
  scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
  scale_y_continuous(breaks=NULL,minor_breaks = NULL)

plot_estLandscape <- 
  tibble(landscape=cluster_index) %>% 
  cbind(expand.grid(x=1:sqrt(size(cluster_index)[2]),y=1:sqrt(size(cluster_index)[2]))) %>%
  ggplot(aes(x,y,fill=landscape)) + 
  theme_bw() +
  geom_tile() + 
  theme(legend.position = 'none') +
  ggtitle('New Landscape')  +
  scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
  scale_y_continuous(breaks=NULL,minor_breaks = NULL)
plot_estLandscape

plot_newCommunity <- 
  tibble(landscape=newlandscape,community=newCommunity) %>% 
  cbind(expand.grid(x=1:sqrt(size(newlandscape)[2]),y=1:sqrt(size(newlandscape)[2]))) %>%
  ggplot(aes(x,y,fill=community)) + 
  theme_bw() +
  geom_tile()  + 
  theme(legend.position = 'none') +
  ggtitle('New Community') +
  scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
  scale_y_continuous(breaks=NULL,minor_breaks = NULL)
plot_SAD <-
  plotSAD(abundance,dbn='LS')+
  theme_bw() +
  ggtitle('original SAD') 
 

plot_newSAD <-
  plotSAD(as.double(newabundance),dbn='LS')+
  theme_bw() +
  ggtitle('new SAD') 
 
gridExtra::grid.arrange(
  plot_Landscape,
  plot_newLandscape,
  plot_Community,
  plot_newCommunity,
  plot_SAD,
  plot_newSAD,
  nrow=3
)

x=as.numeric(abundance)
y=as.numeric(newabundance)
datx=plyr::count(x) %>%
  mutate(cum=(sum(freq)-cumsum(freq)+freq)/sum(freq))
daty=plyr::count(y) %>%
  mutate(cum=(sum(freq)-cumsum(freq)+freq)/sum(freq))
plottwo <-
  ggplot() +
  geom_point(aes(x,cum,color=blue),data=datx) +
  geom_line(aes(x,cum,color=blue),data=datx) +
  geom_point(aes(x,cum,color=red),data=daty) +
  geom_line(aes(x,cum,color=red),data=daty) +
  theme_bw()
plottwo
testresult = dgof::ks.test(x,y)

}

