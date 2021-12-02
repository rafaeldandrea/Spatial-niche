library(tidyverse)
library(magrittr)
library(furrr)
library(parallel)
library(parallelDist) ## function parDist(X) calculates distances between rows of matrix X
library(pdist) ## for function pdist(X, Y), calculates distances between rows of matrices X and Y
library(readxl)
library(pcaMethods)
library(sparr) # for function bivariate.density() in KDE()
library(RPANDA)
library(RANN)


#Modified function to calculate SDP of an adjacency matrix

spects<-function (mat, meth = c("normal"), zero_bound = F) 
{
  skewness <- function(x, na.rm = FALSE) {
    if (is.matrix(x)) 
      apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
      if (na.rm) 
        x <- x[!is.na(x)]
      n <- length(x)
      (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x)) 
      sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
  }
  sigma = 0.1
  gKernel <- function(x) 1/(sigma * sqrt(2 * pi)) * exp(-(x^2)/2 * 
                                                          sigma^2)
  kernelG <- function(x, mean = 0, sd = 1) dnorm(x, mean = mean, 
                                                 sd = sd)
  dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096, 
                   from = min(x) - 3 * sd, to = max(x) + 3 * sd, adjust = 1, 
                   ...) {
    if (has.na <- any(is.na(x))) {
      x <- na.omit(x)
      if (length(x) == 0) 
        stop("no finite or non-missing data!")
    }
    sd <- (if (is.numeric(bw)) 
      bw[1]
      else bw(x)) * adjust
    X <- seq(from, to, len = n)
    M <- outer(X, x, kernel, sd = sd, ...)
    structure(list(x = X, y = rowMeans(M), bw = sd, call = match.call(), 
                   n = length(x), data.name = deparse(substitute(x)), 
                   has.na = has.na), class = "density")
  }
  integr <- function(x, f) {
    if (!is.numeric(x)) {
      stop("The variable of integration \"x\" is not numeric.")
    }
    if (!is.numeric(f)) {
      stop("The integrand \"f\" is not numeric.")
    }
    if (length(x) != length(f)) {
      stop("The lengths of the variable of integration and the integrand do not match.")
    }
    n = length(x)
    integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + 
                                                      f[1:(n - 1)]))
    return(integral)
  }
  if (meth == "standard") {
    e = eigen(graph.laplacian(graph.adjacency(mat, 
                                              weighted = T), normalized = F), symmetric = T, only.values = T)
    if (zero_bound == T) {
      x = subset(e$values, e$values > 0)
    }
    else {
      x = subset(e$values, e$values >= 1)
    }
    d = dens(log(x))
    dsc = d$y/(integr(d$x, d$y))
    principal_eigenvalue <- max(x)
    skewness <- skewness(x)
    peak_height <- max(dsc)
    gaps <- abs(diff(x))
    gapMat <- as.matrix(gaps)
    modalities <- c(1:length(gapMat))
    gapMatCol <- cbind(modalities, gapMat)
    eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
                                                                  2]))
    res <- list(eigenvalues = x, principal_eigenvalue = principal_eigenvalue, 
                asymmetry = skewness, peakedness = peak_height, eigengap = eigenGap[, 
                                                                                    1])
  }
  if (meth == "normal") {
    e = eigen(graph.laplacian(graph.adjacency(mat, 
                                              weighted = T), normalized = T), symmetric = T, only.values = T)
    x = subset(e$values, e$values >= 0)
    d = dens(log(x))
    dsc = d$y/(integr(d$x, d$y))
    principal_eigenvalue <- max(x)
    skewness <- skewness(x)
    peak_height <- max(dsc)
    gaps <- abs(diff(x))
    gapMat <- as.matrix(gaps)
    modalities <- c(1:length(gapMat))
    gapMatCol <- cbind(modalities, gapMat)
    eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
                                                                  2]))
    res <- list(eigenvalues = x, principal_eigenvalue = principal_eigenvalue, 
                asymmetry = skewness, peakedness = peak_height, eigengap = eigenGap[, 
                                                                                    1])
  }
  class(res) <- "spectR"
  return(res)
}

#Take the first census data from BCI
Lx=1000
Ly=500
d_cutoff=20

prefix = 'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci.full'
suffix = '.rdata?raw=true'

filter = dplyr::filter


bci = NULL
census=1
  bci %<>%
    bind_rows(
      get(load(url(paste0(prefix, census, suffix)))) %>%
        as_tibble() %>%
        drop_na(dbh) %>%
        mutate(census = census)  
    )

  baldeck_cutoff = 
    bci %>%
    group_by(sp) %>%
    summarize(
      baldeck = quantile(dbh, .56), 
      .group ='drop'
    )
  
  dat = 
    bci %>%
    inner_join(baldeck_cutoff) %>%
    filter(dbh > baldeck)
  
  
#Adjacency matrix calculation from Rafael's code  

  adjacency_matrix = 
    function(
      dat, 
      d_cutoff, 
      Lx, 
      Ly
    ){
      ## int p(r) * dr from 0 to d*
      cumulative_null_prob_threshold = pi * d_cutoff^2 / (Lx *Ly)
      
      ## abundance cutoff based on E[n] >= 1, where n is the number of pairs
      ## within a circle of radius d*
      abundance_threshold = round(sqrt(1 / cumulative_null_prob_threshold))
      
      abuns = 
        dat %>%
        count(sp) %>%
        arrange(desc(n))
      
      dat_filtered = 
        dat %>%
        inner_join(abuns, by = 'sp') %>%
        filter(n >= abundance_threshold, !is.na(gx), !is.na(gy))
      
      selected_species = unique(dat_filtered$sp)
      
      dat_selected = 
        dat_filtered %>%
        filter(sp %in% selected_species)
      
      sp_pairs = 
        expand_grid(
          sp1 = unique(dat_selected$sp),
          sp2 = unique(dat_selected$sp)
        ) %>%
        filter(sp1 < sp2)
      
      adjacency =
        sp_pairs %>%
        future_pmap_dfr(
          .options = furrr_options(seed = NULL),
          .f = function(sp1, sp2){
            df1 = 
              dat_selected %>% 
              filter(sp == sp1)
            
            n1 = nrow(df1)
            
            df2 = 
              dat_selected %>% 
              filter(sp == sp2)
            
            n2 = nrow(df2)
            
            distances = 
              as.vector(nn2(with(df1,cbind(gx,gy)),with(df1,cbind(gx,gy)),k=20,treetype="kd",searchtype="radius",radius=d_cutoff)$nn.idx)
            
            pairs = sum(distances > 0)
            
            tibble(
              sp1 = sp1, 
              sp2 = sp2, 
              weighted = (pairs / (n1 * n2)) / cumulative_null_prob_threshold,
              binary = 1 * (weighted > 1)
            ) 
          }
        )
      
      adjacency %<>% 
        bind_rows(
          adjacency %>% 
            rename(sp1 = sp2, sp2 = sp1)
        )
      
      A<-adjacency
      A = A %>% 
        complete(
          sp1, 
          sp2, 
          fill = list(weighted = 0, binary = 0)
        ) %>%
        select(-binary) %>% 
        pivot_wider(
          names_from = sp2, 
          values_from = weighted
        ) %>% 
        select(-sp1) %>% 
        as.matrix()  
      
      return(A)
    }
  

  A<-adjacency_matrix(dat,d_cutoff,Lx,Ly)
  
  pr.eigen<-spects(A)$principal_eigenvalue
  skew    <-spects(A)$asymmetry
  peak    <-spects(A)$peakedness
  eigen.gap<-spects(A)$eigengap
  
  mains<-c(pr.eigen,skew,peak,eigen.gap)
  
  library(parallel)
  library(foreach)
  library(iterators)
  library(doParallel)

  
  numcores<-detectCores()
  clust<-makeCluster(min(24,numcores))
  clusterExport(clust,c("dat","adjacency_matrix","spects","d_cutoff","Lx","Ly"))
  registerDoParallel(clust)

  res<-foreach(i=1:100,.combine=rbind,.packages=c("RPANDA","igraph","tidyverse","furrr","RANN"))%dopar%
  {
    nullset<-
      dat%>%mutate(sp=sample(sp))%>%
      adjacency_matrix(d_cutoff=d_cutoff,
                       Lx=Lx,Ly=Ly)%>%spects()
    
    return(c(nullset$principal_eigenvalue,nullset$asymmetry,nullset$peakedness,nullset$eigengap))
  }

  
  res<-rbind(mains,res)
  
  write.csv(res,"nullspect.csv")
