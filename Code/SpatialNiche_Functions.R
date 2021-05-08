

ReorderMatrix <- function(data) {
  order_matrix = ncol(data) + 1 - t(apply(data, 1, rank))
  preferred = t(apply(order_matrix, 1, order))
  mylist = NULL
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      if (!preferred[i, j] %in% mylist) {
        mylist = c(mylist, preferred[i, j])
        break
      }
    }
  }
  if (length(mylist) < ncol(data))
    mylist = c(mylist, setdiff(seq(ncol(data)), mylist))
  data = data[, mylist]
  return(data)
}

## Landscape: vector whose A elements are the soil types of each site ========
Generate_Landscape <- function(seed, landscape_autocor_parm) {
  if (landscape.patches) {
    patches <-
      Patches(
        side = sqrt(A / NC),
        number = NC,
        noise = 0,
        L = sqrt(A)
      )
    set.seed(seed) ## change the parameter 'seed' to obtain different landscapes
    Landscape <-
      sample(NST, replace = TRUE, size = NC)[match(patches, unique(patches))]
  }
  if (landscape.randomfield) {
    ## create a Gaussian random field
    field <-
      RandomField(
        L = sqrt(A),
        rangepar = landscape_autocor_parm,
        sillpar = 1,
        nuggetpar = 0,
        seed = seed
      )
    
    ## simplify the field by boxing sites into one of NST soil types
    Landscape <-
      1 + findInterval(field, quantile(field, seq(NST) / NST), rightmost.closed = TRUE)
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
Generate_Cmatrix <- function(seed = 0, nrow = S, ncol = NST, noise = 0) {
  
    set.seed(seed)
    
    ## calculate matrix coordinates of the preferred resource of each species
    preferred <- pos2coord(coord = cbind(1:nrow, 1 + (1:nrow - 1) %% ncol), dim.mat = c(nrow, ncol))
    
    ## generate C matrix via beta distribution with shape parameter inversely proportional to noise
    C <- matrix(rbeta(nrow * ncol, shape1 = 1, shape2 = 1 / (1e-18 + noise)), nrow, ncol)
    
    ## round values to meaningful quantities
    C <- round(C,5)
    
    ## fix affinity for preferred resources to 1
    C[preferred] = 1
    
    return(C)
  }


## community matrix [sqrt(A) x sqrt(A)] =======
Generate_Community <- function(Landscape, Cmatrix, seed) {
  set.seed(seed)
  R <- table(Landscape)
  expected_N <- as.numeric(Cmatrix %*% R)
  Community <-
    sapply(seq(A),
           function(site) {
             sample(S, size = 1, prob = expected_N * Cmatrix[, Landscape[site]])
           })
  return(Community)
}


## community matrix [sqrt(A) x sqrt(A)] =======
Generate_Reference_Community <- function(Landscape, Cmatrix) {
  Community <- sapply(seq(A), function(site) seq(S)[which.max(Cmatrix[,Landscape[site]])])
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
      ungroup %>%
      arrange(cell, species) %>%
      pivot_wider(names_from = 'species', values_from = 'n', names_prefix = 'species_') %>%
      replace(is.na(.),0) %>%
      select(cell,paste0('species_',1:length(unique(Community))))
  
  return(Tabulation)
}


## estimate Cmatrix matrix =========
Generate_estC <- function(est_Nij, reorder = FALSE) {
  tabulation = as.matrix(est_Nij[ , -1])
  soil_abuns = rowSums(tabulation)
  est_C <- t(tabulation / soil_abuns)
  if(reorder) est_C <- ReorderMatrix(est_C)
  return(est_C)
}


## find the point in a curve with the maximum |2nd derivative| ====
Max2Deriv <- function(x, y) {
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which.max(abs(d2))
  return(indices)
}


## ======== get loglikelihood as a function of No. of cells and No. of soil types =============
Get_Likelihood <- 
  function(
      clustering.method, 
      Community, 
      A = length(Community), 
      S = length(unique(Community)),
      minNST = 3,
      maxNST = 30,
      NC = A / numbers::divisors(sqrt(A))^2, 
      verbose
  ){
    ddply(
      crossing(tibble(NC = NC), tibble(NST = minNST:maxNST)) %>% 
        filter(NC > S & NC <= A/9 & NC > 25), .(NST, NC), function(v) {
          NST <- v$NST
          NC <- v$NC
          if(verbose) print(paste(NST,NC))
          cellgrid <- Generate_Grid(A, NC)
          Tabulation <- Generate_Tabulation(cellgrid, Community)
          abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
          tN <- t(abunds_by_cells)
          tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
          
          ## -- Estimate matrix Nij based on the identification of cells to soil types:
          
          ## Using k-means
          if(clustering.method=='kmeans'){
            clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
            cluster_index=clustermodel$cluster
            clustering.statistic=clustermodel$tot.withinss
            numclusters_cells <- NST
          }
          
          ## Using hierarchical clustering
          if(clustering.method=='hierarchical'){
            dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
            clus_sites <- hclust(dis_sites)
            cluster_index <- cutree(clus_sites, NST)
            clustering.statistic <- NA
            numclusters_cells <- NST
          }
          
          ## Using pam clustering
          if(clustering.method=='pam'){
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
          
          est_C <- Generate_estC(est_Nij)
          Pij <- t(t(N*est_C)/colSums(N*est_C))
          
          
          logprob=log(Pij[, cluster_index])
          loglik <- sum(tN*logprob*is.finite(logprob),na.rm=TRUE)
          
          return(
            data.frame(
              NC = NC,
              numclusters = numclusters_cells,
              loglik = loglik,
              statistic = clustering.statistic
            )
          )
        }
    ) %>%
    mutate(cellsize = sqrt(A / NC)) %>%
    as_tibble
  }



## ========= get est_Nij and est_C via max likelihood ===========
Max_Loglikelihood <- function(df_loglik, clustering.method){
  
  df_maxloglik <-
    df_loglik %>% 
    group_by(NST) %>%
    summarize(
      bestcellsize=cellsize[which.max(loglik)],
      maxloglik=loglik[which.max(loglik)]
    )
  
  nst_continuum <- with(df_maxloglik,seq(min(NST),max(NST),l=1e3))
  
  smoothed_maxloglik <- with(df_maxloglik, predict(loess(maxloglik~NST),nst_continuum))
  
  index <- Max2Deriv(nst_continuum, smoothed_maxloglik)
  
  bestNST <- unique(round(nst_continuum[index]))
  
  bestcellsize <- df_maxloglik %>% filter(NST==bestNST) %>% pull(bestcellsize)
  
  nst=df_maxloglik$NST
  
  NC <- A/bestcellsize^2
  
  NST <- bestNST
  
  cellgrid <- Generate_Grid(A, NC)
  
  Tabulation <- Generate_Tabulation(cellgrid, Community)
  
  abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
  
  tN <- t(abunds_by_cells)
  
  tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
  
  ## Using k-means
  if(clustering.method=='kmeans'){
    clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
    cluster_index=clustermodel$cluster
  }
  
  ## Using hierarchical clustering
  if(clustering.method=='hierarchical'){
    dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
    clus_sites <- hclust(dis_sites)
    cluster_index <- cutree(clus_sites, NST)
  }
  
  ## Using pam clustering
  if(clustering.method=='pam'){
    dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
    clustermodel=cluster::pam(x=dis_sites, k=NST)
    cluster_index <- clustermodel$clustering
  }
  
  est_Nij <-
    Tabulation %>%
    mutate(soiltype = cluster_index) %>%
    select(-cell) %>%
    group_by(soiltype) %>%
    summarize_all(sum) %>%
    filter(soiltype > 0)
  
  est_C <- Generate_estC(est_Nij, reorder = FALSE)
  
  return(
    list(
      Tabulation = Tabulation,
      est_NC = NC,
      est_NST = NST,
      est_Nij = est_Nij,
      est_C = est_C,
      cluster_index = cluster_index
    )
  )
}

## -- Simulate Community Dynamics
Simulate_Dynamics <- {
  function(
    Landscape, 
    Cmatrix, 
    Community,
    annual_turnover_rate,                 ## percentage of the forest that gets replaced per year
    delta_t = 1,                          ## default: one year
    biotime = 0,                          ## initial time, in years
    maxtime = 1e3,                        ## number of years to run for
    interval = 1 / annual_turnover_rate,  ## interval in years for recording current state
    seed
  ){
    stopifnot(is.vector(Community))
    set.seed(seed)
    A = length(Landscape)
    S = length(unique(as.numeric(Community)))
    N = sapply(1:S, function(species) sum(Community == species))
    df_N = N
    names(df_N) = paste0('species_',seq(S))
    total_deaths = 0
    times = 0
    cum_deaths = 0
    rectime = 0
    initial_community = Community
    
    while(biotime <= maxtime) {
      if(mod(biotime,10)==0){
        print(biotime)
      }
      # print(paste(biotime, min(N)), quote = FALSE)
      biotime = biotime + delta_t
      numdeaths = min(2* annual_turnover_rate * A * delta_t, rpois(1, lambda = annual_turnover_rate * A * delta_t))
      vacancies = sample(A, size = numdeaths)
      C0 = Community
      stopifnot(length(unique(C0)) == S)
      Community[vacancies] = sapply(vacancies, function(site) sample(S, size = 1, prob = N * Cmatrix[, Landscape[site]]))
      N = sapply(1:S, function(species) sum(Community == species))
      extinct = which(N == 0)
      while(length(extinct) > 0){
        for(species in extinct){
          original_sites = which(C0 == species)
          restored_site = ifelse(length(original_sites) > 1, sample(original_sites, size = 1), original_sites)
          Community[restored_site] = species
        }
        N = sapply(1:S, function(species) sum(Community == species))
        extinct = which(N == 0)
        
      }  
      total_deaths = total_deaths + numdeaths
      if ((biotime - rectime) > interval) {
        df_N = rbind(df_N, N, deparse.level = 0)
        cum_deaths = c(cum_deaths, total_deaths)
        times = c(times, round(biotime))
        rectime = biotime
        total_deaths = 0
      }
    }
    
    return(
      list(
        Landscape = Landscape, 
        C = Cmatrix, 
        initial_community = initial_community,
        annual_turnover_rate = annual_turnover_rate,                 
        delta_t = delta_t,                          ## default: one year
        biotime = biotime,                          ## initial time, in years
        maxtime = maxtime,                        ## number of years to run for
        interval =interval,  ## interval in years for recording current state
        seed = seed,
        Final_Community = unlist(Community),
        abundances = as_tibble(df_N),
        times = times,
        new_deaths = cum_deaths
      )
    )
  }
}


## calculate difference between true and inferred matrix

# RV index -- cosine between Cmatrix and est_C =========
RV_coefficient <- function(matrix_1, matrix_2){
  stopifnot(is.matrix(matrix_1), is.matrix(matrix_2))
  tr <- function(m) sum(diag(m))
  squarify <- function(m) m %*% t(m)
  center <- function(m) t(t(m) - colMeans(m))
  CC <- squarify(center(matrix_1))
  eCeC <- squarify(center(matrix_2))
  return(list(RV = tr(CC %*% eCeC) / sqrt(tr(CC %*% CC) * tr(eCeC %*% eCeC))))
}
## mean cosine between matrix rows
mean_cosine_rows <- function(matrix){
  C <- matrix
  M <- C %*% t(C)
  T <- diag(M)
  Q <- sqrt(outer(T, T))
  costheta <- as.numeric(gdata::upperTriangle(M / Q))
  return(mean(costheta))
}

## mean cosine between matrix columns
mean_cosine_columns <- function(matrix){
  C <- matrix
  M <- t(C) %*% C
  T <- diag(M)
  Q <- sqrt(outer(T, T))
  costheta <- as.numeric(gdata::upperTriangle(M / Q))
  return(mean(costheta))
}
