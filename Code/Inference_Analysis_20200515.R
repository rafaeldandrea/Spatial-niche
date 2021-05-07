
## ================ Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc <- !(user <- Sys.info()["user"]) %in% c("rafael", "wangdz", "rdandrea-dw9")

## ========== Load functions =========
loadfuns <- {
  lapply(
    list.files(
      path = "~/R_Functions",
      pattern = "[.][R]$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ),
    source
  )
}

## ======== Libraries =========
library(tidyverse)
library(tibble)
library(vegan)
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps
library(fpc) ## for clustering analysis
library(numbers) ## for function divisors()
library(MASS) ## for generalized inverse (left inverse) of Cmatrix when S > NST
library(dbscan) ## for HDBSCAN clustering method
library(pracma) ## for function Mode()
library(plyr) ##
library(gdata) ## for function upperTriangle() called in function mean_cosine_rows()


## =========== Set directory =============
setwd('~/SpatialNiche/Data/20200515--NST20S20-H-normalizednoise-biotime1e6')


## ============ Read data =============
ind = make_index(1)
filename = paste0('scenario_', seq(3, 60, by = 3), '_run_1.RData')[ind]
dat = get(load(filename))
censuses = dat$censuses
Landscape = dat$Landscape

## ========= Calculate Cij based on the true Nij ==================
inferred_C_optimal <- 
  apply(censuses[, -1], 2, function(community){
    Nij = as.matrix(table(community, Landscape))
    presences = rownames(Nij)
    absences = as.character(setdiff(1:20, presences))
    Nij = rbind(Nij, matrix(0, nrow = length(absences), ncol = ncol(Nij)))
    rownames(Nij) = c(presences, absences)
    Nij = Nij[gtools::mixedorder(rownames(Nij)), ]
    as.numeric(Nij / (1e-16 + rowSums(Nij)))
  }) %>%
  as_tibble


## ======== Calculate Cij based on inferred Nij via MLE =============
clustering.method = 'hierarchical'
inferred_C_mle <-
  apply(censuses[, -1], 2, function(Community){
    A = length(Community)
    
    df_loglik <-
      Get_Likelihood(
        clustering.method = clustering.method,
        Community = Community,
        A = length(Community),
        S = length(unique(Community)),
        minNST = 3,
        maxNST = 30,
        NC = A / numbers::divisors(sqrt(A)) ^ 2,
        verbose = !hpcc
      )
  
    # get inferred landscape and affinity matrix via max likelihood
    df_max_loglikelihood <-
      Max_Loglikelihood(
        df_loglik = df_loglik,
        clustering.method = clustering.method,
        override_bestNST = FALSE,
        Community = Community,
        community_size = length(Community)
      )
  
    list2env(df_max_loglikelihood, envir = environment())
    
    presences = rownames(est_C)
    absences = as.character(setdiff(paste0('species_',1:20), presences))
    est_C = rbind(est_C, matrix(0, nrow = length(absences), ncol = ncol(est_C)))
    rownames(est_C) = c(presences, absences)
    est_C = est_C[gtools::mixedorder(rownames(est_C)), ]
    as.numeric(est_C)
  })

## ================== Analysis ===================
true_C = dat$C
S = nrow(true_C)
NST = ncol(true_C)

RV_optimal <-
  apply(inferred_C_optimal, 2, function(num){
    inferred_C = matrix(num, S, NST)
    observed_species = seq(S)[rowSums(inferred_C) > 0]
    RV_coefficient(true_C[observed_species, ], inferred_C[observed_species, ])$RV
  })

RV_mle <-
  sapply(inferred_C_mle, function(num){
    inferred_C = matrix(num, nrow = S)
    observed_species = seq(S)[rowSums(inferred_C) > 0]
    observed_soiltypes = which(colSums(inferred_C) > 0)
    RV_coefficient(
      true_C[observed_species, ], 
      inferred_C[observed_species, observed_soiltypes]
    )$RV
  })

dat$analysis <-
  list(
    niche_index = 1 - mean_cosine_rows(true_C),
    inferred_C_optimal = inferred_C_optimal,
    inferred_C_mle = inferred_C_mle,
    RV_optimal = RV_optimal,
    RV_mle = RV_mle
  )

## =========== SAVE ==============
if(hpcc) save(dat, file = filename)

## ========= PLOTTING ==========
if(FALSE){
  filenames = paste0('scenario_', seq(3, 60, by = 3), '_run_1.RData')
  df = NULL
  for(char in filenames){
    dat = get(load(char))
    ni = dat$analysis$niche_index
    RV_optimal = dat$analysis$RV_optimal
    RV_mle = unlist(dat$analysis$RV_mle)
    df = rbind(
      df, 
      tibble(
        niche_index = ni, 
        mrvoptim = mean(RV_optimal),
        srvoptim = sd(RV_optimal),
        nrvoptim = length(RV_optimal),
        mrvmle = mean(RV_mle),
        srvmle = sd(RV_mle),
        nrvmle = length(RV_mle)
        
      )
    )
  }
}

plot <-
  df %>%
  filter(niche_index < 1) %>%
  ggplot() +
  geom_errorbar(
    aes(
      x = niche_index, 
      ymin = mrvoptim - srvoptim, 
      ymax = mrvoptim + srvoptim
    ),
    color = red,
    width = .02
  ) +
  geom_smooth(
    aes(niche_index, mrvoptim), 
    color = red, 
    span = 5, 
    method = 'loess',
    se = FALSE
  ) +
  geom_point(aes(niche_index, mrvoptim), color = red) +
  geom_errorbar(
    aes(
      x = niche_index, 
      ymin = mrvmle - srvmle, 
      ymax = mrvmle + srvmle
    ),
    color = blue,
    width = .02
  ) +
  geom_smooth(
    aes(niche_index, mrvmle), 
    color = blue, 
    span = 5, 
    method = 'loess',
    se = FALSE
  ) +
  geom_point(aes(niche_index, mrvmle), color = blue) +
  labs(x = 'Niche index', y = 'mean RV') +
  coord_cartesian(ylim=c(0, 1))


nullrv = sapply(seq(1e4), function(x){
  m1 = matrix(runif(400), 20, 20)
  m2 = matrix(runif(400), 20, 20)
  RV_coefficient(m1, m2)$RV
} )

nullcos = sapply(seq(1e5), function(x){
  v1 = runif(400) #; v1 = v1 - mean(v1)
  v2 = runif(400) #; v2 = v2 - mean(v2)
  v1 %*% v2 / sqrt((v1 %*% v1) * (v2 %*% v2))
} )

dfnull = tibble(RV = nullrv)
plotnull <-
  dfnull %>%
  ggplot() +
  geom_histogram(aes(x = RV), bins = 70) +
  geom_vline(aes(xintercept = quantile(RV, c(.5))), color = red) +
  geom_vline(aes(xintercept = quantile(RV, c(.95))), color = red) +
  ggtitle(paste('Null RV quantiles:            50% ->', round(quantile(dfnull$RV, .5),2), 
                '            95% ->', round(quantile(dfnull$RV, .95), 2)))
## ========== EXTRA CODE ============
if(FALSE){
  
  ## -----------------------------
  plot <-
    res %>%
    ggplot(aes(niche_index, RV)) +
    geom_smooth(span = 1.5) +
    geom_point() +
    ylim(c(0,1)) +
    labs(x = 'Niche index', y = 'RV coefficient') +
    ggtitle('Beta distribution, normalized row sums -- inferred C = Nij / Ni')
  
  
  ## ------------------------------
  foo = tibble(
    N = results$N, 
    cor = sapply(1:20, function(i) cor(results$true_C[i, ], results$inferred_C[i, ]))
  )
  
  bar = as_tibble(sapply(1:1e5, function(i) RV_coefficient(matrix(runif(400), 20, 20), matrix(runif(400), 20, 20))$RV))
  bar %>% 
    ggplot() + 
    geom_histogram(aes(value), bins = 200) + 
    xlim(c(0,1)) +
    geom_vline(aes(xintercept = quantile(value, .95)), color = red)
  
  foo %>% ggplot(aes(N, cor)) +
    geom_point() +
    labs(x = 'Abundance', y = 'Correlation between respective rows')
  
  ## ----------------------------
  setwd('~/SpatialNiche/data/20200515--NST20S20-H-normalizednoise-biotime1e6/')
  
  lf = paste0('scenario_', seq(3 ,80, 3), '_run_1.rdata')
  
  for(char in lf){
    dat = get(load(char))
    
    data = list()
    data$C = dat$C
    data$year  = dat$times
    data$Landscape = dat$Landscape
    data$censuses = as_tibble(sapply(dat$evolving_community, function(l) unlist(l)))
    colnames(data$censuses) = data$year
    data$censuses = data$censuses %>% rowid_to_column(var = 'cell')
    
    save(data, file = char)  
  }
  
}
