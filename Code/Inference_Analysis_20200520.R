
calculate_inferred_C_optimal_and_inferred_C_mle = FALSE
calculate_inferred_C_mle_references = FALSE
generate_reference_community = FALSE
calculate_inference_referencecommunities_clusteringdata = FALSE
calculate_RV_references = FALSE
calculate_inference_evolvingcommunity_clusteringdata = FALSE
calculate_RV_evolvingcommunity_noinference = FALSE
calculate_RV_averagedcommunities_noinference = FALSE
calculate_RV_averagedcommunities_gridding = FALSE
analysis_gridding_vs_noinference = FALSE
assess_best_number_of_clusters = TRUE

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
setwd('~/SpatialNiche/Data/20200530--NST20S20-H-allnoises-biotime1e7')


## ============ Read data =============
ind = make_index(1)
filename = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))[ind]
dat = get(load(filename))
censuses = dat$evolving_community
Landscape = dat$Landscape

if(calculate_inferred_C_optimal_and_inferred_C_mle){
  ## ========= Calculate Cij based on the true Nij ==================
  inferred_C_optimal <- 
    sapply(censuses, function(community){
      Nij = as.matrix(table(community, Landscape))
      presences = rownames(Nij)
      absences = as.character(setdiff(1:20, presences))
      Nij = rbind(Nij, matrix(0, nrow = length(absences), ncol = ncol(Nij)))
      rownames(Nij) = c(presences, absences)
      Nij = Nij[gtools::mixedorder(rownames(Nij)), ]
      as.numeric(Nij / (1e-16 + rowSums(Nij)))
    }) %>%
    as_tibble
  
  colnames(inferred_C_optimal) = dat$times

  ## ======== Calculate Cij based on inferred Nij via MLE =============
  clustering.method = 'hierarchical'
  inferred_C_mle <-   ## Note that inferred_C_mle is a list because not all censuses will lead to the same estimated NST
    lapply(censuses, function(Community){
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
}

if(generate_reference_community){
  setwd('~/SpatialNiche/Data/20200520--NST20S20-H-allnoises-biotime1e7')
  filenames = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))
  
  ## ========== Generate Reference Community ============
  for(char in filenames){
    writeLines(char)
    dat = get(load(char))
    dat$reference_communities <-
      sapply(seq(100), function(theseed){
        S = 0
        while(S < 20){
          community = Generate_Community(Landscape = dat$Landscape, Cmatrix = dat$C, seed = theseed)
          S = length(unique(community))
        }
        return(community)
      })
    save(dat, file = char)
  }
}

if(calculate_inferred_C_mle_references){
  clustering.method = 'hierarchical'
  
  dat$inferred_C_mle_references <-   ## Note that inferred_C_mle is a list because not all censuses will lead to the same estimated NST
    apply(dat$reference_communities, 2, function(Community){
      A = length(Community)
      
      df_loglik <-
        Get_Likelihood(
          clustering.method = clustering.method,
          Community = Community,
          minNST = 3,
          maxNST = 30,
          mincellsize = 4,
          maxcellsize = sqrt(A) / 6,
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
    })
  
  save(dat, file = filename)
  
}

if(calculate_inference_referencecommunities_clusteringdata){
  clustering.method = 'hierarchical'
  
  thelist <-  
    apply(dat$reference_communities, 2, function(Community){
      res = Infer_C(Community = Community, true_C = dat$C, mincellsize = 4, maxcellsize = 120 / 6, verbose = !hpcc)
      return(res)
    })
  
  nrows = nrow(thelist[[1]] %>% select(NST, cellsize) %>% unique)
  
  dat$inference_referencecommunities_clusteringdata <-
    thelist %>%
    reduce(full_join) %>%
    mutate(reference_index = rep(seq(ncol(dat$reference_communities)), each = nrows))
      
  if(hpcc) save(dat, file = filename)
}

if(calculate_inference_evolvingcommunity_clusteringdata){
  clustering.method = 'hierarchical'
  
  thelist <-  
    lapply(dat$evolving_community, function(Community){
      res = Infer_C(Community = Community, true_C = dat$C, mincellsize = 4, maxcellsize = 120 / 6, verbose = !hpcc)
      return(res)
    })
  
  nrows = nrow(thelist[[1]] %>% select(NST, cellsize) %>% unique)
  
  dat$inference_evolvingcommunity_clusteringdata <-
    thelist %>%
    reduce(full_join) %>%
    mutate(year = rep(seq(length(dat$evolving_community)), each = nrows))
  
  if(hpcc) save(dat, file = filename)
}

if(calculate_RV_references){
  for(filename in filenames){
    print(filename)
    dat = get(load(filename))
    true_C = dat$C
    S = nrow(true_C)
    NST = ncol(true_C)
    
    foo = if(!is.list(dat$inferred_C_mle_references)) list(dat$inferred_C_mle_references) else dat$inferred_C_mle_references
    
    dat$RV_mle_references <-
      sapply(foo, function(num){
        inferred_C = num
        observed_species = seq(S)[rowSums(inferred_C) > 0]
        observed_soiltypes = which(colSums(inferred_C) > 0)
        RV_coefficient(
          true_C[observed_species, ],
          inferred_C[observed_species, observed_soiltypes]
        )$RV
      })

    save(dat, file = filename)
  }
}

if(calculate_RV_evolvingcommunity_noinference){
  filenames = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))
  for(filename in filenames){
    dat = get(load(filename))
    dat$RV_evolvingcommunity_noinference = sapply(dat$evolving_community,function(Community){
      Nij = table(Community, dat$Landscape)
      Cij = Nij / rowSums(Nij)
      RV = RV_coefficient(Cij, dat$C[match(rownames(Cij), 1:20), ])$RV
    })
    save(dat, file = filename)
  }
}

if(calculate_RV_averagedcommunities_noinference){
  filenames = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))
  for(filename in filenames){
    writeLines(filename)
    dat = get(load(filename))
    Nij = matrix(0,20,20)
    for(i in 2:length(dat$evolving_community)){
      foo = table(dat$evolving_community[[i]], dat$Landscape)
      presences = rownames(foo)
      absences = setdiff(1:20, presences)
      foo = rbind(foo, matrix(0,length(absences), 20))
      rownames(foo) = c(presences, absences)
      foo = foo[gtools::mixedorder(rownames(foo)),]
      Nij = Nij + foo
    } 
    Cij = Nij / rowSums(Nij)
    
    dat$RV_averagedcommunities_noinference <- RV_coefficient(Cij, dat$C)$RV
    
    save(dat, file = filename)
  }
}

if(calculate_RV_averagedcommunities_gridding){
  A = 120 ^ 2
  cellsize = 4
  NC = A / cellsize ^ 2
  NST = 20
  linkage = 'ward.D2'
  cellgrid = Generate_Grid(A = A, NC = NC)
  filenames = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))
  
  for(filename in filenames){
    writeLines(filename)
    dat = get(load(filename))
    cell_by_species = matrix(0, NC, 20)
    
    for(i in 2:length(dat$evolving_community)){
      Community = dat$evolving_community[[i]]
      tabula = Generate_Tabulation(cellgrid = cellgrid, Community = Community)
      presences = match(colnames(tabula)[-1], paste0('species_', 1:20))
      absences = setdiff(1:20, presences)
      if(length(absences) > 0){
        tabula = cbind(tabula, matrix(0, NC, length(absences)))
        colnames(tabula) = c('cell', paste0('species_', presences), paste0('species_', absences))
      }
      tabula = select(tabula, paste0('species_', 1:20))
      cell_by_species = cell_by_species + tabula
    } 
    
    
    tN <- t(cell_by_species)
    tN_normalized <- tN / rowSums(tN) ## species x cells
    dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
    clus_sites <- hclust(dis_sites, method = linkage)
    cluster_index <- cutree(clus_sites, NST)
    
    
    ## tabulation of species into all cells of soil type j
    est_Nij <-
      cell_by_species %>%
      mutate(soiltype = cluster_index) %>%
      group_by(soiltype) %>%
      summarize_all(sum) %>%
      filter(soiltype > 0)
    
    est_C <- Estimate_C3(est_Nij)
    
    dat$est_C_averagedcommunities_gridding <- est_C
    
    dat$RV_averagedcommunities_gridding <- RV_coefficient(est_C, dat$C)$RV
    
    save(dat, file = filename)
  }
}

if(assess_best_number_of_clusters){
  A = 120 ^ 2
  cellsize = 4
  NC = A / cellsize ^ 2
  linkage = 'ward.D2'
  cellgrid = Generate_Grid(A = A, NC = NC)
  filenames = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))
  
  par(mfrow = c(2,10), mar = c(2,2,2,2))
  
  for(filename in filenames[seq(3,60,by=3)]){
    writeLines(filename)
    dat = get(load(filename))
    cell_by_species = matrix(0, NC, 20)
    
    for(i in 2:length(dat$evolving_community)){
      Community = dat$evolving_community[[i]]
      tabula = Generate_Tabulation(cellgrid = cellgrid, Community = Community)
      presences = match(colnames(tabula)[-1], paste0('species_', 1:20))
      absences = setdiff(1:20, presences)
      if(length(absences) > 0){
        tabula = cbind(tabula, matrix(0, NC, length(absences)))
        colnames(tabula) = c('cell', paste0('species_', presences), paste0('species_', absences))
      }
      tabula = select(tabula, paste0('species_', 1:20))
      cell_by_species = cell_by_species + tabula
    } 
    
    tN <- t(cell_by_species)
    tN_normalized <- tN / rowSums(tN) ## species x cells
    dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
   # clus_sites <- hclust(dis_sites, method = linkage)
    
    result <- 
      NbClust::NbClust(
        data = cell_by_species, 
        diss = dis_sites, 
        distance = NULL, 
        min.nc = 3, 
        max.nc = 30, 
        method = 'ward.D2', 
        index = 'sdbw'
      )
    plot(3:30, result$All.index, t='l', main = round(dat$niche_index, 2))
  }
    
  #   idx <- 
  #     sapply(minNST:maxNST, function(NST){
  #       cluster_index <- cutree(clus_sites, NST)
  #       unlist(intCriteria(as.matrix(cell_by_species / rowSums(cell_by_species)), cluster_index, 'all'))
  #     })
  #   index = rownames(idx)
  #   idx = as_tibble(idx)
  #   colnames(idx) = minNST:maxNST
  #   idx = idx %>% mutate(index = index)
  #   
  #   ## tabulation of species into all cells of soil type j
  #   est_Nij <-
  #     cell_by_species %>%
  #     mutate(soiltype = cluster_index) %>%
  #     group_by(soiltype) %>%
  #     summarize_all(sum) %>%
  #     filter(soiltype > 0)
  #   
  #   est_C <- Estimate_C3(est_Nij)
  #   
  #   dat$est_C_averagedcommunities_gridding <- est_C
  #   
  #   dat$RV_averagedcommunities_gridding <- RV_coefficient(est_C, dat$C)$RV
  #   
  #   save(dat, file = filename)
  # }
}

if(analysis_gridding_vs_noinference){
  filenames = gtools::mixedsort(paste0('scenario_',seq(3, 60, by = 3),'_run_1.rdata'))
  res = ddply(tibble(filename = filenames), .(filename), function(v){
    filename = v$filename
    print(filename)
    dat = get(load(filename))
    
    if(!is.null(dat$inference_referencecommunities_clusteringdata)){
      foo_reference <- 
        dat$inference_referencecommunities_clusteringdata %>% 
        mutate(data = 'reference') %>% 
        rename(year = 'reference_index') 
    } else foo_reference = NULL
    if(!is.null(dat$inference_evolvingcommunity_clusteringdata)){
      foo_dynamical <- 
        dat$inference_evolvingcommunity_clusteringdata %>% 
        mutate(data = 'gridding')
    }
    
    foo = rbind(foo_reference, foo_dynamical) %>% mutate(niche_index = dat$niche_index, noise_type = dat$noise_type)
    
    if(!is.null(foo)){
      res <-  
        foo %>% 
        group_by(data, year) %>% 
        #filter(NST == 20) %>%
        filter(RV == max(RV)) %>%
        filter(loglik == max(loglik)) %>%
        ungroup
      return(res)
    }
  }) %>% as_tibble
  
  res_summary <-
    res %>%
    group_by(noise_type, data, niche_index) %>%
    summarize(mRV = mean(RV), sRV = sd(RV)) %>%
    ungroup
  
  df_noinf <- 
    ddply(tibble(filename = filenames), .(filename), function(df){
      filename = df$filename
      dat = get(load(filename))
      mRV = mean(dat$RV_evolvingcommunity_noinference)
      sRV = sd(dat$RV_evolvingcommunity_noinference)
      return(tibble(mRV = mRV, sRV = sRV, niche_index = dat$niche_index))
    }) %>% 
    mutate(data = 'noinference', noise_type = 'normalized') %>%
    select(-filename) %>%
    as_tibble
  
  df_avg_noinf <- 
    ddply(tibble(filename = filenames), .(filename), function(df){
      filename = df$filename
      dat = get(load(filename))
      mRV = dat$RV_averagedcommunities_noinference
      return(tibble(mRV = mRV, sRV = 0, niche_index = dat$niche_index))
    }) %>% 
    mutate(data = 'noinference_avg', noise_type = 'normalized') %>%
    select(-filename) %>%
    as_tibble
  
  df_avg_gridding <- 
    ddply(tibble(filename = filenames), .(filename), function(df){
      filename = df$filename
      dat = get(load(filename))
      mRV = dat$RV_averagedcommunities_gridding
      return(tibble(mRV = mRV, sRV = 0, niche_index = dat$niche_index))
    }) %>% 
    mutate(data = 'gridding_avg', noise_type = 'normalized') %>%
    select(-filename) %>%
    as_tibble
  
  df <-
    res_summary %>%
    union(df_noinf) %>%
    union(df_avg_noinf) %>%
    union(df_avg_gridding)
  
  plot <-
    df %>% 
    filter(data != 'reference', niche_index < 1) %>%
    ggplot(aes(niche_index, mRV, color = data)) + 
    geom_errorbar(aes(ymin = mRV - sRV, ymax = mRV + sRV), width = .02) +
    geom_line() +
    geom_point() + 
    ylim(c(0,1)) +
    labs(x = 'Niche index', y = 'mean RV coefficient', color = 'Analysis') +
    ggtitle('Normalized beta noise')
}



## ========= PLOTTING ==========
if(FALSE){
  filenames = gtools::mixedsort(list.files(pattern = '*_run_1.RData'))
  df = NULL
  for(char in filenames){
    dat = get(load(char))
    ni = dat$analysis$niche_index
    RV_optimal = dat$analysis$RV_optimal
    RV_mle = unlist(dat$analysis$RV_mle)
    RV_mle_ideal_community = dat$analysis$RV_mle_ideal_community
    RV_mle_reference = dat$RV_mle_references
    df = rbind(
      df, 
      tibble(
        scenario = dat$scenario,
        niche_index = ni, 
        mrvoptim = mean(RV_optimal),
        srvoptim = sd(RV_optimal),
        nrvoptim = length(RV_optimal),
        mrvmle = mean(RV_mle),
        srvmle = sd(RV_mle),
        nrvmle = length(RV_mle),
        mrvmle_ideal = RV_mle_ideal_community,
        mrvmle_ref = mean(RV_mle_reference),
        srvmle_ref = sd(RV_mle_reference)
      )
    )
  }
  df <- 
    df %>%
    mutate(noise_type = rep(c('beta','normalized'),20))
}

if(!hpcc){
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
    # geom_errorbar(
    #   aes(
    #     x = niche_index, 
    #     ymin = mrvmle - srvmle, 
    #     ymax = mrvmle + srvmle
    #   ),
    #   color = blue,
    #   width = .02
    # ) +
    # geom_smooth(
    #   aes(niche_index, mrvmle), 
    #   color = blue, 
    #   span = 5, 
    #   method = 'loess',
    #   se = FALSE
    # ) +
    # geom_point(aes(niche_index, mrvmle), color = blue) +
    geom_errorbar(
      aes(
        x = niche_index, 
        ymin = mrvmle_ref - srvmle_ref, 
        ymax = mrvmle_ref + srvmle_ref
      ),
      color = green,
      width = .02
    ) +
    geom_smooth(
      aes(niche_index, mrvmle_ref), 
      color = green, 
      span = 5, 
      method = 'loess',
      se = FALSE
    ) +
    geom_point(aes(niche_index, mrvmle_ref), color = green) +
    labs(x = 'Niche index', y = 'mean RV') +
    coord_cartesian(ylim=c(0, 1)) +
    facet_wrap(~noise_type)
  
  
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
  
}
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
