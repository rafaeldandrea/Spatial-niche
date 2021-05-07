## =============== libraries ==========
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


## ====== data folder ===========
data_date <- '20200409'

count <- dplyr::count
select <- dplyr::select
summarize <- dplyr::summarize
red <- "#DD5144"

## ================ Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc <- !(user <- Sys.info()["user"]) %in% c("rafael", "wangdz", "rdandrea-dw9")

## ========== Load functions =========
loadfuns <- {
  lapply(
    list.files(
      path = switch(
        match(user, c("rafael", "wangdz", "rdandrea", "dw9", "rdandrea-dw9")),
        "~/R_Functions",
        "/Users/wangdz/Documents/UIUC/Odywer/R_Functions",
        "~/R_Functions",
        "/home/dw9/R_Functions",
        "/home/rdandrea/R_Functions"
      ),
      pattern = "[.][R]$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ),
    source
  )
}

## Read scenario from either manual input or command line from Bash file
ind <- make_index(37)

## ========= operation directives ==========
landscape.patches <- FALSE
landscape.randomfield <- TRUE

do.data <- TRUE
do.analysis.sad <- FALSE
do.analysis.ranks <- FALSE
do.analysis.species <- FALSE
do.analysis <- FALSE

scenarios <-
  tibble(
    A = 240 ^ 2,
    clustering.method = c('pam'),
    NST = 80,
    S  = 80,
    landscape_autocor_parm  = 50,
    noise  =  c(
      seq(0, .25, l = 15),
      seq(.26, 1, l = 15),
      seq(1.1, 4, l = 15),
      seq(4.1, 30, l = 15),
      100
    )
  )

## Choose row of 'parameters'data frame to specify scenario
scen = scenarios[ind, ]

## Assign each column of scen to a variable matching its name
list2env(as.list(scen), envir = globalenv())

## number of runs per scenario
runlist <- if(hpcc) 1:100 else 10

if (do.data) {
  
  do.references <- TRUE
  
  for (run in runlist) {
    
    Landscape <- Generate_Landscape(seed = run, landscape_autocor_parm = landscape_autocor_parm)
    
    true_C <- Generate_Cmatrix(seed = run, noise = noise)
    
    rownames(true_C)=paste0('species_',seq(S))
    colnames(true_C)=paste0('resource_',seq(NST))
    
    annual_turnover_rate <- .1
    maxtime <- 1e3
    
    seed_initial_community <- run - 1
    initial_community <- NULL
    while(length(unique(initial_community)) != S){
      seed_initial_community <- seed_initial_community + 1
      initial_community <- Generate_Community(Landscape, true_C, seed = seed_initial_community)
    }
    
    initial_species <- sort(unique(initial_community))
    S <- length(unique(initial_species))
    true_C = true_C[initial_species, ]
    
    com <- 
      Simulate_Dynamics(
        Landscape = Landscape, 
        Cmatrix = true_C, 
        Community = initial_community, 
        annual_turnover_rate = annual_turnover_rate, 
        delta_t = 1,                          
        biotime = 0,                          
        maxtime = maxtime,                      
        interval = 1 / annual_turnover_rate,
        seed = run
      )
    
    Community <- com$Final_Community
    
    R <- as.numeric(table(Landscape))
    
    abundance <- 
      table(Community) %>%
      enframe %>% 
      right_join(tibble(name = as.character(initial_species)), by = 'name') %>% 
      replace_na(list(value = 0)) %>%
      pull(value)
    
    ## ======== get loglikelihood as a function of No. of cells and No. of soil types =============
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
    df_maxloglik <- Max_Loglikelihood(df_loglik = df_loglik, clustering.method = clustering.method)
    list2env(df_maxloglik, envir = globalenv())
    colnames(est_C) <- paste0('resource_',seq(ncol(est_C)))
    
    # infer landscape
    linear_cell_size  <- sqrt(size(Landscape)[2] / est_NC)
    foo <- rep(cluster_index, each = linear_cell_size) # repeat every element by cell linear size
    newlandscape <- unlist(rep(split(foo, ceiling(seq_along(foo) / sqrt(A))), each = linear_cell_size), use.names = FALSE)
    
    ## generate reference communities
    references <- NULL
    # if(do.references){
    #   numref=1e3
    #   references <-
    #     lapply(1:numref, function(seed){
    #       table(Generate_Community(newlandscape, est_C, seed = seed)) %>%
    #       enframe %>%
    #       mutate(seed = seed)
    #     }) %>%
    #     bind_rows %>%
    #     rename(
    #       species = name,
    #       abundance = value
    #     )
    # }
    # 
    if(do.references){
      numref=2
      references <-
        lapply(1:numref, function(seed){
          
            final_community <-
              Simulate_Dynamics(
                Landscape = newlandscape, 
                Cmatrix = est_C, 
                Community = Generate_Reference_Community(newlandscape, est_C), 
                annual_turnover_rate = annual_turnover_rate, 
                delta_t = 1,                          
                biotime = 0,                          
                maxtime = maxtime,                      
                interval = 1 / annual_turnover_rate,
                seed = seed
              )$Final_Community
            
            return(
              table(final_community) %>%
              enframe %>%
              mutate(seed = seed)
            )
        }) %>%
        bind_rows %>%
        rename(
          species = name,
          abundance = value
        )
    }
    
    ## save data
    if (hpcc) {
      datadir <- paste0('~/SpatialNiche/Data/', data_date, '/')
      # datadir <- paste0('/data/rdandrea-dw9/SpatialNiche/Data/', data_date, '/')
      dir.create(datadir, showWarnings = FALSE)
      savename = paste0(datadir, 'scenario_', ind, '_run_', run,'.RData')
      dat <-
        list(
            scenario = ind,
            
            run = run,
            
            noise = noise,
            
            Landscape = Landscape,
            
            C = true_C,
            
            mean_cosine_rows = mean_cosine_rows(true_C),
            
            cv = sd(true_C) / mean(true_C),
            
            RV = RV_coefficient(true_C, est_C)$RV,
            
            est_C = est_C,
            
            best_numcells = est_NC,
            
            est_numsoiltypes = est_NST,
            
            est_Landscape = newlandscape,
            
            observed_community =
              table(Community) %>%
              enframe %>%
              mutate(abundance = as.numeric(value)) %>%
              select(-value) %>%
              rename(species = name),
            
            reference_communities = references
        )
      save(dat, file = savename)
    }
  }
}

if (do.analysis.sad) {
  redo.analysis <- FALSE
  
  setwd(paste0('~/spatialniche/data/', data_date, '/'))
  indices = 1:91
  lf = paste0('data_', indices, '.rdata')
  dat = do.call('rbind', lapply(lf, function(l)
    get(load(l))))
  
  if (redo.analysis) {
    cucconi <-
      dat %>%
      group_by(scenario, run) %>%
      summarize(
        p_value =
          cucconi.test(
            observed_community,
            inferred_community,
            method = 'permutation',
            verbose = FALSE
          )$p.value
      ) %>%
      ungroup
    
    kolsmi <-
      dat %>%
      group_by(scenario, run) %>%
      summarize(
        p_value =
          dgof::ks.test(
            observed_community,
            inferred_community,
            exact = FALSE,
            simulate.p.value = TRUE
          )$p.value
      ) %>%
      ungroup
    
    data <-
      cucconi %>%
      mutate(test = 'cucconi') %>%
      union(kolsmi %>%
              mutate(test = 'ks'))
    
  }
  
  if (!redo.analysis)
    data = get(load('two_sample_sad_test_results.rdata'))
  
  analysis <-
    data %>%
    group_by(test, scenario) %>%
    summarize(success = 100 * mean(p_value >= .05),
              samples = n()) %>%
    left_join(scenarios %>% mutate(scenario = seq(n())), by = 'scenario') %>%
    left_join(dat %>%
                filter(run == 1) %>%
                select(scenario, mean_cosine_rows, cv) %>%
                unique)
  
  plot <-
    analysis %>%
    ggplot(aes(1 - mean_cosine_rows, success)) +
    geom_point() +
    theme_bw() +
    labs(x = 'Degree of Niche Differentiation (1 - mean(cosine(rows(C))))',
         y = 'Percentage of indistuinguishable SADs') +
    ylim(c(0, 100)) +
    theme(aspect.ratio = 1) +
    facet_wrap( ~ test, nrow = 1)
  
  plot_cucconi_hierarchical <-
    analysis %>%
    filter(test == 'cucconi' & clustering.method == 'hierarchical') %>%
    ggplot(aes(1 - mean_cosine_rows, success)) +
    # geom_line() +
    geom_point() +
    theme_bw() +
    labs(x = 'Degree of Niche Differentiation (1 - mean(cosine(rows)))',
         y = 'Percentage of indistuinguishable SADs',
         title = 'Cucconi Test') +
    ylim(c(0, 100)) +
    theme(aspect.ratio = 1)
}

if (do.analysis.ranks) {
  setwd('~/SpatialNiche/Data/20200331/')
  
  redo.analysis <- FALSE
  
  if(redo.analysis){
  
    Normalize <- function(x) x / sum(x)
    
    index_df <- crossing(scenario = 1:91, run = 1:100)
    
    scen <- index_df[ind,]
    
    list2env(scen, envir = .GlobalEnv)
          
    if(!hpcc) print(paste(scenario,run))
  
    dat <- get(load(paste0('scenario_', scenario, '_run_', run, '.RData')))
    
    list2env(dat, envir = .GlobalEnv)
    
    obs_N <- 
      observed_community %>%
      pull(abundance)
    
    exp_N <- sum(observed_community$abundance) * Normalize(as.numeric(est_C %*% table(est_Landscape)))
    
    RV <- RV_coefficient(Cmatrix, est_C)$RV
    
    references <-
      references %>% 
      pivot_wider(names_from = species, values_from = abundance) %>% 
      replace(is.na(.), 0) %>%
      pivot_longer(cols = -seed, names_to = 'species', values_to = 'abundance')
    
    G_i <- 2 * obs_N * log(obs_N / exp_N)
    obs_G <- sum(G_i[is.finite(G_i)])
    obs_X2 <- sum((obs_N - exp_N) ^ 2 / exp_N)
    obs_cucconi <- cucconi.test(obs_N, exp_N, method = 'permutation', verbose = FALSE)$C
    
    ref_gof <-
      references %>%
      mutate(obs_N = abundance) %>%
      group_by(seed) %>%
      mutate(
        G_i = 2 * obs_N * log(obs_N / exp_N),
        X2_i = (obs_N - exp_N) ^ 2 / exp_N
      ) %>%
      summarize(
        G = sum(G_i[is.finite(G_i)]),
        X2 = sum(X2_i),
        cucconi = cucconi.test(obs_N, exp_N, method = 'permutation', verbose = FALSE)$C
      )
    
    res <- 
      tibble(
        scenario = scenario,
        run = run,
        mean_cosine_rows = mean_cosine_rows(C),
        mean_cosine_parms = parms$mean_cosine_rows,
        cv = parms$cv,
        pval_G = mean(ref_gof$G > obs_G),
        pval_X2 = mean(ref_gof$X2 > obs_X2),
        pval_cucconi = mean(ref_gof$cucconi > obs_cucconi),
        RV = RV
      )
    
    ## save results
    if (hpcc) {
      datadir = paste0('~/SpatialNiche/Data/20200331/')
      dir.create(datadir, showWarnings = FALSE)
      savename = paste0(datadir, 'analysis_scenario_', scenario, '_run_', run,'.RData')
      save(res, file = savename)
    }
  }

  if(!redo.analysis){
    
    res <- get(load('~/SpatialNiche/Data/20200331/analysis_data.rdata'))
    
    res_summary <-
      res %>%
      group_by(scenario) %>%
      summarize(
        niche_index = 1 - mean(mean_cosine_rows),
        RV = mean(RV),
        G_matches = 100 * mean(pval_G > .05),
        X2_matches = 100 * mean(pval_X2 > .05),
        SAD_matches = 100 * mean(pval_cucconi > .05)
      )
    
    plot <-
      res_summary %>%
      select(-X2_matches) %>%
      pivot_longer(cols = c(G_matches, SAD_matches), names_to = 'test', values_to = 'test_value') %>%
      pivot_longer(cols = c(niche_index, RV), names_to = 'index', values_to = 'index_value') %>%
      ggplot(aes(index_value, test_value)) +
      geom_point() +
      facet_wrap(index~test, scales = 'fixed') +
      theme_bw()
    
    gridExtra::grid.arrange(plot)
  }
}

if(do.analysis.species){
  setwd('~/SpatialNiche/Data/20200403/')
  
  Normalize <- function(x) x / sum(x)
  
  index_df <- crossing(scenario = 1:91, run = 1:100)
  
  foo <-
    ddply(index_df,.(scenario,run),
          function(v){
  
            list2env(v, envir = environment())
            
            filename <- paste0('scenario_', scenario, '_run_', run, '.RData')
            
            if(file.exists(filename)){
              
                writeLines(paste(scenario, run))
              
                dat <- get(load(filename))
              
                list2env(dat, envir = environment())
                
                obs_N <- 
                  observed_community %>%
                  pull(abundance)
                
                exp_N <- sum(observed_community$abundance) * Normalize(as.numeric(est_C %*% table(est_Landscape)))
              
                return(
                  tibble(
                    obs_N = obs_N, 
                    exp_N = exp_N,
                    mean_cosine = parms$mean_cosine_rows
                  )
                )
            }
          }
      ) %>% 
    as_tibble %>%
    mutate(
      obs_N = as.numeric(obs_N),
      niche_index = 1 - mean_cosine)
  
  res <- res %>% mutate(niche_index = 1 - mean_cosine_rows)
  chosen_values <- sapply(quantile(res$scenario, 1:12 / 12), function(.) DescTools::Closest(res$scenario,.)[1])
  
  bar <-
    foo %>%
    group_by(scenario) %>%
    mutate(niche_index_mean = mean(niche_index)) %>%
    ungroup
  
  plot <-
    foo %>%
    filter(scenario==40) %>%
    mutate(niche_index = round(niche_index, 2)) %>%
    ggplot(aes(exp_N, obs_N)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    facet_wrap(~scenario)
  
  gridExtra::grid.arrange(plot)
}

if(do.analysis){
  setwd('~/SpatialNiche/Data/20200408/')
  
  Normalize <- function(x) x / sum(x)
  
  index_df <- crossing(scenario = 1:61, run = 1:100)
  
  foo <-
    ddply(index_df,.(scenario,run),
          function(v){
            
            list2env(v, envir = environment())
            
            filename <- paste0('scenario_', scenario, '_run_', run, '.RData')
            
            if(file.exists(filename)){
              
              writeLines(paste(scenario, run))
              
              dat <- get(load(filename))
              
              list2env(dat, envir = environment())
              
              obs_N <- 
                observed_community %>%
                pull(abundance)
              
              exp_N <- as.numeric(est_C %*% table(est_Landscape))
              
              return(
                tibble(
                  obs_N = as.numeric(obs_N), 
                  exp_N = exp_N,
                  mean_cosine_rows = mean_cosine_rows,
                  cv = cv,
                  RV = RV,
                  niche_index = 1 - mean_cosine_rows,
                  cor_n = cor(obs_N, exp_N)
                )
              )
            }
          }
    ) %>% 
    as_tibble
  
  bar <-
    foo %>% 
    group_by(scenario) %>%
    summarize_at(c('cor_n','niche_index','RV','cv'),mean) %>%
    inner_join(
      foo %>%
      group_by(scenario) %>%
      summarize_at('obs_N', var),
      by = 'scenario'
    )
  
  res <-
    foo %>%
    ddply(.(scenario, run, niche_index), function(v){
      list2env(v, envir = environment())
      print(paste(unique(scenario), unique(run)))
      gtest <- try(DescTools::GTest(x = obs_N, p = exp_N/sum(exp_N)), silent = TRUE)
      cucconi <- cucconi.test(obs_N, exp_N, verbose = FALSE)
      return(
        tibble(
          G = 2 * sum(obs_N * log(obs_N / exp_N)),
          G.stat = gtest$statistic,
          G.pval = gtest$p.value,
          cucconi.stat = cucconi$C,
          cucconi.pval = cucconi$p.value
        )
      )
    }) %>%
    as_tibble
  
  res2 <-
    foo %>%
    ddply(.(scenario, run, niche_index), function(v){
      list2env(v, envir = environment())
      print(paste(unique(scenario), unique(run)))
      gtest <- try(DescTools::GTest(x = obs_N, p = exp_N/sum(exp_N)), silent = TRUE)
      return(
        tibble(
          G = 2 * sum(obs_N * log(obs_N / exp_N)),
          G.stat = gtest$statistic,
          G.pval = gtest$p.value
        )
      )
    }) %>%
    as_tibble %>%
    left_join(res)
  
  res_sum1 <-
    res2 %>%
    group_by(scenario) %>%
    summarize_at(c('niche_index','G.stat','cucconi.stat'),mean) %>%
    pivot_longer(cols = c(G.stat,cucconi.stat), names_to ='test')
  
  res_sum2 <-
    res2 %>% 
    group_by(scenario) %>% 
    summarize(
      niche_index = mean(niche_index), 
      cucconi.signif = mean(cucconi.pval < .05),
      G.signif = mean(G.pval < .05)
    ) %>%
    pivot_longer(cols = c('G.signif','cucconi.signif'), names_to = 'test') %>%
    mutate(Indistinguishability = 100 * (1 - value))
  
  chosen_ni <- sapply(seq(0,1,l=12), function(.) DescTools::Closest(foo$niche_index,.)[1])
  chosen_scenarios <- foo %>% filter(niche_index %in% chosen_ni) %>% pull(scenario) %>% unique()
  
  plot_test_stats <-
    res_sum1 %>%
    ggplot(aes(niche_index, value)) +
    geom_point() +
    facet_wrap(~test, scales = 'free') +
    scale_y_log10()
  
  plot_significance <-
    res_sum2 %>%
    ggplot(aes(niche_index, Indistinguishability)) +
    geom_point() +
    facet_wrap(~test)
  
  plot_n <-
    foo %>% 
    filter(niche_index %in% chosen_ni) %>%
    mutate(niche_index = round(niche_index,2)) %>%
    ggplot(aes(exp_N,obs_N)) + 
    geom_point() + 
    facet_wrap(~niche_index) + 
    geom_abline(slope = 1, intercept = 0, color = red) + 
    scale_x_log10() + 
    scale_y_log10()
  
  plot_niche_index <-
    bar %>%
    ggplot(aes(scenario,niche_index)) +
    geom_point()
  
  plot_cv <-
    bar %>%
    ggplot(aes(scenario,cv)) +
    geom_point()
  
  plot_niche_index_cv <-
    bar %>%
    ggplot(aes(niche_index, cv)) +
    geom_point()
  
  plot_var_n <-
    bar %>%
    ggplot(aes(niche_index, obs_N)) +
    geom_point() 
  
  plot_cor <-
    bar %>%
    ggplot(aes(niche_index,cor_n)) + 
    geom_point() 
  
  plot_rv <-
    bar %>%
    ggplot(aes(niche_index, RV)) +
    geom_point() 
  
  
  gridExtra::grid.arrange(
    plot_niche_index,
    plot_cv,
    plot_niche_index_cv,
    plot_var_n,
    plot_cor,
    plot_rv,
    nrow = 2
  )
    
  par(mfrow=c(3,4),mar=c(2,2,2,2))
  for(i in rev(chosen_scenarios)){
    image(get(load(paste0('scenario_',i,'_run_1.rdata')))$C)
  }
  
}

