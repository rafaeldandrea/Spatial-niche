library(tidyverse)
library(magrittr)
library(RandomFields)
library(furrr)
library(parallel)

do.data = 0
do.plots = 1

if(do.data){
  plan(multisession, workers = detectCores() - 1)
  
  source('~/SpatialNiche/R Scripts/cluster_analysis.R')
  
  generate_landscape =
    function(
      Lx,
      Ly,
      quadrat_length,
      rangepar,
      sillpar,
      nuggetpar,
      num_soiltypes,
      seed
    ) {
      stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
      RFoptions(seed = seed)
      
      stress = 
        RFsimulate(
          RMgauss(
            scale = rangepar + 1e-16,
            var = 2 * sillpar * (1 - nuggetpar)
          ) + 
            RMtrend(mean = 0) + 
            RMnugget(var = 2 * sillpar * nuggetpar),
          x = seq(Lx / quadrat_length),
          y = seq(Ly / quadrat_length)
        )@data$variable1
      
      landscape = 
        expand_grid(
          y = seq(Ly / quadrat_length),
          x = seq(Lx / quadrat_length)
        ) %>%
        mutate(
          soiltype = 
            factor(
              cut(
                stress, 
                num_soiltypes, 
                labels = FALSE
              )
            )
        )
      
      return(landscape)
    }
  
  recruitment = 
    function(
      soiltype, 
      species_group,
      species_abundance,
      theta
    ){
      prob_mat = 
        outer(
          soiltype, 
          species_group, 
          function(s, g) ifelse(s == g, theta, 1)
        )
      
      species = 
        apply(
          prob_mat, 
          1, 
          function(probs){
            sample(
              seq_along(probs), 
              size = 1, 
              prob = probs * species_abundance
            )
          } 
        )
      
      return(species)
    }
  
  build_recruits = 
    function(
      coords, 
      landscape, 
      species_group,
      species_abundance,
      theta
    ){
      df =
        coords %>%
        mutate(
          x = 
            cut(
              gx, 
              breaks = seq(0, Lx, quadrat_length), 
              labels = FALSE
            ),
          y = 
            cut(
              gy, 
              breaks = seq(0, Ly, quadrat_length), 
              labels = FALSE
            )
        ) %>%
        left_join(
          landscape, 
          by = c('x', 'y')
        ) %>%
        mutate(
          species = 
            soiltype %>%
            recruitment(
              species_group = species_group,
              species_abundance = species_abundance,
              theta = theta
            )
        ) %>%
        return()
    }
  
  ## average intercensus mortality and births are each 
  ## typically ~= 10% of the number of living trees
  bci_all = readRDS('~/BCI/all_data.rds')
  
  censuses = 
    expand_grid(
      census1 = 1:7,
      census2 = 1:7
    ) %>%
    filter(census2 == census1 + 1)
  
  ndeaths = 
    censuses %>%
    pmap_dbl(
      .f = function(census1, census2){
        a = 
          bci_all %>% 
          filter(
            census == census1, 
            status == 'A',
            dbh >= 100
          )
        
        d = 
          bci_all %>% 
          filter(
            census == census2, 
            status == 'D'
          )
        
        length(intersect(a$treeID, d$treeID)) %>%
          return
      }
    )
  
  nbirths = 
    censuses %>%
    pmap_dbl(
      .f = function(census1, census2){
        a1 = 
          bci_all %>% 
          filter(
            census == census1, 
            status == 'A',
            dbh < 100
          )
        
        a2 = 
          bci_all %>% 
          filter(
            census == census2, 
            status == 'A',
            dbh >= 100
          )
        
        length(intersect(a1$treeID, a2$treeID)) %>%
          return
      }
    )
  
  mean_births = mean(c(nbirths, ndeaths))
  sd_births = sd(nbirths)
  
  mean_deaths = mean(c(nbirths, ndeaths))
  sd_deaths = sd(ndeaths)
  
  
  bci =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
      )
    )
  
  results = readRDS('~/SpatialNiche/Data/20210712/all_bci_censuses.rds')
  
  bci %<>% 
    filter(sp %in% intersect(sp, results$name))
  
  abuns = 
    bci %>%
    count(census, sp)
  
  number_of_species = 
    abuns %>% 
    filter(n >= 40) %>%
    pull(sp) %>% 
    unique() %>% 
    length()
  
  census_turnover = 
    bci %>%
    count(census) %>%
    mutate(
      mean = mean(n),
      turnover = n - mean(n)
    )
  
  average_community_size = unique(census_turnover$mean)
  
  parameters =
    expand_grid(
      average_community_size = average_community_size,
      nspecies = number_of_species,
      nsoiltypes = 4,
      ncensuses = 100,
      Lx = 1000,
      Ly = 500,
      quadrat_length = 20,
      rangepar = 20, 
      sillpar = 1, 
      nuggetpar = .001, 
      seed = 0, 
      theta = 2
    )
  
  parms = parameters[1, ]
  
  list2env(as.list(parms), envir = environment())
  
  census = 0
  
  species_list = seq(nspecies)
  
  species_group = 
    rep(seq(nsoiltypes), each = nspecies / nsoiltypes)
  
  delta = nspecies - length(species_group)
  if(delta > 0){
    species_group = c(rep(1, delta), species_group)
  }
  
  landscape = 
    generate_landscape(
      Lx = Lx, 
      Ly = Ly, 
      quadrat_length = quadrat_length, 
      rangepar = rangepar, 
      sillpar = sillpar, 
      nuggetpar = nuggetpar, 
      seed = seed, 
      num_soiltypes = nsoiltypes
    )
  
  community = 
    tibble(
      gx = runif(average_community_size, min = 0, max = Lx),
      gy = runif(average_community_size, min = 0, max = Ly)
    ) %>%
    build_recruits(
      landscape = landscape,
      species_group = species_group,
      species_abundance = 1,
      theta = theta
    )
  
  data = 
    tibble(
      census = census
    ) %>%
    bind_cols(community)
  
  census_deaths = 
    rnorm(
      ncensuses, 
      mean = mean_deaths,
      sd = sd_deaths
    )
  
  census_births = 
    rnorm(
      ncensuses, 
      mean = mean_births,
      sd = sd_births
    )
  
  while(census < ncensuses){
    census = census + 1
    ndeaths = census_deaths[census]
    nbirths = census_births[census]
    
    abundances = 
      community %>%
      count(species) %>%
      right_join(
        tibble(
          species = seq(nspecies)
        ),
        by = 'species'
      ) %>%
      replace_na(list(n = 0))
    
    births = 
      tibble(
        gx = runif(nbirths, min = 0, max = Lx),
        gy = runif(nbirths, min = 0, max = Ly)
      ) %>%
      build_recruits(
        landscape = landscape,
        species_group = species_group,
        species_abundance = abundances$n,
        theta = theta
      )
    
    community %<>%
      sample_n(size = nrow(community) - ndeaths) %>%
      bind_rows(births)
    
    if(census %% 10 == 0){
      data %<>%
        bind_rows(
          tibble(
            census = census
          ) %>%
            bind_cols(community)  
        )   
    }
    
    
  }
  
  result = 
    data %>%
    pull(census) %>%
    unique() %>%
    future_map_dfr(
      .options = furrr_options(seed = TRUE),
      .f = function(cen){
        dat = 
          data %>%
          filter(census == cen) %>%
          rename(sp = species)
        
        A = 
          adjacency_matrix(
            dat, 
            autolinked_species_only = TRUE, 
            d_cutoff = 10, 
            d_step = 1e-5, 
            Lx = 1000, 
            Ly = 500
          )
        
        clusters = 
          cluster_analysis(
            A = A, 
            algorithm = 'louvain', 
            weighted = FALSE, 
            self_loops = FALSE
          )
        
        clusters$result %>%
          mutate(census = cen) %>%
          return
      }
    )
  
  
}

if(do.plots){
  
  theme_set(theme_bw())
  theme_update(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  plot_ngroups = 
    result %>%
    mutate(number_of_groups = factor(number_of_groups)) %>%
    group_by(nsoiltypes, number_of_groups) %>%
    count() %>%
    ungroup() %>%
    mutate(p = n / sum(n) * 100) %>%
    ggplot(aes(number_of_groups, p, fill = number_of_groups)) +
    geom_col(color = 'black') +
    theme(legend.position = 'none') +
    facet_wrap(~nsoiltypes, labeller = label_both) +
    labs(x = 'number of groups', y = 'probability (%)')
  
  plot_modularity =
    result %>% 
    mutate(ngroups = factor(number_of_groups)) %>%
    ggplot(aes(ngroups, modularity, fill = ngroups)) + 
    geom_boxplot() +
    facet_wrap(~nsoiltypes, labeller = label_both) +
    theme(legend.position = 'none') +
    labs(x = 'number of groups')
  
  gridExtra::grid.arrange(plot_ngroups, plot_modularity, nrow = 1 )
  
}
