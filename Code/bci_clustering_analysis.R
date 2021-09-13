library(tidyverse)
library(magrittr)
library(furrr)
library(parallel)
library(parallelDist)
library(readxl)
library(pcaMethods)


theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)


do.clustering.analysis = 1
do.recruitment.analysis = 0
do.nutrient.analysis = 0
do.trait.analysis = 0
do.paper.figures = 0

do.data = 0
do.plots = 1
fdp = 'bci'

filter = dplyr::filter

data_directory = 'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/'
suffix = '?raw=true'

read_datafile = 
  function(filename){
    readRDS(url(paste0(data_directory, fdp, filename, suffix)))
  }

census_data = read_datafile('_census_data.rds')
cluster_data = read_datafile('_clustering_analysis.rds')
kde_full = read_datafile('_inferred_soiltypes.rds')
soiltype_v_nutrients = read_datafile('_C5_soiltype_vs_nutrients.rds')
nutrients = read_datafile('_nutrient_data.rds')
if(fdp == 'bci') trait_data_raw = read_datafile('_trait_data.rds')


save_date = gsub('-', '', Sys.Date())
save_directory = paste0('~/SpatialNiche/Data/', save_date, '/')


L = ifelse(fdp == 'bci', 1000, 500)

if(do.clustering.analysis){
  
  if(do.data){
    mypc = (Sys.info()['sysname'] == 'Windows') || (Sys.info()['sysname'] =='Darwin')
    
    cores = if(mypc) 4 else detectCores() - 10
    plan(multisession, workers = cores)
    save_date = gsub('-', '', Sys.Date())
    save_directory = paste0('~/SpatialNiche/Data/', save_date, '/')
    
    source('https://github.com/rafaeldandrea/Spatial-niche/raw/Lap/Code/clustering_functions.R')
  
    if(fdp == 'bci'){
      Lx = 1000
      bci =
        readRDS(
          url(
            'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
          )
        )
      filename = paste0(save_directory, 'bci_clustering_analysis.rds')
      thecensus=1:7
    }
    
    if(fdp == 'lap'){
      Lx = 500
      bci  = read_all_laplanada()
      filename = paste0(save_directory, 'lap_clustering_analysis.rds')
      thecensus = 1:2
    }
    
    parameters = 
      expand_grid(
        thecensus = thecensus,
        algorithm = c('louvain', 'walktrap'),
        d_cutoff = seq(10, 30, by = 2),
        self_loops = FALSE,
        d_step = 1e-5,
        Lx = Lx,
        Ly = 500,
        autolinked_species_only = TRUE,
        weighted = TRUE,
        seed = 0:10
      )
    
    if(fdp %in% c('bci', 'lap')){
      parameters %<>%
        filter(algorithm == 'louvain', seed == 0)  
    }
      
    
    chunk_size = cores
    
    for(piece in seq(nrow(parameters) / chunk_size)){
      
      indices = chunk_size * (piece - 1) + seq(chunk_size)
      
      parms = parameters[indices, ]
      
      fdp_analyzed = 
        parms %>%
        future_pmap_dfr(
          .f = function(
            thecensus,
            algorithm,
            d_cutoff,
            self_loops,
            d_step,
            Lx,
            Ly,
            autolinked_species_only,
            weighted,
            seed
          ){
            
            dat =
              bci %>%
              filter(census == thecensus)
            
            if(seed > 0){
              dat %<>%
                mutate(sp = sample(sp))
            }
            
            result = 
              dat %>%
              adjacency_matrix(
                autolinked_species_only = autolinked_species_only, 
                d_cutoff = d_cutoff, 
                d_step = d_step, 
                Lx = Lx, 
                Ly = Ly
              ) %>%
              cluster_analysis(
                algorithm = algorithm,
                weighted = weighted,
                self_loops = self_loops
              ) %>%
              pluck('result') %>%
              mutate(
                census = thecensus,
                d_cutoff = d_cutoff,
                autolinked_species_only = autolinked_species_only,
                d_step = d_step,
                seed = seed
              ) %>%
              return()
          },
          .options = furrr_options(seed = TRUE)
        ) %>%
        rename(sp = name)
      
      dir.create(save_directory, showWarnings = FALSE)
      

      if(file.exists(filename)){
        fdp_analyzed = 
          readRDS(filename) %>%
          bind_rows(fdp_analyzed)
      }
      
      saveRDS(fdp_analyzed, file = filename)
    }
    
    cluster_data = readRDS(filename)
    
    adults = 
      census_data %>%
      filter(dbh >= 100) %>%
      select(census, treeID, sp, gx, gy, dbh)
    
    recruits = 
      adults %>%
      filter(census > 1) %>%
      group_by(treeID) %>%
      slice_min(census) %>%
      ungroup %>%
      mutate(recruit = TRUE)
    
    trees =
      adults %>% 
      left_join(
        recruits, 
        by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
      ) %>%
      replace_na(list(recruit = FALSE))
    
    combined_data =
      trees %>%
      inner_join(
        cluster_data, 
        by = c('census', 'sp')
      )
    
    reference =
      combined_data %>%
      inner_join(
        combined_data %>%
          select(census, d_cutoff, number_of_groups) %>%
          unique() %>%
          slice_max(number_of_groups, with_ties = FALSE)
      ) %>%
      count(group) %>%
      mutate(reference = rev(rank(n))) %>%
      full_join(
        combined_data %>%
          inner_join(
            combined_data %>%
              select(census, d_cutoff, number_of_groups) %>%
              unique() %>%
              slice_max(number_of_groups, with_ties = FALSE)
          )
      ) %>%
      select(treeID, reference)
    
    z_full = NULL
    for(thecensus in sort(unique(combined_data$census))){
      
      for(thed_cutoff in sort(unique(combined_data$d_cutoff))){
        
        z = NULL
        
        x =
          combined_data %>%
          filter(
            census == thecensus,
            d_cutoff == thed_cutoff
          ) %>%
          inner_join(
            reference,
            by = 'treeID'
          ) %>%
          count(
            census,
            d_cutoff,
            group,
            reference
          ) %>%
          mutate(group = factor(group)) %>%
          complete(group, nesting(census, d_cutoff, reference)) %>%
          replace_na(list(n = 0)) %>%
          arrange(desc(n))
        
        if(nrow(x) > 0){
          z %<>%
            bind_rows(
              tibble(
                census = thecensus,
                d_cutoff = thed_cutoff,
                group = x$group[1],
                reference = x$reference[1]
              )
            )
          
          for(i in 2:nrow(x)){
            
            g = x$group[i]
            r = x$reference[i]
            
            if(!g %in% z$group & !r %in% z$reference){
              z %<>%
                bind_rows(
                  tibble(
                    census = thecensus,
                    d_cutoff = thed_cutoff,
                    group = g,
                    reference = r
                  )
                )
              
              z_full %<>%
                bind_rows(z)
            }
          }
        }
      }
    }
    
    z_full %<>%
      unique()
    
    cluster_data_consistent =
      cluster_data %>%
      inner_join(z_full) %>%
      mutate(group = reference) %>%
      select(-reference)
    
    saveRDS(cluster_data_consistent, filename)
    
  }
  
  if(do.plots){
    data = readRDS('/Users/wangdz/Downloads/lap_clustering_analysis.rds')
    
    summary = 
      data %>%
      filter(seed == 0) %>%
      select(-c(name, group)) %>%
      unique() %>%
      group_by(
        algorithm, 
        weighted, 
        self_loops,
        d_cutoff, 
        autolinked_species_only, 
        d_step
      ) %>%
      summarize_at(
        c('modularity', 'number_of_groups'), 
        list(mean = mean, sd = sd, se = function(x) sd(x) / sqrt(length(x)))
      ) %>%
      ungroup()
    
    null_summary = 
      data %>%
      filter(seed > 0) %>%
      select(-c(name, group)) %>%
      unique() %>%
      group_by(
        algorithm, 
        weighted, 
        self_loops,
        d_cutoff, 
        autolinked_species_only, 
        d_step,
        census
      ) %>%
      summarize_at(
        c('modularity', 'number_of_groups'), 
        list(mean = mean, sd = sd, se = function(x) sd(x) / sqrt(length(x)))
      ) %>%
      ungroup()
    
    newdata = 
      data %>%
      filter(seed == 0) %>%
      select(-c(name, group)) %>%
      unique() %>%
      left_join(null_summary) %>%
      mutate(z = (modularity - modularity_mean) / modularity_sd)
    
    plot_bars = 
      summary %>%
      ggplot(aes(d_cutoff, number_of_groups_mean, fill = modularity_mean)) +
      geom_errorbar(
        aes(
          x = d_cutoff, 
          ymin = number_of_groups_mean - 2 * number_of_groups_se,
          ymax = number_of_groups_mean + 2 * number_of_groups_se
        )
      ) +
      geom_col(color = 'black') +
      facet_grid(weighted ~ algorithm, labeller = label_both) +
      scale_fill_gradientn(colors = terrain.colors(100))
    
    plot_colors =
      data %>%
      filter(seed == 0) %>%
      mutate(ngroups = factor(number_of_groups)) %>%
      ggplot(aes(d_cutoff, census, fill = ngroups)) +
      geom_tile(color = 'black') +
      facet_grid(weighted ~ algorithm, labeller = label_both)
    
    gridExtra::grid.arrange(plot_bars, plot_colors, ncol = 2)
    
    plot_effect_size_colors =
      newdata %>%
      ggplot(aes(d_cutoff, census, fill = z)) +
      geom_tile(color = 'black') +
      facet_grid(weighted ~ algorithm, labeller = label_both) +
      scale_fill_gradientn(colors = terrain.colors(100))
    
    plot_effect_size_lines =
      newdata %>%
      group_by(weighted, algorithm, d_cutoff) %>%
      summarize_at('z', list(mean = mean, sd = sd)) %>%
      ungroup() %>%
      ggplot(aes(d_cutoff, mean)) +
      geom_line() +
      facet_grid(weighted ~ algorithm, labeller = label_both)
    
    gridExtra::grid.arrange(plot_effect_size_colors, plot_effect_size_lines, ncol = 2)
    
  }
  
}

## Definition of theta := P(recruit | match) / P(recruit | !match)
if(do.recruitment.analysis){
  
  library(sparr) # for function bivariate.density() in KDE()
  filter = dplyr::filter
  
  ## Determine whether working on SeaWulf (SBU hpc) or personal computer
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 4
  plan(multisession, workers = cores)
  
  adults = 
    census_data %>%
    filter(dbh >= 100) %>%
    select(census, treeID, sp, gx, gy, dbh)
  
  recruits = 
    adults %>%
    filter(census > 1) %>%
    group_by(treeID) %>%
    slice_min(census) %>%
    ungroup %>%
    mutate(recruit = TRUE)
  
  trees =
    adults %>% 
    left_join(
      recruits, 
      by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
    ) %>%
    replace_na(list(recruit = FALSE))
  
  combined_data =
    trees %>%
    inner_join(
      cluster_data, 
      by = c('census', 'sp')
    )
  
  reference =
    combined_data %>%
    inner_join(
      combined_data %>%
        select(census, d_cutoff, number_of_groups) %>%
        unique() %>%
        slice_max(number_of_groups, with_ties = FALSE)
    ) %>%
    count(group) %>%
    mutate(reference = rev(rank(n))) %>%
    full_join(
      combined_data %>%
        inner_join(
          combined_data %>%
            select(census, d_cutoff, number_of_groups) %>%
            unique() %>%
            slice_max(number_of_groups, with_ties = FALSE)
        )
    ) %>%
    select(treeID, reference)

  z_full = NULL
  for(thecensus in sort(unique(combined_data$census))){

    for(thed_cutoff in sort(unique(combined_data$d_cutoff))){

      z = NULL

      x =
        combined_data %>%
        filter(
          census == thecensus,
          d_cutoff == thed_cutoff
        ) %>%
        inner_join(
          reference,
          by = 'treeID'
        ) %>%
        count(
          census,
          d_cutoff,
          group,
          reference
        ) %>%
        mutate(group = factor(group)) %>%
        complete(group, nesting(census, d_cutoff, reference)) %>%
        replace_na(list(n = 0)) %>%
        arrange(desc(n))

      if(nrow(x) > 0){
        z %<>%
          bind_rows(
            tibble(
              census = thecensus,
              d_cutoff = thed_cutoff,
              group = x$group[1],
              reference = x$reference[1]
            )
          )

        for(i in 2:nrow(x)){

          g = x$group[i]
          r = x$reference[i]

          if(!g %in% z$group & !r %in% z$reference){
            z %<>%
              bind_rows(
                tibble(
                  census = thecensus,
                  d_cutoff = thed_cutoff,
                  group = g,
                  reference = r
                )
              )

            z_full %<>%
              bind_rows(z)
          }
        }
      }
    }
  }

  z_full %<>%
    unique()

  combined_data_consistent =
    combined_data %>%
    left_join(z_full) %>%
    mutate(group = reference) %>%
    select(-reference)

  cluster_data =
    cluster_data %>%
    inner_join(z_full) %>%
    mutate(group = reference) %>%
    select(-reference)


  parms = 
    combined_data_consistent %>%
    select(
      Census = census, 
      Algorithm = algorithm,
      Seed = seed,
      D_cutoff = d_cutoff,
      Group = group
    ) %>%
    unique %>%
    filter(
      Algorithm == 'louvain',
      Seed == 0
    )
  
  
  KernelDensityEstimation = 
    function(gx, gy, Lx = L, Ly = 500, quadrat_length = 20, ...){
      
      evalpoints =
        expand_grid(
          x = quadrat_length / 2 + seq(0, Lx - quadrat_length, by = quadrat_length),
          y = quadrat_length / 2 + seq(0, Ly - quadrat_length, by = quadrat_length)
        ) %>%
        arrange(y, x)
      
      dens = 
        bivariate.density(
          ppp(gx, gy, xrange = c(0, Lx), yrange = c(0, Ly)), 
          h0 = quadrat_length, 
          xy = 
            list(
              x = evalpoints$x,
              y = evalpoints$y
            )
        )
      
      evalpoints %>%
        mutate(
          density = as.numeric(t(dens$z$v))
        ) %>%
        return
    }
  
  KDE = 
    function(Census, Algorithm, Seed, D_cutoff, Group, .data){
      df = 
        .data %>% 
        filter(
          census == Census,
          algorithm == Algorithm,
          seed == Seed,
          d_cutoff == D_cutoff,
          group == Group
        ) %>%
        select(gx, gy) %>%
        unique()
      
      result = 
        KernelDensityEstimation(
          gx = df$gx,
          gy = df$gy
        ) %>%
        mutate(
          census = Census,
          algorithm = Algorithm,
          seed = Seed,
          d_cutoff = D_cutoff,
          group = Group
        ) %>%
        return
    }
  
  
  kde_full = 
    parms %>%
    future_pmap_dfr(
      .f = KDE,
      .data = combined_data_consistent,
      .options = furrr_options(seed = TRUE)
    )
  
  soiltype = 
    kde_full %>%
    group_by(algorithm, seed, d_cutoff, census, x, y) %>%
    slice_max(density, with_ties = FALSE) %>%
    ungroup %>%
    rename(soiltype = group) %>%
    select(-density)
  
  kde_full %<>%
    left_join(
      soiltype,
      by = c('x', 'y', 'census', 'algorithm', 'seed', 'd_cutoff')
    ) %>% 
    mutate(fdp = fdp)
  
  # saveRDS(kde_full, file = kde_filename)
  
  plot_soiltypes = 
    kde_full %>% 
    mutate(soiltype = factor(soiltype)) %>% 
    ggplot(aes(x, y, fill = soiltype)) + 
    geom_tile() + 
    facet_grid(census ~ d_cutoff) + 
    theme(aspect.ratio = ifelse(fdp == 'bci', .5, 1))
  
  plot_soiltypes %>%
    show()
  
  
  rec_df = 
    recruits %>%
    mutate(
      x = seq(10, L - 10, 20)[cut(gx, seq(20, L, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
    ) %>%
    inner_join(
      kde_full %>%
        select(
          labelling_census = census, 
          x, 
          y, 
          algorithm, 
          seed, 
          d_cutoff, 
          soiltype, 
          fdp
        ) %>%
        unique(),
      by = c('x', 'y')
    ) %>%
    inner_join(
      cluster_data %>%
        rename(labelling_census = census)
    )
  
  Pm_df = 
    kde_full %>% 
    select(x, y, d_cutoff, soiltype, census) %>% 
    unique() %>% 
    count(census, d_cutoff, soiltype) %>% 
    mutate(Pm = n / (L * 500 / 20 ^ 2)) %>% 
    ungroup() %>%
    rename(labelling_census = census)
  
  res =
    rec_df %>%
    group_by(
      labelling_census,
      d_cutoff,
      group
    ) %>%
    summarize(
      recruits = n(),
      matches = sum(group == soiltype),
      .groups = 'drop'
    ) %>%
    left_join(
      Pm_df %>%
        rename(group = soiltype)
    ) %>%
    mutate(
      theta = ((1 - Pm) / Pm) / (recruits / matches - 1)
    )
  
  res_summary = 
    res %>%
    group_by(labelling_census, d_cutoff) %>%
    mutate(weight = recruits / sum(recruits)) %>%
    summarize(
      theta_mean = weighted.mean(theta, weight),
      theta_se = sqrt(sum(weight * (theta - theta_mean) ^ 2)),
      .groups = 'drop'
    )
  
  res_summary = 
    res %>%
    group_by(labelling_census, d_cutoff) %>%
    mutate(weight = recruits / sum(recruits)) %>%
    summarize(
      theta_mean = weighted.mean(theta, weight),
      theta_se = sqrt(sum(weight * (theta - theta_mean) ^ 2)),
      .groups = 'drop'
    ) %>%
    group_by(d_cutoff) %>%
    summarize(
      theta_mean = mean(theta_mean),
      theta_se = sqrt(sum(theta_se ^ 2)) / n(),
      .groups = 'drop'
    )
  
  plot_histograms = 
    res %>% 
    ggplot(aes(theta)) + 
    geom_histogram() + 
    facet_wrap(~ d_cutoff)
  
  plot_bars = 
    res_summary %>%
    ggplot(aes(d_cutoff, theta_mean)) +
    geom_hline(yintercept = 1, color = 'grey') +
    geom_errorbar(
      aes(
        x = d_cutoff, 
        ymin = theta_mean, 
        ymax = theta_mean + 2 * theta_se
      )
    ) +
    geom_col(fill = 'plum4') +
    labs(x = 'Distance cutoff', y = 'Theta estimate') +
    theme(aspect.ratio = 1) +
    scale_y_continuous(breaks = 0:6)
  
  plot_histograms %>%
    show()
  
  plot_bars %>%
    show
  
}

if(do.nutrient.analysis){
  
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 4
  plan(multisession, workers = cores)
  
  if(do.data){
    
    library(caret)
    library(C50)
    library(readxl)
    
    if (fdp=='bci'){
      nutrients %<>%
        pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
        group_by(nutrient) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
    }
    
    if (fdp == 'lap') {
      nutrients %<>%
        mutate(gx = gx + 10, gy = gy + 10) %>%  #move from left lower to center
        pivot_longer(-c(gx, gy), names_to = 'nutrient') %>%
        group_by(nutrient) %>%
        rename(x = gx) %>%
        rename(y = gy) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
    }
    
    
    nutrients_wide = 
      nutrients %>%
      select(-value) %>%
      pivot_wider(names_from = nutrient, values_from = standardized)
    
    
    kde = 
      kde_full %>%
      select(-c(group, density)) %>%
      unique()
    
    dtf =
      nutrients_wide %>%
      full_join(kde, by = c('x', 'y'))
    
    parms = 
      kde %>%
      select(
        Census = census, 
        Algorithm = algorithm, 
        Seed = seed, 
        D_cutoff = d_cutoff, 
        Fdp = fdp
      ) %>%
      unique()
    
    soiltype_v_nutrients = 
      parms %>%
      future_pmap_dfr(
        .options = furrr_options(seed = TRUE),
        .f = function(Census, Algorithm, Seed, D_cutoff, Fdp){ 
          data = 
            dtf %>% 
            filter(
              census == Census,
              algorithm == Algorithm,
              seed == Seed,
              d_cutoff == D_cutoff,
              fdp == Fdp
            ) %>%
            {if (fdp == 'bci') select(.,Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, `N(min)`, pH, soiltype)
              else .} %>% 
            {if (fdp == 'lap')select(.,Al, Ca, Cu, Fe, K, Mg, Mn, P, Zn, pH, soiltype)
              else .}%>%
            mutate(soiltype = factor(soiltype, levels = unique(soiltype)))
          
          set.seed(Seed)
          
          C5_model = 
            train(
              soiltype ~ ., 
              data = data, 
              method = 'C5.0',
              trControl = trainControl(method = 'repeatedcv', repeats = 10),
              metric = 'Kappa'
            )
          
          df = 
            C5_model$results %>%
            slice_max(Kappa, n = 1) %>%
            as_tibble %>%
            mutate(
              census = Census,
              algorithm = Algorithm,
              seed = Seed,
              d_cutoff = D_cutoff,
              fdp = Fdp
            ) %>%
            return()
        }
      )
    
    saveRDS(soiltype_v_nutrients, file = paste0(save_directory, fdp, '_C5_soiltype_vs_nutrients.rds'))
    
  }
  
  
  if(do.plots){
    
    soiltype_v_nutrients_summary = 
      soiltype_v_nutrients %>% 
      group_by(d_cutoff) %>% 
      summarize(
        mean = mean(Kappa), 
        se = sd(Kappa) / sqrt(n()), 
        .groups = 'drop'
      )
    
    plot_bars = 
      soiltype_v_nutrients_summary %>% 
      ggplot(aes(d_cutoff, mean)) + 
      geom_errorbar(
        aes(x = d_cutoff, ymin = mean - 2 * se, ymax = mean + 2 * se), 
        width = .5
      ) +
      geom_col(fill = 'plum4') + 
      labs(x = 'Distance cutoff', y = 'Cohen\'s Kappa') +
      coord_cartesian(ylim = c(0, 1)) +
      theme(aspect.ratio = 1)+
      ggtitle('Soil nutrients VS inferred soil types')
    
    
    if (fdp=='bci'){
      nutrients = 
        read_excel(
          nutrients_filename, 
          sheet = 2
        ) %>%
        pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
        group_by(nutrient) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
      
      kde_full = readRDS(kde_filename)
    }
    
    if (fdp == 'lap') {
      nutrients = 
        read_excel(
          nutrients_filename,
          sheet = 4
        ) %>%
        #move from left lower to center
        mutate(gx = gx + 10, gy = gy + 10) %>%
        pivot_longer(-c(gx, gy), names_to = 'nutrient') %>%
        group_by(nutrient) %>%
        rename(x = gx) %>%
        rename(y = gy) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
      kde_full = readRDS(kde_filename)
    }
    
    cor_analysis =
      nutrients %>%
      select(x, y, nutrient, standardized) %>%
      inner_join(
        kde_full %>%
          select(-soiltype) %>%
          rename(soiltype = group), 
        by = c('x', 'y')
      ) %>%
      group_by(
        census, 
        algorithm,
        seed,
        d_cutoff,
        fdp,
        nutrient,
        soiltype
      ) %>%
      summarize(
        cor = cor(density, standardized, use = 'complete.obs'), 
        p.value = cor.test(density, standardized)$p.value, 
        p.adj = p.adjust(p.value, method = 'hochberg'),
        significant = (p.adj < .05),
        cor_sig = ifelse(significant == TRUE, cor, NA),
        .groups = 'drop'
      ) 
    
    ranked_soiltypes = 
      cor_analysis %>%
      group_by(
        algorithm,
        seed,
        fdp,
        d_cutoff,
        census,
        soiltype
      ) %>%
      summarize(
        order = -mean(cor),
        .groups = 'drop'
      ) %>%
      group_by(
        algorithm,
        seed,
        fdp,
        d_cutoff,
        census
      ) %>%
      mutate(ranked_soiltype = factor(rank(order))) %>%
      ungroup()
    
    plot_nutrient_soiltype_correlations =
      cor_analysis %>%
      left_join(
        ranked_soiltypes,
        by = c("census", "algorithm", "seed", "d_cutoff", "fdp", "soiltype")
      ) %>% 
      filter(d_cutoff == 20) %>%
      ggplot(aes(nutrient, cor_sig, group = ranked_soiltype, fill = ranked_soiltype)) + 
      geom_col(position = 'dodge') + 
      facet_grid(census ~ ranked_soiltype, labeller = label_both) +
      labs(fill = 'group', y = 'Pearson correlation coefficient') +
      ggtitle('Distance cutoff = 20 m')
    
    
    plot_kde = 
      kde_full %>% 
      left_join(
        ranked_soiltypes,
        by = c("census", "algorithm", "seed", "d_cutoff", "fdp", "soiltype")
      ) %>%
      filter(d_cutoff %in% c(6, 8, 10, 16, 20, 30)) %>% 
      ggplot(aes(x, y, fill = ranked_soiltype)) + 
      geom_tile() + 
      facet_grid(d_cutoff ~ census, labeller = label_both) + 
      theme(aspect.ratio = .5)
    
    plot_bars %>%
      show
    
    plot_nutrient_soiltype_correlations %>%
      show
    
    plot_kde %>%
      show
    
  }
  
}  

if(do.trait.analysis){
  
  
  size_traits = 
    c(
      "DBH_AVG", 
      "HEIGHT_AVG", 
      "DIAM_AVG"
    )
  
  leaf_traits = 
    c(
      "LEAFAREA_AVD",
      "LEAFTHCK_AVD",
      "LMADISC_AVD",
      "LMALEAF_AVD",
      "LMALAM_AVD",
      "LDDISC_AVD",
      "LDMC_AVD",
      "LEAFAREA_AVI",
      "LEAFTHCK_AVI",
      "LMADISC_AVI",
      "LMALEAF_AVI",
      "LMALAM_AVI" ,
      "LDDISC_AVI" ,
      "LDMC_AVI"  ,
      "AVG_LAMTUF" ,
      "AVG_VEINTUF"
    )
  
  seed_traits = 
    c(
      "FRUIT_FRSH",
      "FRUIT_DRY" ,
      "DSPR_FRESH",
      "DSPR_DRY",
      "SEED_FRESH",
      "SEED_DRY"
    )
  
  wood_traits = 
    c(
      "SG60C_AVG",
      "SG100C_AVG" 
    )
  
  vital_traits =
    c(
      "RGR_10",      
      "RGR_50",
      "RGR_100",
      "MORT_10",
      "MORT_100"
    )
  
  traitlist = 
    list(
      vital = vital_traits, 
      leaf = leaf_traits, 
      seed = seed_traits, 
      wood = wood_traits, 
      size = size_traits
    )
  
  
  foo = 
    trait_data_raw %>%
    mutate(sp = tolower(`SP$`)) %>%
    pivot_longer(-c(1:6, sp), names_to = 'trait') %>%
    filter(!is.na(value)) %>%
    filter(!str_detect(trait, '_N')) %>%
    filter(!str_detect(trait, 'N_')) %>%
    filter(!str_detect(trait, '_SE')) %>%
    filter(!str_detect(trait, 'SEM_')) %>%
    filter(value > 0) %>%
    mutate(logged_value = log(value))
  
  normality = 
    foo %>%
    group_by(trait) %>%
    summarize(
      normal = shapiro.test(value)$p.value > .05,
      lognormal = shapiro.test(logged_value)$p.value > .05,
      .groups = 'drop'
    )
  
  trait_data =
    foo %>%
    left_join(normality, by = 'trait') %>%
    mutate(standardized = ifelse(normal, value, logged_value)) %>%
    group_by(trait) %>%
    mutate(standardized = scale(standardized)[, 1]) %>%
    ungroup %>%
    select(sp, trait, standardized)
  
  pca_data = 
    seq_along(traitlist) %>%
    map_dfr(
      .f = function(i){
        subdat = 
          trait_data %>%
          filter(trait %in% traitlist[[i]]) %>%
          select(sp, trait, standardized) %>%
          pivot_wider(names_from = trait, values_from = standardized) 
        
        subpca = 
          subdat %>%
          pca(method = 'ppca', scale = 'none', center = FALSE)
        
        tibble(
          sp = subdat$sp,
          trait_type = names(traitlist[i]),
          pc1 = subpca@scores[, 1],
          pc2 = subpca@scores[, 2]
        ) %>%
          return
      }
    ) %>%
    left_join(
      trait_data %>%
        select(sp) %>%
        unique,
      by = 'sp'
    )
  
  trait_tibble = 
    tibble(
      trait_type = rep(names(traitlist), lengths(traitlist)), 
      trait = unlist(traitlist)
    )
  
  pc1_sign = 
    trait_data %>%
    inner_join(trait_tibble) %>%
    inner_join(pca_data) %>%
    group_by(trait, trait_type) %>%
    summarize(
      cor_sign = sign(cor(pc1, standardized)),
      .groups = 'drop'
    ) %>% 
    group_by(trait_type) %>% 
    summarize(
      cor_sign = round(mean(cor_sign)), 
      .groups = 'drop'
    )
  
  pca_data %<>%
    inner_join(pc1_sign) %>%
    mutate(pc1 = pc1 * cor_sign)
  
  ## Read cluster data and name groups by correlation with soil nutrients
  # cluster_data = 
  #   readRDS('~/SpatialNiche/Data/20210722/bci_clustering_analysis.rds') %>%
  #   rename(sp = name) %>%
  #   filter(
  #     algorithm == 'louvain',
  #     weighted == TRUE,
  #     seed == 0
  #   ) %>%
  #   mutate(
  #     group = 
  #       factor(
  #         group, 
  #         levels = sort(unique(group))
  #       )
  #   )
  
 nutrients %<>%
   pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
   group_by(nutrient) %>%
   mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
   ungroup
 
 all_data = 
   census_data %>%
   filter(dbh >= 100) %>%
   inner_join(cluster_data) %>%
   mutate(
     x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
     y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
   ) %>%
   inner_join(
     kde_full %>%
       select(x, y, census, d_cutoff, soiltype, fdp) %>%
       unique()
   ) %>%
   filter(census == 7, d_cutoff == 20) %>%
   inner_join(nutrients) %>%
   select(census, d_cutoff, sp, gx, gy, x, y, group, soiltype, nutrient, standardized) %>%
   inner_join(pca_data)
 
 correlations = 
   all_data |> 
   group_by(trait_type, nutrient) |> 
   summarize(
     cor = cor(standardized, pc1), 
     pval = cor.test(standardized, pc1)$p.value,
     .groups = 'drop'
    ) 
 
 cor_analysis =
   nutrients %>%
   select(x, y, nutrient, standardized) %>%
   inner_join(
     kde_full %>%
       select(-soiltype) %>%
       rename(soiltype = group), 
     by = c('x', 'y')
   ) %>%
   group_by(
     census, 
     algorithm,
     seed,
     d_cutoff,
     fdp,
     nutrient,
     soiltype
   ) %>%
   summarize(
     cor = cor(density, standardized, use = 'complete.obs'), 
     p.value = cor.test(density, standardized)$p.value, 
     p.adj = p.adjust(p.value, method = 'hochberg'),
     significant = (p.adj < .05),
     cor_sig = ifelse(significant == TRUE, cor, NA),
     .groups = 'drop'
   ) 
 
 mean_correlation = 
   cor_analysis |> 
   filter(census == 7, d_cutoff == 20, !nutrient %in% c('Al', 'pH')) |> 
   group_by(soiltype) |> 
   summarize(mean_cor = mean(cor), .groups = 'drop') |>
   mutate(
     group_ID = 
       factor(
         soiltype, 
         levels = soiltype[order(mean_cor, decreasing = TRUE)]
       )
   )
 
 all_data %<>%
   inner_join(
     mean_correlation |>
       select(group = soiltype, group_ID)
   )
 
 cor_analysis %<>%
   inner_join(
     mean_correlation |>
       select(soiltype, group_ID)
   )  
 
 
 plot_correlations = 
   correlations |> 
   filter(pval <= .05) |>
   ggplot(aes(trait_type, nutrient, fill = cor)) + 
   geom_tile() + 
   scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue')
 
 plot_correlations |>
   show()
 
 plot_group_densities = 
   all_data |> 
   filter(trait_type == 'vital', nutrient == 'Fe') |> 
   ggplot(aes(gx, gy)) + 
   geom_density_2d_filled() + 
   theme(aspect.ratio = .5) + 
   facet_wrap(~ group_ID)
 
 plot_trait_violins_trees = 
   all_data |> 
   filter(nutrient == 'Fe') |> 
   ggplot(aes(group_ID, pc1, fill = group_ID)) + 
   geom_violin(draw_quantiles = .5) + 
   facet_wrap(~ trait_type)
 
 plot_trait_violins_species = 
   all_data |> 
   filter(nutrient == 'Fe') |> 
   select(sp, group_ID, pc1, trait_type) |>
   unique() |>
   ggplot(aes(group_ID, pc1, fill = group_ID)) + 
   geom_violin(draw_quantiles = .5) + 
   facet_wrap(~ trait_type)
   
  plot_group_nutrient_correlations = 
    cor_analysis |> 
    filter(census == 7, d_cutoff == 20) |> 
    ggplot(aes(nutrient, cor, fill = group_ID)) + 
    geom_col() + 
    facet_wrap(~ group_ID)
  
  
  ranked_groups = 
    cor_analysis %>%
    group_by(
      algorithm,
      seed,
      fdp,
      d_cutoff,
      census,
      soiltype
    ) %>%
    summarize(
      order = -mean(cor),
      .groups = 'drop'
    ) %>%
    group_by(
      algorithm,
      seed,
      fdp,
      d_cutoff,
      census
    ) %>%
    mutate(
      ranked_group = factor(rank(order))
    ) %>%
    ungroup() %>%
    rename(group = soiltype)
  
  data = 
    cluster_data %>%
    mutate(group = factor(group)) %>%
    # inner_join(ranked_groups) %>%
    inner_join(pca_data) %>%
    rename(ranked_group = group)
  
  dcut = 20
  plot_violins = 
    data %>%
    filter(d_cutoff == dcut) %>%
    ggplot(aes(ranked_group, pc1, group = ranked_group, fill = ranked_group)) +
    geom_violin() + 
    facet_grid(trait_type ~ census, labeller = label_both) +
    theme(legend.position = 'none') +
    ggtitle(paste0('Distance cutoff = ', dcut))
  
  plot_violins %>%
    show
  
  PairwiseWilcox = function(pc1, ranked_group, ...){
    
    test = pairwise.wilcox.test(pc1, ranked_group)
    
    p.value = test$p.value
    
    pairs = 
      expand_grid(
        rn = as.numeric(rownames(p.value)), 
        cn = as.numeric(colnames(p.value))
      ) %>% 
      filter(rn > cn) %>% 
      arrange(cn, rn) %>%
      mutate(name = paste0(rn, cn)) %>% 
      pull(name)
    
    tibble(
      pair = pairs, 
      p.value = as.numeric(p.value[!is.na(p.value)])
    ) %>%
      return()
  }
  
  parms = 
    data %>%
    select(
      census,
      d_cutoff,
      algorithm,
      weighted,
      self_loops,
      number_of_groups,
      autolinked_species_only,
      d_step,
      trait_type
    ) %>%
    unique()
  
  pairwise = NULL
  for(i in seq(nrow(parms))){
    suppressMessages(
      foo <- 
        with(data %>% inner_join(parms[i, ]), PairwiseWilcox(pc1, ranked_group)) %>%
        bind_cols(parms[i, ])
    )
    
    pairwise = rbind(pairwise, foo)
  }
  
  pairwise %<>%
    mutate(significant = p.value < .05)
  
  plot_intergroup_significance = 
    pairwise %>% 
    group_by(d_cutoff) %>% 
    summarize(
      sig = 100 * mean(significant), 
      .groups = 'drop'
    ) %>% 
    ggplot(aes(d_cutoff, sig)) + 
    geom_point() + 
    geom_line() +
    labs(
      x = 'Distance cutoff', 
      y = 'Percentage of significantly different groups'
    ) +
    ggtitle('t-test: traits of inferred niches') +
    theme(aspect.ratio = 1) +
    geom_hline(yintercept = 5, col = 'red')
  
  plot_intergroup_significance %>%
    show
  
  if(!exists('combined_data')){
    census_data = readRDS(census_filename)
    
    cluster_data = readRDS(cluster_filename) %>%
      rename(sp = name) %>%
      filter(
        algorithm == 'louvain',
        weighted == TRUE,
        seed == 0
      ) %>%
      mutate(
        group = 
          factor(
            group, 
            levels = sort(unique(group))
          )
      )
    
    adults = 
      census_data %>%
      filter(dbh >= 100) %>%
      select(census, treeID, sp, gx, gy, dbh)
    
    recruits = 
      adults %>%
      filter(census > 1) %>%
      group_by(treeID) %>%
      slice_min(census) %>%
      ungroup %>%
      mutate(recruit = TRUE)
    
    trees =
      adults %>% 
      left_join(
        recruits, 
        by = c('census', 'treeID', 'sp', 'gx', 'gy', 'dbh')
      ) %>%
      replace_na(list(recruit = FALSE))
    
    combined_data =
      trees %>%
      inner_join(
        cluster_data, 
        by = c('census', 'sp')
      )
    
    parms = 
      combined_data %>%
      select(
        Census = census, 
        Algorithm = algorithm,
        Seed = seed,
        D_cutoff = d_cutoff,
        Group = group
      ) %>%
      unique %>%
      filter(
        Algorithm == 'louvain',
        Seed == 0
      )
  }
  
  
  y = 
    combined_data %>%
    mutate(
      x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
    ) %>%
    inner_join(data) %>%
    inner_join(nutrients)
  
  
  
}  

if(do.C5.trait.analysis){
  
  dtf =
    all_data %>% 
    select(sp, trait_type, pc1, group_ID) %>% 
    unique %>% 
    pivot_wider(names_from = trait_type, values_from = pc1) %>%
    arrange(sp)
  
  abundances =
    all_data %>% 
    select(sp, gx, gy) %>% 
    unique() %>%
    count(sp)
  
  C5_model_trait = 
    train(
      group_ID ~ ., 
      data = dtf, 
      method = 'C5.0',
      na.action = na.pass,
      trControl = trainControl(method = 'repeatedcv', repeats = 10),
      metric = 'Kappa'
    )
  
}

## ============ OLD CODE ============

if(FALSE){
  fn = function(sub){
    tbl = 
      expand_grid(
        tibble(
          sp1 = sub$name, 
          gr1 = sub$group
        ),
        tibble(
          sp2 = sub$name,
          gr2 = sub$group
        )
      ) %>%
      mutate(same = 1 * (gr1 == gr2)) %>%
      return()
  }
  
  parms = 
    expand_grid(
      algorithm = c('louvain', 'walktrap'),
      weighted = c(TRUE, FALSE),
      d_cutoff = seq(6, 30, by  = 2),
      census = 1:7
    )
  
  res =
    seq(nrow(parms)) %>%
    future_map_dfr(
      .f = function(index){
        sub = 
          parms[index, ] %>%
          inner_join(
            data %>% filter(seed == 0), 
            by = c('algorithm', 'weighted', 'd_cutoff', 'census'))
        
        tbl = 
          expand_grid(
            tibble(
              sp1 = sub$name, 
              gr1 = sub$group
            ),
            tibble(
              sp2 = sub$name,
              gr2 = sub$group
            )
          ) %>%
          bind_cols(parms[index,]) %>%
          mutate(same = 1 * (gr1 == gr2))
      },
      .options = furrr_options(seed = TRUE)
    )
  
  summary_res = 
    res %>%
    mutate(
      sp1 = factor(sp1, levels = sort(unique(data$name)), ordered = TRUE),
      sp2 = factor(sp2, levels = sort(unique(data$name)), ordered = TRUE)
    ) %>%
    filter(sp1 > sp2) %>%
    group_by(
      algorithm,
      d_cutoff,
      weighted, 
      sp1,
      sp2
    ) %>%
    summarize(
      sum = sum(same), 
      .groups = 'drop'
    )
  
  plot_same = 
    summary_res %>% 
    group_by(algorithm, d_cutoff, weighted) %>% 
    count(sum) %>% 
    ungroup() %>%
    mutate(sum = factor(sum)) %>%
    ggplot(aes(sum, n)) + 
    geom_col() +
    facet_grid(algorithm + weighted ~ d_cutoff)
  
  
}

## I ended up not using this function
bayes = 
  function(
    matches, 
    recruits, 
    nsoiltypes, 
    alpha, 
    beta, 
    theta_min = 0, 
    theta_max = 5, 
    dtheta = .01
  ){
    Theta = seq(theta_min, theta_max, by = dtheta) 
    Phi = Theta / (1 + Theta)
    prior = dbeta(Phi, shape1 = alpha, shape2 = beta)
    likelihood = 
      sapply(
        Theta, 
        FUN = 
          function(theta)
            exp(
              sum(
                log(
                  dbinom(
                    matches, 
                    size = recruits, 
                    prob = theta / (theta + nsoiltypes - 1)
                  )
                ),
                na.rm = TRUE
              )
            )
      )
    
    df = 
      tibble(
        phi = Phi,
        theta = Theta,
        prior = prior,
        posterior_theta = prior * likelihood / sum(prior * likelihood * dtheta),
        posterior_phi = posterior_theta * (1 + Theta) ^ 2
      )
    
    return(df)
  }
