library(tidyverse)
library(magrittr)
library(furrr)
library(parallel)
library(parallelDist) ## function parDist(X) calculates distances between rows of matrix X
library(pdist) ## for function pdist(X, Y), calculates distances between rows of matrices X and Y
library(readxl)
# library(pcaMethods)
library(sparr) # for function bivariate.density() in KDE()


theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)


do.clustering.analysis = 0
do.kde.analysis = 0
do.cfa.analysis = 1
do.C5.analysis = 0
do.recruitment.analysis = 0
do.nutrient.analysis = 0
do.trait.analysis = 0
do.paper.figures = 0
do.network.analysis = 0


do.data = 1
do.plots = 0
fdp = 'bci'

filter = dplyr::filter

data_directory = 'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/'
suffix = '?raw=true'

read_datafile = 
  function(filename, manual = FALSE){
    if(manual){
      readRDS(filename)
    }else{
      readRDS(url(paste0(data_directory, fdp, filename, suffix)))  
    } 
  }

# cluster_data = readRDS('~/SpatialNiche/Data/20211022/bci_clustering_analysis_consistent.rds')
# kde_full = readRDS('~/SpatialNiche/Data/20211004/bci_inferred_soiltypes.rds')

census_data = read_datafile('_census_data.rds')
# cluster_data = read_datafile('_clustering_analysis.rds')
# kde_full = read_datafile('_inferred_soiltypes.rds')
soiltype_v_nutrients = read_datafile('_C5_soiltype_vs_nutrients.rds')
nutrients = read_datafile('_nutrient_data.rds')
if(fdp == 'bci') trait_data_raw = read_datafile('_trait_data.rds')


save_date = gsub('-', '', Sys.Date())
save_directory = paste0('~/SpatialNiche/Data/', save_date, '/')


L = ifelse(fdp == 'bci', 1000, 500)

if(do.clustering.analysis){
  
  if(do.data){
    seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
    office = as.logical(Sys.info()['user'] == 'Rafael D\'Andrea')
    laptop = as.logical(Sys.info()['user'] == 'rdand')
    
    cores = max(4, if(seawulf) detectCores() - 2)
    plan(multisession, workers = cores)
    
    save_date = gsub('-', '', Sys.Date())
    save_directory = paste0('~/SpatialNiche/Data/', save_date, '/')
    
    source('https://github.com/rafaeldandrea/Spatial-niche/raw/main/Code/clustering_functions.R')
    
    if(fdp == 'bci'){
      Lx = 1000
      filename = paste0(save_directory, 'bci_clustering_analysis.rds')
      
      bci = 
        readRDS(
          url(
            'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bcifullcensus.rds?raw=true'
          )
        )
      
      baldeck_cutoff = 
        bci |>
        group_by(sp) |>
        summarize(
          baldeck = quantile(dbh, .56), 
          .groups ='drop'
        )
      
      dat = 
        bci |>
        inner_join(baldeck_cutoff) |>
        filter(dbh > baldeck) |>
        mutate(fdp = 'bci')
      
    }
    
    if(fdp == 'lap'){
      Lx = 500
      bci  = read_all_laplanada()
      filename = paste0(save_directory, 'lap_clustering_analysis.rds')
      
      baldeck_cutoff = 
        bci |>
        group_by(sp) |>
        summarize(
          baldeck = quantile(dbh, .56), 
          .groups ='drop'
        )
      
      dat = 
        bci |>
        inner_join(baldeck_cutoff) |>
        filter(dbh > baldeck) |>
        mutate(fdp = 'lap')
    }
    
    parameters = 
      expand_grid(
        thecensus = 
          dat |>
          pull(census) |>
          unique(),
        algorithm = 'louvain',
        d_cutoff = c(10, 20, 30),
        Lx = Lx,
        Ly = 500,
        weighted = TRUE,
        seed = 0:100
      ) |>
      filter(thecensus == 8 | d_cutoff == 20) |>
      filter(thecensus == 8 | seed == 0)
    
    fdp_analyzed = 
      parameters %>%
      future_pmap_dfr(
        .f = function(
          thecensus,
          algorithm,
          d_cutoff,
          Lx,
          Ly,
          weighted,
          seed
        ){
          
          data =
            dat %>%
            filter(census == thecensus)
          
          if(seed > 0){
            data %<>%
              mutate(sp = sample(sp))
          }
          
          result = 
            data %>%
            adjacency_matrix(
              d_cutoff = d_cutoff, 
              Lx = Lx, 
              Ly = Ly
            ) %>%
            cluster_analysis(
              algorithm = algorithm,
              weighted = weighted
            ) %>%
            pluck('result') %>%
            mutate(
              census = thecensus,
              d_cutoff = d_cutoff,
              seed = seed
            ) %>%
            return()
        },
        .options = furrr_options(seed = TRUE)
      ) %>%
      rename(sp = name)
    
    if(seawulf) dir.create(save_directory, showWarnings = FALSE)
    
    
    # if(file.exists(filename)){
    #   fdp_analyzed = 
    #     readRDS(filename) %>%
    #     bind_rows(fdp_analyzed)
    # }
    
    if(seawulf) saveRDS(fdp_analyzed, file = filename)
    
    cluster_data = readRDS(filename)
    filename_consistent = paste0(save_directory, 'bci_clustering_analysis_consistent.rds')
    
    x = 
      cluster_data %>% 
      filter(seed == 0) %>% 
      select(census, d_cutoff, sp, group, number_of_groups)
    
    x0 = 
      x %>%
      filter(census == 8, d_cutoff == 20)
    
    res = NULL
    for(rg in 1:4){
      foo = 
        x %>%
        group_by(census, d_cutoff, group) %>%
        summarize(
          int = length(intersect(sp, x0 %>% filter(group == rg) %>% pull(sp))),
          .groups = 'drop'
        ) %>%
        group_by(census, d_cutoff) %>%
        slice_max(int, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(reference_group = rg) %>%
        select(-int)
      
      res = 
        res %>%
        bind_rows(foo)
      
    }
    
    consistent_cluster_data = 
      cluster_data %>% 
      full_join(res) %>% 
      replace_na(list(reference_group = 5)) %>%
      select(-group) %>%
      rename(group = reference_group)
    
    saveRDS(consistent_cluster_data, filename_consistent)
    
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

if(do.network.analysis){
  bci =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
      )
    )
  
  source('https://github.com/rafaeldandrea/Spatial-niche/raw/main/Code/clustering_functions.R')
  
  dat = 
    bci |>
    filter(census == 7)
  
  network = 
    dat %>%
    adjacency_matrix_parallelDist(
      autolinked_species_only = FALSE, 
      d_cutoff = 20, 
      d_step = 1e-5, 
      Lx = 1000, 
      Ly = 500
    ) %>%
    cluster_analysis(
      algorithm = 'louvain',
      weighted = TRUE,
      self_loops = FALSE
    )
  
  result = 
    network %>%
    pluck('result')
  
}

if(do.kde.analysis){
  
  #label the folder "bci_dryad.zip"
  
  raw.files = paste0('https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci.tree',1:8,'.rdata?raw=true')
  dir<-getwd()
  setwd(paste0(dir,"/","bci_dryad"))
  raw.files<-list.files()
  #Load data from DRYAD datasets: Census 7
  mydat<- lapply(raw.files, function(x) {
    load(file = url(x))
    get(ls()[ls()!= "filename"])
  })
  setwd(dir)
  all <- do.call("rbind", mydat)%>%tibble()
  bci<-all%>%mutate(census=rep(1:8,sapply(mydat, nrow)))
  bci<-bci%>%select(sp,gx,gy,dbh,census)%>%drop_na()
  
  
  prefix = 'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci.full'
  suffix = '.rdata?raw=true'
  bci = NULL
  for(census in 1:7){
    bci_raw = 
      get(load(url(paste0(prefix, census, suffix)))) %>%
      as_tibble() %>%
      drop_na(dbh)
    
    bci_sp = 
      bci_raw %>% 
      group_by(sp) %>% 
      summarize(
        baldeck_cutoff = quantile(dbh, .44),
        baldeck_n = sum(dbh > baldeck_cutoff),
        .groups = 'drop'
      ) %>% 
      inner_join(
        bci_raw %>% 
          count(sp)
      )
    
    bci = 
      bci %>%
      bind_rows(
        bci_raw %>%
          inner_join(bci_sp) %>%
          filter(dbh >= baldeck_cutoff) %>%
          select(sp, gx, gy) %>%
          mutate(census = census)  
      )
  }
  
  census_data = bci
  
  consistent_cluster_data = readRDS('~/SpatialNiche/Data/20211022/bci_clustering_analysis_consistent.rds')
    
  data = 
    census_data %>% 
    select(sp, census, gx, gy) %>%
    inner_join(
      consistent_cluster_data %>%
        filter(seed == 0) %>%
        select(sp, census, d_cutoff, group)
    ) %>%
    select(census, d_cutoff, gx, gy, group)
  
  
  kde_full = 
    expand_grid(
      census = unique(data$census),
      d_cutoff = unique(data$d_cutoff)
    ) %>%
    future_pmap_dfr(
      .f = KDE,
      .data = data,
      .options = furrr_options(seed = NULL)
    )
  
  kde_full %<>%
    inner_join(
      kde_full %>%
        group_by(census, d_cutoff, x, y) %>%
        slice_max(density, n = 1, with_ties = FALSE) %>%
        rename(soiltype = group) %>%
        ungroup() %>%
        select(-density)
    ) %>%
    mutate(fdp = fdp) %>%
    mutate(soiltype = factor(soiltype, levels = c(4, 2, 3, 1, 5)))
      
}

## Definition of theta := P(recruit | match) / P(recruit | !match)
if(do.recruitment.analysis){
  
  library(sparr) # for function bivariate.density() in KDE()
  filter = dplyr::filter
  
  ## Determine whether working on SeaWulf (SBU hpc) or personal computer
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 4
  plan(multisession, workers = cores)
  
  bci = NULL
  for(census in 1:7){
    bci_raw = 
      get(load(url(paste0(prefix, census, suffix)))) %>%
      as_tibble() %>%
      drop_na(dbh)
    
    bci_sp = 
      bci_raw %>% 
      group_by(sp) %>% 
      summarize(
        baldeck_cutoff = quantile(dbh, .44),
        baldeck_n = sum(dbh > baldeck_cutoff),
        .groups = 'drop'
      ) %>% 
      inner_join(
        bci_raw %>% 
          count(sp)
      )
    
    bci = 
      bci %>%
      bind_rows(
        bci_raw %>%
          inner_join(bci_sp) %>%
          filter(dbh >= baldeck_cutoff) %>%
          select(sp, gx, gy, dbh) %>%
          mutate(census = census)  
      )
  }
  
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
      cluster_data %>%
        filter(seed == 0), 
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
    mutate(group = factor(group)) %>%
    left_join(z_full) %>%
    mutate(group = reference) %>%
    select(-reference)
  
  cluster_data =
    cluster_data %>%
    mutate(group = factor(group)) %>%
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
    
    if(fdp=='bci'){
      nutrients %<>%
        pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
        group_by(nutrient) %>%
        mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
        ungroup
    }
    
    if(fdp == 'lap') {
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
      select(census, algorithm, d_cutoff, fdp) %>%
      unique()
    
    indices = 
      if(seawulf){
        seq(nrow(parms)) 
      }else{
        which(parms$census == 7 & parms$d_cutoff == 20) 
      } 
      
    soiltype_v_nutrients = 
      indices %>%
      future_map_dfr(
        .options = furrr_options(seed = TRUE),
        .f = function(index){ 
          data = 
            dtf %>% 
            inner_join(parms[index, ]) %>%
            {
              if (fdp == 'bci') 
                select(.,Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, `N(min)`, pH, soiltype)
              else .
            } %>% 
            {
              if (fdp == 'lap')select(.,Al, Ca, Cu, Fe, K, Mg, Mn, P, Zn, pH, soiltype)
              else .
            } %>%
            mutate(soiltype = factor(soiltype, levels = unique(soiltype)))
          
          C5_model = 
            train(
              soiltype ~ ., 
              data = data, 
              method = 'C5.0',
              trControl = trainControl(method = 'repeatedcv', repeats = 10),
              metric = 'Kappa'
            )
          
          C5_model$results %>%
            as_tibble %>%
            bind_cols(parms[index, ]) %>%
            return()
        }
      )
    
    if(seawulf){
      dir.create(save_directory, showWarnings = FALSE)
      saveRDS(
        soiltype_v_nutrients, 
        file = paste0(save_directory, fdp, '_C5_soiltype_vs_nutrients.rds')
      )
    }
    
    
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

## confirmatory factor analysis
if(do.cfa.analysis){
  
  library(lavaan)
  library(GGally)
  
  plants = 
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/20211205/20211205_kde_full.rds?raw=true'
      )
    ) |>
    filter(
      census == 8,
      d_cutoff == 20
    ) |>
    select(
      x, y, group, density
    ) |>
    pivot_wider(names_from = group, values_from = density)
  
  names(plants) = c('x', 'y', 'g1', 'g2', 'g3', 'g4')
  
  elevation = 
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bcielevation.rds?raw=true'
      )
    ) |> 
    mutate(
      x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)], 
      y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
    ) |> 
    replace_na(list(x = 10, y = 10)) |> 
    group_by(x, y) |> 
    summarize(elevation = mean(elev), .groups = 'drop') 
  
  water = 
    read.table(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
      ),
      header = TRUE
    ) |>
    as_tibble() |>
    mutate(
      x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)], 
      y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
    ) |> 
    replace_na(list(x = 10, y = 10)) |> 
    group_by(x, y) |> 
    summarize(water = mean(swp), .groups = 'drop') 
  
  nutrients =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bci_nutrient_data.rds?raw=true'
      )
    )
  
  pc_nutrients = 
    nutrients |> 
    select(B, Ca, Cu, Fe, K, Mg, Mn, N, `N(min)`, Zn) |> 
    pcaMethods::pca(nPcs = 1, scale = 'uv', center = TRUE)
  
  pc_cations = 
    nutrients |> 
    select(B, Ca, Cu, Fe, K, Mg, Mn, N, `N(min)`, Zn, Al, P) |> 
    pcaMethods::pca(nPcs = 1, scale = 'uv', center = TRUE)
  
  pc_chemistry = 
    nutrients |> 
    select(B, Ca, Cu, Fe, K, Mg, Mn, N, `N(min)`, Zn, Al, P, pH) |> 
    pcaMethods::pca(nPcs = 1, scale = 'uv', center = TRUE)
  
  data = 
    nutrients |>
    full_join(water) |>
    full_join(elevation) |>
    full_join(plants) |>
    bind_cols(
      pc_nutrients@scores |> 
        as_tibble() |> 
        transmute(nutrients = scale(PC1)[,1])
    ) |>
    bind_cols(
      pc_cations@scores |> 
        as_tibble() |> 
        transmute(cations = scale(PC1)[,1])
    ) |>
    bind_cols(
      pc_chemistry@scores |> 
        as_tibble() |> 
        transmute(chemistry = scale(PC1)[,1])
    ) |> 
    mutate(
      across(
        Al:chemistry, 
        function(x) scale(x)[,1]
      )
    )
  
  x = psych::fa.parallel(data |> select(Al:water)) # suggests 4 components
  
  pc_soil = 
    data |>
    select(Al:water) |>
    pcaMethods::pca(nPcs = 4, scale = 'uv', center = TRUE)
  
  pc_soil@loadings |> 
    as_tibble() |> 
    mutate(feature = rownames(pc_soil@loadings)) |> 
    ggplot(aes(abs(PC1), abs(PC2), label = feature)) + 
    geom_text()
    
  
  # pc_loadings =
  #   pc@loadings |>
  #   as_tibble() |>
  #   mutate(nutrient = rownames(pc@loadings)) |>
  #   pivot_longer(-nutrient, names_to = 'component')

  # pc_loadings|>
  #   ggplot(aes(nutrient, abs(value), fill = nutrient)) +
  #   geom_col() +
  #   facet_wrap(~component) +
  #   theme(legend.position = 'none')
  #
  # pc_loadings|>
  #   ggplot(aes(component, abs(value), fill = component)) +
  #   geom_col() +
  #   facet_wrap(~nutrient) +
  #   theme(legend.position = 'none')
  
  models = 
    list(
       
        '
          g1 + g2 + g3 + g4 ~ nutrients + water + P + Al + pH
        ',
        
        '
          nutrients + P ~ g1 + g2 + g3 + g4 + water + pH
        ',
      
        '
          nutrients + P ~ g1 + g2 + g3 + g4
        ',
        
        
        
        '
          g1 + g2 + g3 + g4 ~ cations + water + pH
        ',
        
        '
          cations ~ g1 + g2 + g3 + g4 + water + pH
        ',
        
        '
          cations ~ g1 + g2 + g3 + g4
        ',
        
        
        
        '
          g1 + g2 + g3 + g4 ~ chemistry + water
        ',
        
        '
          chemistry ~ g1 + g2 + g3 + g4 + water
        ',
        
        '
          chemistry ~ g1 + g2 + g3 + g4
        ',
        
        
        '
          correlated_nutrients =~ B + Ca + Cu + Fe + K + Mg + Mn + Zn + N + `N(min)`
          g1 + g2 + g3 + g4 ~ correlated_nutrients + Al + P + pH + water
        ',

        '
          correlated_nutrients =~ B + Ca + Cu + Fe + K + Mg + Mn + Zn + N + `N(min)`
          correlated_nutrients + P ~ g1 + g2 + g3 + g4 + pH + water
        ',
        
        '
          LV1 =~ B + Ca + Cu + Fe + K + Mg + Mn + `N(min)` + Zn
          LV2 =~ P + N + elevation
          LV3 =~ water + Al
          LV4 =~ pH
          
          g1 + g2 + g3 + g4 ~ LV1 + LV2 + LV3 + LV4
        ',
        
        '
          LV1 =~ B + Ca + Cu + K + Mg + `N(min)` + Zn
          LV2 =~ Fe + Mn + pH
          LV3 =~ P + water + Al
          LV4 =~ N
          
          g1 + g2 + g3 + g4 ~ LV1 + LV2 + LV3 + LV4
        '
       
    )
  
  names(models) = paste0('model', 1:length(models))
  
  fit = lapply(models, cfa, data = data)
  
  measures = sapply(fit, fitMeasures)
  measures = 
    measures |>
    as_tibble() |>
    mutate(measure = rownames(measures)) |>
    select(measure, everything()) |>
    filter(measure %in% c('chisq', 'cfi', 'tli', 'rfi', 'aic','rmsea', 'rmr', 'srmr')) |>
    pivot_longer(-measure, names_to = 'model')
  
  plot_measures = 
    measures |> 
    ggplot(aes(value, model, fill = model)) + 
    geom_col() + 
    facet_wrap(~measure, scales = 'free') +
    theme(legend.position = 'none')
  
  coefs =
    # matrix(coef(fit[[10]])[1:20], 4, ) |>
    matrix(coef(fit[[10]])[11:30], 4, ) |>
    as_tibble() |>
    mutate(group = paste0('g', 1:4)) |>
    rename(nutrients = V1, water = V2, P = V3, Al = V4, pH = V5) |>
    select(group, everything())

  plot_pairwise =
    coefs |>
    ggpairs(columns = 2:ncol(coefs)) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')

  plot_normality = 
    data |> 
    pivot_longer(-c(x,y), names_to = 'index') |> 
    ggplot(aes(value)) + 
    facet_wrap(~index) + 
    geom_histogram(binwidth = .1) + 
    geom_line(
      aes(x, y), 
      data = 
        tibble(
          x = seq(-5,5,.01), 
          y = 1250 *.1 * dnorm(x)
        ), 
      color = 'red'
    ) +
    ylab('Density')
  
  plot1 =
    data |> 
    select(x, y, g1, g2, g3, g4) |> 
    pivot_longer(-c(x, y), names_to = 'group', values_to = 'density') |>
    ggplot(aes(x, y, fill = density)) +
    geom_tile() +
    theme(aspect.ratio = .5) +
    scale_fill_gradientn(colors = terrain.colors(2000)) +
    facet_wrap(~group)
  
  plot2 = 
    plants |> 
    pivot_longer(
      -c(x, y), 
      names_to = 'group', 
      values_to = 'density'
    ) |> 
    group_by(x, y) |> 
    slice_max(density) |> 
    ungroup() |> 
    ggplot(aes(x, y, fill = group)) + 
    geom_tile() + 
    theme(aspect.ratio = .5)
}

## test C5.0 on Gaussian random field
if(do.C5.analysis){
  library(RandomFields)
  library(gstat)
  library(caret)
  
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 4
  plan(multisession, workers = cores)
  
  analysis = 
    function(
      number_of_nulls,
      kde_full,
      nutrients,
      clusters_run,
      save.result = clusters_run
    ){
      
      dtf_wide = 
        kde_full %>%
        group_by(x, y) %>%
        slice_max(density, with_ties = FALSE) %>%
        ungroup() %>%
        inner_join(nutrients)
      
      if(fdp == 'bci'){
        soil_dtf = 
          dtf_wide %>%
          select(Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, `N(min)`, pH, soiltype)
      }
      
      if(fdp == 'lap'){
        soil_dtf = 
          dtf_wide %>%
          select(Al, Ca, Cu, Fe, K, Mg, Mn, P, Nmin, pH, soiltype)
      }
      
      RandomField = 
        function(
          Lx, 
          Ly, 
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed = seed, 
          plot = FALSE
        ) {
          stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
          
          RFoptions(seed = seed)
          
          stress = 
            RFsimulate(
              RMgauss(
                scale = rangepar + 1e-16,
                var = 2 * sillpar * (1 - nuggetpar)
              ) + RMtrend(mean = 0) + RMnugget(var = 2 * sillpar * nuggetpar),
              x = 1:Lx,
              y = 1:Ly
            )@data$variable1
          
          if(plot){
            plot(
              raster::raster(matrix(stress, Lx, Ly)),
              las = 1,
              xlab = 'x-coordinate',
              ylab = 'y-coordinate'
            )
          }
          
          return(stress)
          
        }
      
      Gaussian_vgm_optim = 
        function(parms, sample.vgm){
          range = parms[1]
          sill = parms[2]
          nugget = parms[3]
          dist = c(0, sample.vgm$dist)
          observed = c(0, sample.vgm$gamma)
          predicted = (sill - nugget) * (1 - exp(-dist ^ 2 / range ^ 2)) + nugget * (dist > 0)
          sum_squared_errors = sum((observed - predicted) ^ 2)
          return(sum_squared_errors)
        } 
      
      sample.vgm = 
        variogram(
          soiltype ~ 1, 
          data = dtf_wide, 
          locations = ~ x + y, 
          width = 1
        )
      
      fitted.vgm = 
        optim(
          par = c(range = 20, sill = .6, nugget = 20), 
          fn = Gaussian_vgm_optim,
          sample.vgm = sample.vgm,
          method = "L-BFGS-B",
          lower = c(1, 1, 0)
        )$par
      
      sill = as.numeric(fitted.vgm[2])		## semivariance at the landscape level --> regional heterogeneity
      nugget = as.numeric(fitted.vgm[3])	## semivariance at distance = 1	--> local uniformity
      range = as.numeric(fitted.vgm[1])		## distance at which the semivariance reaches 63% of the sill
      range95 = sqrt(3) * range   ## distance at which the semivariance reaches 95% of the sill
      
      randomfield_test = 
        function(
          seed, 
          Lx,
          Ly,
          date,
          .data, 
          .range95, 
          .sill, 
          .func
        ){
          set.seed(seed)
          
          test_data = 
            expand_grid(
              y = seq(10, Ly - 10, by = 20),
              x = seq(10, Lx - 10, by = 20)
            ) %>%
            mutate(
              field = 
                .func(
                  Lx = Lx / 20, 
                  Ly = Ly / 20, 
                  rangepar = range95, 
                  sillpar = sill, 
                  nuggetpar = .015, 
                  seed = seed
                ),
              soiltype = 
                cut(
                  field, 
                  breaks = quantile(
                    field, 
                    p = c(
                      0, 
                      dtf_wide %>% 
                        count(soiltype) %>% 
                        mutate(p = cumsum(n) / sum(n)) %>% 
                        pull(p)
                    )
                  ), ## these quantiles match the relative prevalence 
                  ## of each soiltype in the real data
                  labels = FALSE,
                  include.lowest = TRUE
                ) %>%
                as.factor
            )
          
          C5_model =
            train(
              soiltype ~ .,
              data = .data %>% mutate(soiltype = test_data$soiltype),
              method = 'C5.0',
              trControl = trainControl(method = 'repeatedcv', repeats = 10),
              metric = 'Kappa'
            )
          
          res = 
            bind_cols(
              seed = seed,
              C5_model$results
            ) %>%
            as_tibble
          
          return(res)
          
        }
      
      results = 
        1:number_of_nulls %>%
        future_map_dfr(
          .f = randomfield_test,
          Lx = ifelse(fdp == 'bci', 1000, 500),
          Ly = 500,
          date = date,
          .data = soil_dtf,
          .sill = sill,
          .range95 = range95,
          .func = RandomField,
          .options = furrr_options(seed = NULL)
        )
      
      if(save.result){
        
        save_date = gsub('-', '', Sys.Date())
        save_directory = paste0('~/SpatialNiche/Data/', save_date)
        dir.create(save_directory, showWarnings = FALSE)
        saveRDS(
          results,
          file = paste0(save_directory, '/', fdp, '_C5_randomfields.rds')
        ) 
      }
      
      return(results)
      
    }
  
  results = 
    analysis(
      number_of_nulls = ifelse(seawulf, 100, 1),
      kde_full = kde_full,
      nutrients = nutrients,
      clusters_run = seawulf
    )
  
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
  
  nutrients %<>%
    pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
    group_by(nutrient) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup
  
  all_data = 
    census_data %>%
    filter(dbh >= 100) %>%
    inner_join(
      cluster_data %>% 
        filter(seed == 0)
    ) %>%
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
    all_data %>% 
    group_by(trait_type, nutrient) %>% 
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
    cor_analysis %>% 
    filter(census == 7, d_cutoff == 20, !nutrient %in% c('Al', 'pH')) %>% 
    group_by(soiltype) %>% 
    summarize(mean_cor = mean(cor), .groups = 'drop') %>%
    mutate(
      group_ID = 
        factor(
          soiltype, 
          levels = soiltype[order(mean_cor, decreasing = TRUE)]
        )
    )
  
  all_data %<>%
    inner_join(
      mean_correlation %>%
        select(group = soiltype, group_ID)
    )
  
  cor_analysis %<>%
    inner_join(
      mean_correlation %>%
        select(soiltype, group_ID)
    )  
  
  species_group_v_trait =
    all_data %>% 
    select(sp, trait_type, pc1, group_ID) %>% 
    unique 
  
  species_group_v_trait_wide = 
    species_group_v_trait %>% 
    pivot_wider(names_from = trait_type, values_from = pc1) %>%
    arrange(sp)
  
  abundances =
    all_data %>% 
    select(sp, gx, gy) %>% 
    unique() %>%
    count(sp)
  
  group_levels = 
    species_group_v_trait %>% 
    pull(group_ID) %>% 
    levels()
  
  plot_correlations = 
    correlations %>% 
    filter(pval <= .05) %>%
    ggplot(aes(trait_type, nutrient, fill = cor)) + 
    geom_tile() + 
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue')
  
  plot_correlations %>%
    show()
  
  plot_group_densities = 
    all_data %>% 
    filter(trait_type == 'vital', nutrient == 'Fe') %>% 
    ggplot(aes(gx, gy)) + 
    geom_density_2d_filled() + 
    theme(aspect.ratio = .5) + 
    facet_wrap(~ group_ID)
  
  plot_trait_violins_trees = 
    all_data %>% 
    filter(nutrient == 'Fe') %>% 
    ggplot(aes(group_ID, pc1, fill = group_ID)) + 
    geom_violin(draw_quantiles = .5) + 
    facet_wrap(~ trait_type)
  
  plot_trait_violins_species = 
    all_data %>% 
    filter(nutrient == 'Fe') %>% 
    select(sp, group_ID, pc1, trait_type) %>%
    unique() %>%
    ggplot(aes(group_ID, pc1, fill = group_ID)) + 
    geom_violin(draw_quantiles = .5) + 
    facet_wrap(~ trait_type)
  
  plot_group_nutrient_correlations = 
    cor_analysis %>% 
    filter(census == 7, d_cutoff == 20) %>% 
    ggplot(aes(nutrient, cor, fill = group_ID)) + 
    geom_col() + 
    facet_wrap(~ group_ID)
  
  plot_trait_space = 
    species_group_v_trait_wide %>%
    ggplot(aes(vital, wood, color = group_ID)) +
    geom_point(size = 4) +
    theme(aspect.ratio = 1)
  
  
  table = 
    species_group_v_trait %>%
    filter(trait_type %in% c('vital', 'wood', 'leaf')) %>%
    group_by(trait_type) %>% 
    summarize(
      group_ID_x = rep(group_levels[-1], 3),
      group_ID_y = rep(group_levels[-length(group_levels)], each = 3),
      pval = as.numeric(pairwise.wilcox.test(pc1, group_ID)$p.value),
      .groups = 'drop'
    )

  
  
}  

if(do.paper.figures){
  
  do.fig1 = 1
  do.fig2 = 1
  do.fig3 = 1
  do.fig4 = 1
  do.table = 1
  do.SI.figs = 1
  
  selection = tibble(census = 7, d_cutoff = 20)
  
  plotting_colors = 
    c(
      red = "#DD5144", 
      green = "#1DA462", 
      blue = "#4C8BF5", 
      yellow = "#FFCD46"
    )
  
  group_levels =
    nutrients %>%
    pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
    group_by(nutrient) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup %>%
    select(x, y, nutrient, standardized) %>%
    inner_join(
      kde_full %>%
        inner_join(selection) %>%
        select(-soiltype) %>%
        rename(soiltype = group), 
      by = c('x', 'y')
    ) %>%
    group_by(
      nutrient,
      soiltype
    ) %>%
    summarize(
      cor = cor(density, standardized, use = 'complete.obs'), 
      .groups = 'drop'
    ) %>% 
    filter(!nutrient %in% c('Al', 'pH')) %>% 
    group_by(soiltype) %>% 
    summarize(mean_cor = mean(cor), .groups = 'drop') %>%
    mutate(
      group_ID = 
        factor(
          soiltype, 
          levels = soiltype[order(mean_cor, decreasing = TRUE)]
        )
    )
  
  # data = 
  #   'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bci_figure_data.rds?raw=true' %>%
  #   url() %>%
  #   readRDS() %>%
  #   mutate(group = factor(group, levels = levels(group_levels$group_ID)))
  
  data = 
    readRDS('~/SpatialNiche/Data/20211004/bci_figure_data.rds') %>%
    mutate(group = factor(group, levels = levels(group_levels$group_ID)))
  
  species_names = 
    '~/BCI/BCITRAITS_20101220.xlsx' %>%
    readxl::read_excel() %>%
    mutate(
      sp = tolower(`SP$`),
      genus = substr(`GENUS$`, start = 1, stop = 1),
      genus = paste0(genus, '.'),
      species = `SPECIES$`,
      name = paste(genus, species)
    ) %>%
    select(sp, name) %>%
    mutate(name = ifelse(sp == 'swars1', 'S. simplex_var1', name)) %>%
    mutate(name = ifelse(sp == 'swars2', 'S. simplex_var2', name))
  
  species_by_group = 
    cluster_data %>%
    filter(seed == 0) %>%
    inner_join(selection) %>%
    select(sp1 = sp, group) %>%
    mutate(
      group = factor(group, levels = levels(group_levels$group_ID)),
      sp1 = factor(sp1, levels = sp1[order(group)])
    )
  
  
  
  if(do.fig1){
    ## Figure 1A. Geographic location of trees
    dtf1a = 
      census_data %>%
      select(sp, gx, gy, dbh) %>%
      filter(dbh >= 100) %>%
      inner_join(
        cluster_data %>%
          inner_join(selection) %>%
          select(sp, group)
      ) %>%
      mutate(
        group = factor(group, levels = levels(group_levels$group_ID))
      ) %>%
      drop_na()
    
    fig1a = 
      dtf1a %>%
      ggplot(aes(gx, gy, color = group)) +
      geom_point() +
      theme(aspect.ratio = .5) +
      labs(x = 'x coordinate (m)', y = 'y coordinate (m)') +
      scale_color_manual(
        breaks = unique(dtf1a$group), 
        values = as.character(plotting_colors)
      ) +
      theme(legend.position = 'none')
    
    
    ## Figure 1B. Adjacency matrix
    source('~/SpatialNiche/R_Scripts/clustering_functions.R')
    
    adjacency_tibble =
      census_data %>%
      filter(dbh >= 100)%>%
      inner_join(selection) %>%
      adjacency_matrix(
        autolinked_species_only = TRUE, 
        d_cutoff = selection$d_cutoff, 
        d_step = 1e-5, 
        Lx = 1000, 
        Ly = 500
      ) %>%
      pluck('adjacency') %>%
      inner_join(species_by_group) %>%
      mutate(
        sp1 = factor(sp1, levels = levels(species_by_group$sp1)),
        sp2 = factor(sp2, levels = levels(species_by_group$sp1)),
        binary = factor(binary, levels = c(0, 1))
      )
    
    plotted_names = 
      species_names %>%
      right_join(
        species_by_group %>% 
          rename(sp = sp1)
      ) %>%
      replace_na(list(name = 'T. integerrima')) %>%
      arrange(name) %>%
      mutate(name = factor(name, levels = name[order(group)])) %>%
      arrange(name)
    
    fig1b = 
      adjacency_tibble %>%
      ggplot(aes(sp1, sp2, fill = binary)) +
      geom_tile() +
      theme(aspect.ratio = 1) +
      labs(x = '', y = 'species') +
      theme(legend.position = 'none') +
      scale_fill_manual(values = c('white', 'black')) +
      theme(axis.text.x = element_blank()) +
      scale_y_discrete(labels = plotted_names$name) +
      theme(axis.text.y = element_text(color = plotting_colors[plotted_names$group]))
    
    fig1 = cowplot::plot_grid(fig1a, fig1b, labels = 'AUTO')
  }
  
  if(do.fig2){
    
    fig2a = 
      data %>%
      filter(trait_type == 'vital', nutrient == 'Al') %>%
      select(gx, gy, group) %>%
      mutate(group_name = names(plotting_colors[group])) %>%
      mutate(
        group_name = 
          factor(
            group_name, 
            levels = unique(group_name[order(group)])
          )
      ) %>%
      ggplot(aes(gx, gy)) +
      geom_density_2d_filled() +
      facet_wrap(~group_name) +
      theme(
        aspect.ratio = .5,
        legend.position = 'drop'
      ) +
      labs(x = 'x coordinate (m)', y = 'y coordinate (m)')
    
    dtf2b = 
      kde_full %>%
      inner_join(selection) %>%
      select(x, y, soiltype) %>%
      unique() %>%
      mutate(
        soiltype = 
          factor(
            soiltype, 
            levels = levels(group_levels$group_ID)
          )
      ) %>%
      arrange(soiltype)
    
    fig2b = 
      dtf2b %>%
      ggplot(aes(x, y, fill = soiltype)) +
      geom_tile() +
      scale_fill_manual(
        breaks = unique(dtf2b$soiltype), 
        values = as.character(plotting_colors)
      ) +
      theme(
        legend.position = 'none',
        aspect.ratio = .5
      ) +
      labs(x = 'x coordinate (m)', y = 'y coordinate (m)')
    
    fig2 = 
      cowplot::plot_grid(
        fig2a, 
        fig2b,
        labels = 'AUTO'
      )
      
  }
  
  if(do.fig3){
    
    dtf3a = 
      data %>%
      filter(trait_type == 'vital') %>%
      select(x, y, nutrient, standardized) %>%
      unique() %>%
      pivot_wider(names_from = nutrient, values_from = standardized)
    
    pca_model = 
      dtf3a %>%
      select(-c(x, y)) %>%
      pca(method = 'ppca', scale = 'none', center = FALSE)
    
    dtf3a %<>%
      bind_cols(
        pc1 = pca_model@scores[, 1],
        pc2 = pca_model@scores[, 2]
      ) %>%
      inner_join(
        kde_full %>%
          inner_join(selection) %>%
          select(x, y, soiltype) %>%
          unique() %>%
          mutate(
            soiltype = 
              factor(
                soiltype, 
                levels = levels(group_levels$group_ID)
              )
          ) %>%
          arrange(soiltype)
      ) %>%
      arrange(soiltype)
    
    fig3a = 
      dtf3a %>%
      ggplot(aes(pc1, pc2, color = soiltype)) +
      geom_point() +
      theme(
        aspect.ratio = 1,
        legend.position = 'none'
      ) +
      scale_color_manual(
        breaks = unique(dtf3a$soiltype), 
        values = as.character(plotting_colors)
      ) +
      labs(x = 'nutrient PC1', y = 'nutrient PC2')
    
    dtf3b = 
      kde_full %>%
      inner_join(selection) %>%
      mutate(group = factor(group, levels = levels(group_levels$group_ID))) %>%
      inner_join(
        data %>%
          filter(trait_type == 'vital')
      ) %>%
      group_by(group, nutrient) %>%
      summarize(
        correlation = cor(density, standardized),
        p.value = cor.test(density, standardized)$p.value,
        signif = p.value <= .05,
        .groups = 'drop'
      ) %>%
      arrange(group) %>%
      mutate(group_name = names(plotting_colors[group])) %>%
      mutate(
        group_name = 
          factor(
            group_name, 
            levels = unique(group_name[order(group)])
          )
      )
    
    fig3b = 
      dtf3b %>%
      ggplot(aes(nutrient, correlation, fill = group_name)) +
      geom_col(color = 'black') +
      facet_wrap(~ group_name) +
      scale_fill_manual(
        breaks = unique(dtf3b$group_name), 
        values = as.character(plotting_colors)
      ) +
      theme(
        legend.position = 'none',
        aspect.ratio = 1
      )
    
    fig3 = cowplot::plot_grid(fig3a, fig3b, labels = 'AUTO')
      
  }
  
  if(do.fig4){
    dtf4 = 
      data %>%
      filter(
        nutrient == 'P',
        trait_type %in% c('vital', 'leaf', 'wood')
      ) %>%
      select(sp, group, trait_type, pc1) %>%
      unique() %>%
      pivot_wider(names_from = trait_type, values_from = pc1) %>%
      arrange(group)
    
    fig4 = 
      dtf4 %>%
      ggplot(aes(vital, wood, color = group)) +
      geom_point(size = 4) +
      theme(
        aspect.ratio = 1,
        legend.position = 'none'
      ) +
      scale_color_manual(
        breaks = unique(dtf4$group),
        values = as.character(plotting_colors)
      ) +
      labs(x = 'vital rates PC1', y = 'wood density PC1')
    
  }
  
  if(do.table){
    table = 
      data %>%
      filter(
        nutrient == 'P',
        trait_type %in% c('vital', 'leaf', 'wood')
      ) %>%
      select(sp, group, trait_type, pc1) %>%
      unique() %>%
      mutate(group_name = names(plotting_colors)[group]) %>%
      group_by(trait_type) %>% 
      summarize(
        group_x = rep(names(plotting_colors)[-1], 3),
        group_y = rep(names(plotting_colors)[-length(names(plotting_colors))], each = 3),
        pval = as.numeric(pairwise.wilcox.test(pc1, group)$p.value),
        .groups = 'drop'
      ) %>%
      pivot_wider(names_from = group_y, values_from = pval)
  }
  
  if(do.SI.figs){
    
    ## spatial density of nutrients
    sfig_nutrients = 
      data %>%
      filter(trait_type == 'vital') %>%
      ggplot(aes(x, y, fill = standardized)) +
      geom_tile() +
      facet_wrap(~nutrient, nrow = 3) +
      scale_fill_gradientn(colors = terrain.colors(1000)) +
      theme(
        aspect.ratio = .5, 
        legend.position = 'none'
      ) +
      labs(x = 'x coordinate (m)', y = 'y coordinate (m)')
    
    ## violin plots
    sfig_violins_by_trees = 
      data %>%
      filter(
        nutrient == 'P', 
        trait_type %in% c('vital', 'wood', 'leaf')
      ) %>%
      ggplot(aes(group, pc1, fill = group)) +
      geom_violin() +
      facet_wrap(~ trait_type, scales = 'free') +
      scale_fill_manual(
        breaks = c(2, 4, 1, 3),
        values = as.character(plotting_colors)
      ) +
      theme(
        legend.position = 'none',
        aspect.ratio = 1,
        axis.text.x = element_blank()
      ) +
      labs(x = '', y = 'Trait PC1')
    
    sfig_violins_by_species = 
      data %>%
      filter(
        nutrient == 'P', 
        trait_type %in% c('vital', 'wood', 'leaf')
      ) %>%
      select(group, pc1, trait_type) %>%
      unique() %>%
      ggplot(aes(group, pc1, fill = group)) +
      geom_violin() +
      facet_wrap(~ trait_type, scales = 'free') +
      scale_fill_manual(
        breaks = c(2, 4, 1, 3),
        values = as.character(plotting_colors)
      ) +
      theme(
        legend.position = 'none',
        aspect.ratio = 1,
        axis.text.x = element_blank()
      ) +
      labs(x = '', y = 'Trait PC1')
    
    sfig_violins = 
      cowplot::plot_grid(
        sfig_violins_by_trees,
        sfig_violins_by_species,
        nrow = 2,
        labels = 'AUTO'
      ) 
    
    ## correlation between traits and trait type PC1
    sfig_traits_v_PC1 = 
      'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bci_trait_data_processed.rds?raw=true' %>%
      url() %>%
      readRDS() %>% 
      filter(trait_type %in% c('vital', 'wood', 'leaf')) %>% 
      ggplot(aes(standardized, pc1)) + 
      geom_point() + 
      facet_wrap(trait_type ~ trait, ncol = 8) + 
      theme(aspect.ratio = 1) +
      labs(
        x = 'standardized trait value',
        y = 'trait type PC1'
      )
    
    ## correlations between PC1 of each trait type
    dtfsfig = 
      data %>% 
      filter(
        nutrient == 'P', 
        trait_type %in% c('vital', 'wood', 'leaf')
      ) %>% 
      select(sp, trait_type, pc1) %>% 
      unique() %>% 
      left_join(
        x = .,
        y = .,
        by = 'sp'
      ) %>% 
      filter(trait_type.x > trait_type.y)
    
    sfig_trait_type_PC1 = 
      dtfsfig %>%
      ggplot(aes(pc1.x, pc1.y)) + 
      geom_smooth(method = 'lm') +
      geom_point() + 
      facet_grid(trait_type.x ~  trait_type.y, scales = 'free') + 
      theme(aspect.ratio = 1) +
      labs(x = 'trait type PC1', y = 'trait type PC1')
    
    correlations_pvalues = 
      dtfsfig %>% 
      group_by(trait_type.x, trait_type.y) %>% 
      summarize(
        cor = cor(pc1.x, pc1.y), 
        pval = cor.test(pc1.x, pc1.y)$p.value, 
        .groups = 'drop'
      )
    
    ## number of groups per census x d_cutoff
    sfig_ngroups = 
      cluster_data %>% 
      filter(seed == 0) %>%
      filter(d_cutoff >= 20) %>%
      group_by(census, d_cutoff) %>% 
      summarize(
        ngroups = unique(number_of_groups), 
        .groups = 'drop'
      ) %>%
      mutate(`# groups` = factor(ngroups)) %>%
      ggplot(aes(d_cutoff, census, fill = `# groups`)) +
      geom_tile(color = 'black') +
      theme(aspect.ratio = 1) +
      labs(x = 'distance cutoff')
    
    ## soil-type raster per census x d_cutoff
    plotting_colors5 = 
      c(
        red = "#DD5144",
        plum = 'plum4',
        green = "#1DA462",
        blue = "#4C8BF5",
        yellow = "#FFCD46")
      
    sfig_soiltypes = 
      kde_full %>%
      filter(d_cutoff >= 20) %>%
      select(x, y, census, d_cutoff, soiltype) %>%
      unique() %>%
      # mutate(
      #   soiltype = 
      #     factor(
      #       soiltype, 
      #       levels = c(1, 4, 2, 5, 3)
      #     )
      # ) %>%
      ggplot(aes(x, y, fill = soiltype)) +
      geom_tile() +
      facet_grid(d_cutoff ~ census) +
      theme(
        aspect.ratio = .5,
        legend.position = 'none'
      ) +
      labs(x = 'x coordinate (m)', y = 'y coordinate (m)') #+
      # scale_fill_manual(
      #   breaks = c(1, 4, 2, 5, 3), 
      #   values = as.character(plotting_colors5)
      # )
    
    ## correlation between group density and nutrient concentration per census x d_cutoff
    groupnames = c('red', 'purple', 'green', 'blue', 'yellow')
    
    dtf_sfig = 
      data %>%
      select(x, y, nutrient, standardized) %>%
      unique() %>%
      inner_join(
        kde_full %>%
          select(census, d_cutoff, x, y, group, density)
      ) %>%
      filter(d_cutoff >= 20) %>%
      group_by(nutrient, group, d_cutoff) %>%
      summarize(
        cor = cor(standardized, density), 
        .groups = 'drop'
      )
    
    dtf_sfig %<>%
      mutate(
        group = 
          factor(
            group,
            levels = 
              dtf_sfig %>% 
              filter(d_cutoff == 20) %>% 
              group_by(group) %>% 
              summarize(
                mean_cor = mean(cor), 
                .groups = 'drop'
              ) %>% 
              arrange(desc(mean_cor)) %>% 
              pull(group)
          )
      ) %>%
      mutate(group = groupnames[group]) %>%
      mutate(group = factor(group, levels = groupnames))
      
    
    sfig_bars = 
      dtf_sfig %>% 
      ggplot(aes(nutrient, cor, fill = group)) + 
      geom_col(color = 'black') + 
      facet_grid(group ~ d_cutoff) +
      theme(
        aspect.ratio = 1,
        legend.position = 'none'
      ) +
      scale_fill_manual(
        breaks = groupnames, 
        values = as.character(plotting_colors5)
      ) +
      ylab('correlation')
    
    ## average theta across censuses by d_cutoff
    recruits = 
      census_data %>%
      filter(dbh >= 100) %>%
      select(census, treeID, sp, gx, gy, dbh) %>%
      filter(census > 1) %>%
      group_by(treeID) %>%
      slice_min(census) %>%
      ungroup %>%
      mutate(recruit = TRUE)
    
    rec_df = 
      recruits %>%
      mutate(
        x = seq(10, 990, 20)[cut(gx, seq(20, 1000, 20), labels = FALSE)],
        y = seq(10, 490, 20)[cut(gy, seq(20,  500, 20), labels = FALSE)]
      ) %>%
      inner_join(
        kde_full %>%
          select(
            labelling_census = census, 
            x, 
            y, 
            d_cutoff, 
            soiltype, 
            fdp
          ) %>%
          unique(),
        by = c('x', 'y')
      ) %>%
      inner_join(
        cluster_data %>%
          filter(seed == 0) %>%
          rename(labelling_census = census)
      )
    
    prior_df = 
      kde_full %>% 
      select(x, y, d_cutoff, soiltype, census) %>% 
      unique() %>% 
      count(census, d_cutoff, soiltype) %>% 
      mutate(Pm = n / (1000 * 500 / 20 ^ 2)) %>% 
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
        prior_df %>%
          rename(group = soiltype)
      ) %>%
      mutate(
        theta = ((1 - Pm) / Pm) / (recruits / matches - 1)
      ) %>%
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
    
    sfig_theta_bars = 
      res %>%
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
      theme(aspect.ratio = 1)
    
    ## modularity in data vs permutated data (where species groups are permutated)
    ## (BCI census = 7, d_cutoff = 20)
    dat = 
      "C:/Users/rdand/Google Drive/GitHub/Spatial-niche/Data/Manuscript-Data/bci_clustering_analysis_census7_dcutoff20.rds" %>%
      readRDS()
    
    all_censuses = "C:/Users/rdand/Google Drive/GitHub/Spatial-niche/Data/Manuscript-Data/bci_clustering_analysis.rds" %>%
      readRDS()
    
    null = 
      dat %>% 
      select(
        algorithm, 
        weighted, 
        seed, 
        modularity, 
        census, 
        d_cutoff
      ) %>% 
      unique() %>% 
      filter(seed > 0) %>% 
      group_by(weighted) %>% 
      summarize(
        mu = mean(modularity), 
        sigma = sd(modularity), 
        .groups = 'drop'
      )
    
     sfig_modularity_effect_size = 
      all_censuses %>% 
       filter(d_cutoff >= 20) %>%
      mutate(z = (modularity - null$mu) / null$sigma) %>% 
      ggplot(aes(factor(d_cutoff), z)) + 
      geom_boxplot(fill = 'plum4') +
      labs(
        x = 'distance cutoff (m)', 
        y = 'modularity effect size'
      ) +
       theme(aspect.ratio = 1)
     
     value = 
       all_censuses %>% 
       mutate(z = (modularity - null$mu) / null$sigma) %>% 
       filter(census == 7, d_cutoff == 20) %>%
       select(algorithm, weighted, census, d_cutoff, modularity, z) %>%
       unique()
    

  }
  
  fignames = 
    c(
      paste0('fig',1:4), 
      ls(pattern = 'sfig_*')[-c(1, 2)]
    )
  
  # for(figname in fignames){
  #   fig = get(figname)
  #   png(
  #     file = 
  #       paste0(
  #         'c:/users/rdand/Google Drive/GitHub/Spatial-niche/Manuscript/Figures/',
  #         figname,
  #         '.png'
  #       ),
  #         width = 1000,
  #         height = 500 
  #   )
  #   fig %>% show
  #   dev.off()
  # }
  # 
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

## I ended up not using this function, written for the recruitment test 
## where I find theta as a function of d_cutoff
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


## Using C5.0 on traits returns a Cohen's Kappa of 0.2 for BCI census = 7, d_cutoff = 20
## If abundance-weighing, Kappa becomes 1, presumably because the classifier learns to 
## find each species based on their traits
if(FALSE){
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
  
}
