## New script

library(wordspace) ## for function dist.matrix()
library(seriation) ## for function seriate()
library(caret) ## for function confusionMatrix()
library(ggplotify) ## to convert from base plot to grob
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
## a graph then find communities in the graph
library(gtools)  ## for function mixedsort
library(tidyverse)
library(furrr)  ## for parallel computing
library(parallel)  
library(magrittr)
# library(ks)       ## for kernel density estimation function kde() 

do.data = 0
do.plots.20210709 = 0
do.plots.20210712 = 0
do.bci.census.analysis = 0

## pdf of the distance d between two random points in a 
## square of side a, from Johan Philip 2007
euclidean2d_distance_probability = function(d, a){
  if(d > a * sqrt(2)) return(0)
  f =
    ifelse(
      d < a,
      -4 * d / a ^ 3 + pi / a ^ 2 + d ^ 2 / a ^ 4,
      -2 / a ^ 2 + 
        4 / a ^ 2 * asin(a / d) + 
        4 / a ^ 3 * sqrt(d ^ 2 - a ^ 2) - 
        pi / a ^ 2 - 
        d ^ 2 / a ^ 4
    )
  return(2 * d * f)
} 

## pdf of the distance d between two random points in a 
## rectangle of sides a and b > a, from Johan Philip 2007
euclidean2d_distance_probability_ab = function(d, a, b){
  if(a == b) return(euclidean2d_distance_probability(d, a))
  s = d ^ 2
  a2 = a ^ 2
  b2 = b ^ 2
  if(d < 0 | s > a2 + b2) return(0)
  if(a > b){
    c = a
    a = b
    b = c
  }
  if(s <= a2){
    f = - 2 * d / a2 / b - 2 * d / a / b2 + pi / a / b + s / a2 / b2
  }
  if(s > a2 & s <= b2){
    f = - 2 * d / a2 / b - 1 / b2 + 2 / a / b * asin(a / d) + 2 / a2 / b * sqrt(s - a2)
  }
  if(s > b2){
    f = 
      - 1 / b2 + 2 / a / b * asin(a / d) + 2 / a2 / b * sqrt(s - a2) - 
      1 / a2 + 2 / a / b * asin(b / d) + 2 / a / b2 * sqrt(s - b2) -
      pi / a / b - s / a2 / b2
  }
  return(2 * d * f)
}

## randomize adjacency matrix while keeping total number of edges and diagonal intact
randomize_adjacency_matrix = function(nvertices, nedges = 1) {
  edges.max = nvertices * (nvertices - 1) / 2
  
  stopifnot(length(nvertices) == 1 && 1 <= nvertices)
  stopifnot(0 <= nedges & nedges <= edges.max)
  
  index.edges = lapply(list(1:(nvertices - 1)), function(k) rep(k * (k + 1) / 2, nvertices - k)) 
  index.edges = index.edges[[1]] + 1:edges.max
  graph.adjacency = matrix(0, ncol = nvertices, nrow = nvertices)
  graph.adjacency[sample(index.edges, nedges)] = 1
  
  return(graph.adjacency + t(graph.adjacency) + diag(nvertices))
}

read_bci = function(census){
  
  df = 
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
      )
    )
  
  cen = census
  
  if(cen > 0){
    
    df %>%
      filter(census == cen) %>%
      select(-census) %>%
      return
  
  } else{
    
    df %>%
      select(-census) %>%
      unique() %>%
      return
    
  }
  
}

read_laplanada = function(census){
  
  address = 
    paste0(
      '~/SpatialNiche/Data/relaplanadacensusdata/Censo ',
      census,
      '_La Planada.txt'
    )
  
  lap = 
    read.table(
      address, 
      header = TRUE
    ) %>%
    as_tibble %>%
    filter(dbh >= 100) %>%
    select(sp, gx, gy)
  
}

adjacency_matrix = 
  function(
    dat, 
    autolinked_species_only, 
    d_cutoff, 
    d_step, 
    Lx, 
    Ly
  ){
  
  distances = seq(0, d_cutoff, by = d_step) 
  
  pdexp_fdp = 
    distances %>%
    map_dbl(
      euclidean2d_distance_probability_ab,
      a = Ly,
      b = Lx
    )
  
  ## int p(r) * dr from 0 to d*
  cumulative_null_prob_threshold = sum(pdexp_fdp * d_step) 
  
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
    filter(n >= abundance_threshold)
  
  ## selects species that are more correlated with *themselves* than
  ## a hypothetical random pair of species with equal abundance
  if(autolinked_species_only){
    selected_species =
      unique(dat_filtered$sp) %>%
      map_dfr(
        .f = function(char){
          
          df = 
            dat_filtered %>%
            filter(sp == char)
          
          d = dist(cbind(df$gx, df$gy))
          
          nn12_positive = d[d > 0]
          
          connected =
            mean(nn12_positive <= d_cutoff) -
            cumulative_null_prob_threshold
          
          return(tibble(sp = char, keep = connected > 0))
        }
      ) %>%
      filter(keep == TRUE) %>%
      pull(sp)
    
  }else{
    selected_species = unique(dat_filtered$sp)
  }
  
  
  xy = 
    expand_grid(
      sp1 = selected_species,
      sp2 = selected_species
    ) %>%
    filter(sp2 >= sp1)
  
  sp_dtf = 
    xy %>%
    future_pmap_dfr(
      .f = function(sp1, sp2){
        df1 = 
          dat_filtered %>%
          filter(sp == sp1)
        
        df2 = 
          dat_filtered %>%
          filter(sp == sp2)
        
        nn12 = 
          dist.matrix(
            df1 %>% 
              select(gx, gy) %>%
              as.matrix,
            df2 %>% 
              select(gx, gy) %>%
              as.matrix,
            method = 'euclidean'
          ) %>%
          as.numeric
        
        nn12_positive = nn12[nn12 > 0]
        
        connection = mean(nn12_positive <= d_cutoff)
        
        return(tibble(sp1, sp2, connection))
      },
      .options = furrr_options(seed = TRUE)
    ) %>%
    mutate(
      sp1 = factor(sp1), 
      sp2 = factor(sp2)
    ) %>%
    arrange(sp1, sp2)
  
  res = 
    sp_dtf %>%
    bind_rows(
      sp_dtf %>% 
        filter(sp1 != sp2) %>%
        rename(sp1 = sp2, sp2 = sp1)
    ) %>%
    unique
  
  res_matrix =
    res %>% 
    pivot_wider(
      id_cols = sp1, 
      names_from = sp2, 
      values_from = connection,
      values_fill = NA
    ) %>%
    select(-1) %>% 
    as.matrix 
  
  adjacency_matrix = 1 * (res_matrix > cumulative_null_prob_threshold)
  
  return(
    list(
      dat = dat,
      parms = 
        tibble(
          distance_cutoff = d_cutoff,
          distance_step = d_step,
          Lx = Lx,
          Ly = Ly,
          abundance_cutoff = abundance_threshold
        ),
      abundances = abuns,
      binary_adjacency_matrix = adjacency_matrix,
      weighted_adjacency_matrix = res_matrix
    )
  )
}

cluster_analysis = 
  function(
    A, 
    algorithm,
    weighted,
    self_loops
  ){
    
    if(!weighted){
      adjacency_matrix = A$binary_adjacency_matrix
      weight_parameter = NULL
    } 
    
    if(weighted){
      adjacency_matrix = A$weighted_adjacency_matrix
      weight_parameter = TRUE
    } 
    
    graph = 
      graph_from_adjacency_matrix(
        adjacency_matrix, 
        weighted = weight_parameter,
        mode = 'undirected',
        diag = self_loops
      )  
    
    fn = get(paste0('cluster_', algorithm))
    communities = fn(graph)
    
    result = 
      membership(communities) %>%
      enframe %>%
      rename(group = value) %>%
      mutate(
        algorithm = algorithm, 
        weighted = weighted,
        self_loops = self_loops,
        group = factor(group),
        modularity = modularity(communities),
        number_of_groups = length(communities)
      ) 
    
    return(
      list(
        data = A,
        adjacency_matrix = adjacency_matrix,
        parms = 
          tibble(
            algorithm, 
            weighted , 
            self_loops
          ),
        graph = graph,
        communities = communities,
        result = result
      )
    )
    
  }
    
wrapper = 
  function(
    fdp,
    census,
    autolinked_species_only,
    weighted,
    self_loops,
    algorithm,
    d_cutoff,
    d_step,
    Lx,
    Ly
  ){
    if(fdp == 'bci'){
      dat = read_bci(census)
    }
    
    if(fdp == 'laplanada'){
      dat = read_laplanada(census)
    }
    
    A = adjacency_matrix(dat, autolinked_species_only, d_cutoff, d_step, Lx, Ly)
    result = cluster_analysis(A, algorithm, weighted, self_loops)$result
    
    result %>%
      mutate(
        fdp = fdp,
        census = census,
        autolinked_species_only = autolinked_species_only,
        d_cutoff = d_cutoff,
        d_step = d_step,
        Lx = Lx,
        Ly = Ly
      ) %>%
      return
    
  }

if(do.data){
  
  plan(multisession, workers = detectCores() - 1)
  
  parameters = 
    expand_grid(
      fdp = 'bci',
      census = 0:7,
      autolinked_species_only = c(TRUE, FALSE),
      weighted = c(TRUE, FALSE),
      algorithm = 'louvain',
      self_loops = c(TRUE, FALSE),
      d_cutoff = 5:20,
      d_step = 1e-5,
      Lx = 1000,
      Ly = 500
    )
  
  results = 
    parameters %>%
    future_pmap_dfr(
      wrapper,
      .options = furrr_options(seed = TRUE)
    )
  
  saveRDS(results, file = '~/SpatialNiche/Data/20210712/all_bci_censuses.rds')
  
}

if(do.plots.20210709){
  
  results = readRDS('~/SpatialNiche/Data/20210709/all_bci_censuses.rds')
  
  theme_set(theme_bw())
  theme_update(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  plot = 
    results %>% 
    group_by(census, method, d_cutoff) %>% 
    summarize(number_of_groups = length(unique(group))) %>% 
    ungroup() %>% 
    mutate(number_of_groups = factor(number_of_groups)) %>% 
    ggplot(aes(d_cutoff, census, fill = number_of_groups)) + 
    geom_tile(color = 'black') + 
    facet_wrap(~method) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    labs(fill = '# colors')
  
  plot_modularity = 
    results %>% 
    group_by(census, method, d_cutoff) %>% 
    summarize(modularity = unique(modularity)) %>% 
    ungroup() %>% 
    ggplot(aes(d_cutoff, census, fill = modularity)) + 
    geom_tile(color = 'black') + 
    facet_wrap(~method) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    scale_fill_gradientn(colors = terrain.colors(100))
  
  gridExtra::grid.arrange(plot, plot_modularity)
}


if(do.plots.20210712){
  
  results = readRDS('~/SpatialNiche/Data/20210712/all_bci_censuses.rds')
  
  results %<>%
    mutate(
      number_of_groups = ifelse(modularity > 0, number_of_groups, NA),
      modularity = ifelse(modularity > 0, modularity, NA)
    )
  
  theme_set(theme_bw())
  theme_update(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  plot_ngroups_unweighted = 
    results %>% 
    filter(weighted == FALSE) %>%
    select(-c(name, group)) %>%
    unique() %>%
    mutate(number_of_groups = factor(number_of_groups)) %>% 
    ggplot(aes(d_cutoff, census, fill = number_of_groups)) + 
    geom_tile(color = 'black') + 
    facet_grid(self_loops ~ autolinked_species_only, labeller = label_both) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    labs(fill = '# colors') +
    ggtitle('BCI - Binary adjacency matrix - Louvain algorithm')
  
  plot_ngroups_weighted = 
    results %>% 
    filter(weighted == TRUE) %>%
    select(-c(name, group)) %>%
    unique() %>%
    mutate(number_of_groups = factor(number_of_groups)) %>% 
    ggplot(aes(d_cutoff, census, fill = number_of_groups)) + 
    geom_tile(color = 'black') + 
    facet_grid(self_loops ~ autolinked_species_only, labeller = label_both) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    labs(fill = '# colors') +
    ggtitle('BCI - Quantitative adjacency matrix - Louvain algorithm')
  
  plot_modularity_unweighted = 
    results %>% 
    filter(weighted == FALSE) %>%
    select(-c(name, group)) %>%
    unique() %>%
    ggplot(aes(d_cutoff, census, fill = modularity)) + 
    geom_tile(color = 'black') + 
    facet_grid(self_loops ~ autolinked_species_only, labeller = label_both) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    scale_fill_gradientn(colors = terrain.colors(100)) +
    ggtitle('BCI - Binary adjacency matrix - Louvain algorithm')
  
  plot_modularity_weighted = 
    results %>% 
    filter(weighted == TRUE) %>%
    select(-c(name, group)) %>%
    unique() %>%
    ggplot(aes(d_cutoff, census, fill = modularity)) + 
    geom_tile(color = 'black') + 
    facet_grid(self_loops ~ autolinked_species_only, labeller = label_both) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    scale_fill_gradientn(colors = terrain.colors(100)) +
    ggtitle('BCI - Quantitative adjacency matrix - Louvain algorithm')
  
  gridExtra::grid.arrange(
    plot_ngroups_unweighted, 
    plot_ngroups_weighted,
    plot_modularity_unweighted,
    plot_modularity_weighted
  )
  
  plot_mod_vs_ngroups = 
    results %>%
    filter(
      weighted == FALSE, 
      autolinked_species_only == TRUE, 
      self_loops == FALSE
    ) %>%
    select(-c(name, group)) %>%
    unique() %>%
    ggplot(aes(number_of_groups, modularity)) +
    geom_jitter(height = 0, width = .1) +
    facet_wrap(~d_cutoff, labeller= label_both)
  
  plot_modularity_vs_dcutoff = 
    results %>% 
    filter(weighted == FALSE) %>%
    select(-c(name, group)) %>%
    unique() %>%
    mutate(census = factor(census)) %>%
    ggplot(aes(d_cutoff, modularity, group = census, color = census)) + 
    geom_line() +
    geom_point() +
    facet_grid(self_loops ~ autolinked_species_only, labeller = label_both) + 
    theme(aspect.ratio = 1) +
    xlab('distance cutoff') +
    ggtitle('BCI - Binary adjacency matrix - Louvain algorithm')
  
}

if(do.bci.census.analysis){
  data =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
      )
    )
  
  results = readRDS('~/SpatialNiche/Data/20210712/all_bci_censuses.rds')
  
  data %<>% 
    filter(sp %in% intersect(data$sp, results$name))
  
  census_turnover = 
    data %>%
    count(census) %>%
    mutate(
      mean = mean(n),
      turnover = n - mean(n)
    )
  
  census_anomaly =
    data %>%
    count(census, sp) %>%
    left_join(
      data %>%
        count(sp) %>%
        rename(mean = n) %>%
        mutate(mean = mean / 7),
      by = 'sp'
    ) %>%
    mutate(anomaly = 100*(n / mean - 1))
  
  plot_anomaly =
    census_anomaly %>%
    filter(mean > 40) %>%
    mutate(census = factor(census)) %>%
    ggplot(aes(census, abs(anomaly), fill = census)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = .1, alpha = .2) +
    geom_hline(yintercept = 0, color = rgb(153/255, 0, 0)) +
    theme(legend.position = 'none') +
    ggtitle('BCI -- 100 * (abundance_census / abundance_mean - 1)') +
    coord_cartesian(ylim = c(0, 75))
  
  summary = 
    data %>%
    count(census, sp) %>%
    group_by(sp) %>%
    summarize_at(
      'n', 
      list(
        mean = mean, 
        sd = sd,
        range = function(x) max(x) - min(x)
      )
    ) %>%
    ungroup() %>%
    mutate(cv = sd / mean)
  
  plot_summary = 
    summary %>%
    filter(mean > 40) %>%
    mutate(sp = factor(sp, levels = sp[order(mean)])) %>%
    ggplot(aes(mean, sp)) +
    geom_col(fill = 'grey80') +
    geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), color = 'red')
  
  plot_range = 
    summary %>%
    filter(mean > 40) %>%
    mutate(sp = factor(sp, levels = sp[order(mean)])) %>%
    ggplot(aes(range, sp)) +
    geom_col()
  
  gridExtra::grid.arrange(
    plot_summary + theme(aspect.ratio = 3), 
    plot_range + theme(aspect.ratio = 3),
    nrow = 1
  )
}
