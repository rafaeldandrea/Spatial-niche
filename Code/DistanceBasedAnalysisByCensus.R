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
library(ks)       ## for kernel density estimation function kde() 

plan(multisession, workers = detectCores() - 1)
load('~/R_Functions/colors.rdata')


run.analysis = 0
wrangle.results = 0
wrangle.confusion = 0
do.fdp = 1
do.network.validation = 0

fdp = 'bci'

## ================ Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 
## 0 if running on my PC.
hpcc = !(user <- Sys.info()['user']) %in% c('rafael', 'rdand')

## pdf of the distance d between two random points in a 
## square of side a, from Johan Philip 2007
euclidean2d_distance_probability <- function(d, a){
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
euclidean2d_distance_probability_ab <- function(d, a, b){
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
randomize_adjacency_matrix <- function(nvertices, nedges = 1) {
  edges.max = nvertices * (nvertices - 1) / 2
  
  stopifnot(length(nvertices) == 1 && 1 <= nvertices)
  stopifnot(0 <= nedges & nedges <= edges.max)
  
  index.edges = lapply(list(1:(nvertices - 1)), function(k) rep(k * (k + 1) / 2, nvertices - k)) 
  index.edges = index.edges[[1]] + 1:edges.max
  graph.adjacency = matrix(0, ncol = nvertices, nrow = nvertices)
  graph.adjacency[sample(index.edges, nedges)] = 1
  
  return(graph.adjacency + t(graph.adjacency) + diag(nvertices))
}

## numerical confirmation that euclidean2d_distance_probability(d, a)
## provides the correct probabilities. 
## Notice that the integration for the cumulative prob is over dr, not 2 * pi * r * dr
## I currently don't understand why.
if(FALSE){
  test.dists = 
    dist.matrix(
      matrix(runif(4000), 2000), 
      matrix(runif(1000), 500), 
      method = 'euclidean'
    )
  plot(sort(test.dists), seq_along(test.dists) / length(test.dists), las = 1, t='l')
  x = seq(0, sqrt(2), l = 1000)
  p = sapply(x, function(d) euclidean2d_distance_probability(d, 1))
  F = cumsum(c(0, diff(x)) * p)
  lines(x, F, col = colors$red)  
}

distances = 0:ceiling(sqrt(2) * 144)/144
pdexp = sapply(distances, function(d) euclidean2d_distance_probability(d, 1))

cumulative_null_prob = cumsum(pdexp * c(0, diff(distances)))

setwd('~/SpatialNiche/Data/20200709-80scen-A20kS220-H-betanoise-biotime2e5/')
datafiles = paste0('scenario_', 1:80,'.RData')

## Read array argument from command line
scenario = as.numeric(commandArgs(TRUE))[1]
if(is.na(scenario)) scenario = 19

datafile = datafiles[scenario]
dat = get(load(datafile))

if(run.analysis){
  censuses = dat$evolving_community
  linear_community_dimension = sqrt(length(censuses[[1]]))
  
  xymap = 
    tidyr::crossing(
      x = seq(linear_community_dimension) / linear_community_dimension, 
      y = seq(linear_community_dimension) / linear_community_dimension
    ) %>% 
    arrange(y) %>% 
    rowid_to_column(var = 'site')
  
  for(thecensus in 1:100){
    jobID = (scenario - 1) * 100 + thecensus
    
    filename = paste0('distance_based_analysis_job_', jobID, '.RData')
    
    try(
      {
        community =
          xymap %>%
          bind_cols(species = censuses[[thecensus]])
        
        species_list = sort(unique(community$species))
        # species_list = 1:10
  
        
        res = 
          tidyr::crossing(
            scenario = scenario,
            census = thecensus,
            sp1 = species_list, 
            sp2 = species_list
          ) %>% 
          filter(sp2 >= sp1) %>%
          ddply(
            .(scenario, census, sp1, sp2), 
            function(df){
              sp1 = df$sp1
              sp2 = df$sp2
              
              if(!hpcc) writeLines(paste(sp1, sp2))
              
              df1 = 
                community %>%
                filter(species == sp1)
              
              df2 = 
                community %>%
                filter(species == sp2)
              
              if(nrow(df1) * nrow(df2) * cumulative_null_prob[15] > 1){
                nn12 = 
                  dist.matrix(
                    df1 %>% 
                      select(x, y) %>%
                      as.matrix,
                    df2 %>% 
                      select(x, y) %>%
                      as.matrix,
                    method = 'euclidean'
                  ) %>%
                  as.numeric
                
                nn12_positive = nn12[nn12 > 0]
                connected = sapply(distances, function(d) mean(nn12_positive <= d)) - cumulative_null_prob
              }
            }
          ) %>%
          as_tibble 
          # rename(fraction_close_pairs = V1)
        
        if(hpcc) save(res, file = filename)
      },
      silent = TRUE
    )
    
  }
  
}

if(wrangle.results){
  
  do.plots = FALSE
 
  censuses = 1:100
  
  jobIDs = 100 * (scenario - 1) + censuses
  
  dtf = 
    tibble(jobID = jobIDs) %>%
    ddply(.(jobID), function(df){
      jobID = df$jobID
      readname = paste0('distance_based_analysis_job_', jobID, '.RData')
      if(file.exists(readname)) res = get(load(readname))
    }) %>% 
    as_tibble
  
  # res = 
  #   dtf %>%
  #   bind_rows(
  #     dtf %>%
  #       rename(sp1 = sp2, sp2 = sp1)
  #   ) %>% 
  #   unique %>%
  #   group_by(scenario, sp1, sp2) %>%
  #   summarize(
  #     fraction_close_pairs = mean(fraction_close_pairs), 
  #     .groups = 'drop'
  #  )
  
  res = 
    dtf %>%
    select(-census, -jobID) %>%
    group_by(scenario, sp1, sp2) %>%
    summarize_all(list(function(x){ mean(x, na.rm = TRUE)})) %>%
    ungroup
  
  res = 
    res %>%
    bind_rows(
      res %>% 
        filter(sp1 != sp2) %>%
        rename(sp1 = sp2, sp2 = sp1)
    ) %>%
    unique
  
  species_list = unique(res$sp1)
  
  # res_matrix =
  #   res %>% 
  #   pivot_wider(
  #     id_cols = sp1, 
  #     names_from = sp2, 
  #     values_from = fraction_close_pairs,
  #     values_fill = 0
  #   ) %>%
  #   select(-1) %>% 
  #   as.matrix 
  
  stats = NULL
  for(chosen_distance_in_pixels in 1:144){
    # chosen_distance_in_pixels = 1
    selected_column_from_res = 4 + chosen_distance_in_pixels
    
    if(!hpcc){
      res_matrix =
        res %>% 
        pivot_wider(
          id_cols = sp1, 
          names_from = sp2, 
          values_from = selected_column_from_res,
          values_fill = NA
        ) %>%
        select(-1) %>% 
        as.matrix 
    }
    
    if(hpcc){
      res_matrix =
        res %>% 
        pivot_wider(
          id_cols = sp1, 
          names_from = sp2, 
          values_from = selected_column_from_res,
          values_fill = list(NA)
        ) %>%
        select(-1) %>% 
        as.matrix 
    }
    
    
    seriated_order = 
      res_matrix %>% 
      seriate(method = 'PCA') %>% 
      get_order
    
    seriated_species = species_list[seriated_order]
    
    res_long = 
      res_matrix %>%
      as_tibble %>%
      cbind(sp1 = species_list) %>%
      pivot_longer(-sp1, names_to = 'sp2', values_to = 'clustering') %>%
      mutate(
        true_similarity = as.numeric((dat$C %*% t(dat$C))[species_list, species_list])
      ) %>%
      mutate(
        true_similarity_binary = round(true_similarity / max(true_similarity))
      ) %>%
      mutate(
        sp1 = factor(sp1, levels = seriated_species),
        sp2 = factor(sp2, levels = seriated_species)
      ) %>% 
      mutate(
        empirical_clustering = factor(1 * (clustering > 0), levels = c(0, 1)),
        reference_clustering = factor(true_similarity_binary, levels = c(0, 1))
      ) %>%
      mutate(
        clustering = clustering / cumulative_null_prob[1 + chosen_distance_in_pixels]
      )
    
    confusion = 
      caret::confusionMatrix(
        data = res_long$empirical_clustering, 
        reference = res_long$reference_clustering
      )
    
    accuracy = confusion$overall['Accuracy']
    
    if(do.plots){
      plot1 = 
        res_long %>% 
        ggplot(aes(sp1, sp2, fill = clustering)) + 
        geom_tile() +
        labs(
          x = 'species 1',
          y = 'species 2',
          fill = 'obs/exp pairs - 1'
        )
      
      plot2 = 
        res_long %>% 
        ggplot(aes(sp1, sp2, fill = empirical_clustering)) + 
        geom_tile() +
        labs(
          x = 'species 1',
          y = 'species 2',
          fill = 'inferred sim'
        )
      
      plot3 = 
        res_long %>%
        select(sp1, sp2, reference_clustering) %>%
        unique %>%
        ggplot(aes(sp1, sp2, fill = reference_clustering)) +
        geom_tile() +
        labs(
          x = 'species 1',
          y = 'species 2',
          fill = 'actual sim'
        )
      
      plot4 = 
        res_long %>%
        ggplot(aes(clustering, true_similarity_binary)) +
        geom_point() +
        geom_smooth(
          method = 'glm',
          formula = 'y ~ x',
          method.args = list(family = 'binomial'),
          se = FALSE
        ) +
        theme(aspect.ratio = 1) +
        labs(
          x = 'proportion of close pairs relative to expectations',
          y = 'actual similarity (rounded to 0 or 1)'
        )
      
      gridExtra::grid.arrange(
        plot1, 
        plot2, 
        plot3,
        plot4,
        nrow = 2,
        top = 
          paste(
            '  N Soil Types = ', c(5, 10, 20, 50)[1 + (scenario - 1) %/% 20], 
            '  Niche connected = ', round(1 - mean_cosine_rows(dat$C), 2),
            '  Sensitivity = ', round(confusion$byClass['Sensitivity'], 2),
            '  Specificity = ', round(confusion$byClass['Specificity'], 2),
            '  Balanced Accuracy = ', round(confusion$byClass['Balanced Accuracy'], 2)
          )
      )  
    }
    
    stats = 
      stats %>%
      rbind(
        tibble(
          scenario = scenario,
          n_soil_types = c(5, 10, 20, 50)[1 + (scenario - 1) %/% 20],
          niche_index = round(1 - mean_cosine_rows(dat$C), 2),
          threshold_distance_in_pixels = chosen_distance_in_pixels,
          sensitivity = confusion$byClass['Sensitivity'],
          specificity = confusion$byClass['Specificity'],
          balanced_accuracy = confusion$byClass['Balanced Accuracy']
        )  
      )  
  }
  
  if(hpcc) save(stats, file = paste0('distance_based_results_scenario_', scenario, '.RData'))
  
}

if(wrangle.confusion){
  filenames = paste0('distance_based_results_scenario_', 1:80, '.RData')
  bar = NULL
  for(char in filenames){
    try(
      {
        foo = get(load(char))
        bar = rbind(bar, foo)
      },
      silent = TRUE
    )
  }
  
  stats = 
    bar %>%
    mutate(
      nominal_niche_index = seq(.05, 1, by = .05)[1 + (scenario - 1) %% 20]
    )
  
  plot_stats = 
    stats %>%
    filter(
      niche_index > .5,
      threshold_distance_in_pixels <= 50
    ) %>%
    pivot_longer(
      c(sensitivity, specificity, balanced_accuracy),
      names_to = 'statistic',
      values_to = 'value'
    ) %>%
    ggplot(aes(threshold_distance_in_pixels * 100 / 144, value, group = statistic, color = statistic)) + 
    geom_line(size = 1) + 
    facet_wrap(n_soil_types ~ nominal_niche_index, ncol = 10) +
    theme(legend.position = 'top') +
    xlab('thresold distance relative to plot size (%)')
  
  plot_best_distance =
    stats %>%
    filter(niche_index > .5) %>%
    group_by(n_soil_types, nominal_niche_index) %>%
    slice_max(balanced_accuracy) %>%
    mutate(n_soil_types = as.factor(n_soil_types)) %>%
    ggplot(
      aes(
        nominal_niche_index, 
        threshold_distance_in_pixels * 100 / 144, 
        group = n_soil_types, 
        color = n_soil_types
      )
    ) +
    geom_line(size = 1) +
    geom_point()+
    labs(
      x = 'nominal niche connected',
      y = 'threshold distance relative to plot size (%)'
    )
  
  plot_peak_accuracy =
    stats %>%
    filter(niche_index > .5) %>%
    group_by(n_soil_types, nominal_niche_index) %>%
    slice_max(balanced_accuracy) %>%
    mutate(n_soil_types = as.factor(n_soil_types)) %>%
    ggplot(
      aes(
        nominal_niche_index, 
        balanced_accuracy, 
        group = n_soil_types, 
        color = n_soil_types
      )
    ) +
    geom_line(size = 1) +
    geom_point()+
    labs(
      x = 'nominal niche connected',
      y = 'balanced accuracy'
    )
  
  gridExtra::grid.arrange(
    plot_best_distance,
    plot_peak_accuracy,
    nrow = 1
  )
}

if(do.fdp){
  
  distance_threshold = 11 ## d* in meters
  distance_step = 1e-5 ## dr
  distances_fdp = seq(0, distance_threshold, by = distance_step) 
  
  Ly = 500
  if(fdp == 'bci') Lx = 1000
  if(fdp == 'lap') Lx = 500
  
  # p(r)
  pdexp_fdp = 
    sapply(
      distances_fdp, function(d) 
        euclidean2d_distance_probability_ab(d, a = Ly, b = Lx)
    )
  
  ## int p(r) * dr from 0 to d*
  cumulative_null_prob_threshold = sum(pdexp_fdp * distance_step) 
  
  ## abundance cutoff based on E[n] >= 1, where n is the number of pairs
  ## within a circle of radius d*
  abundance_threshold = round(sqrt(1 / cumulative_null_prob_threshold))
  
  if(FALSE){
    bci = NULL
    for(i in 1:7){
      census_file = paste0('~/bci/bci.full',i,'.rdata')
      bci = 
        bci %>%
        rbind(
          get(load(census_file)) %>%
            as_tibble %>%
            filter(dbh > 100) %>%
            filter(!is.na(gx)) %>%
            filter(!is.na(gy)) %>%
            mutate(census = i)
        )
    }
    
    abuns = 
      bci %>%
      group_by(sp) %>%
      count %>%
      arrange(desc(n)) %>%
      ungroup
    
    bci = 
      bci %>%
      select(treeID, sp, gx, gy, census, dbh) %>%
      unique %>%
      mutate(sp = ordered(factor(sp, levels = rev(abuns$sp))))
    
    save(bci, file = '~/SpatialNiche/Data/BCI/bci_coords_allcensuses.RData')
    
  }
  
  set.seed(scenario)
  
  lap = 
    read.table(
      '~/SpatialNiche/Data/relaplanadacensusdata/Censo 2_La Planada.txt', 
      header = TRUE
    ) %>%
    as_tibble %>%
    filter(dbh >= 100) %>%
    select(sp, gx, gy)
  
  bci = 
    get(
      load(
        url(
          'https://github.com/rafaeldandrea/BCI/blob/master/bci.full7.rdata?raw=true'
        )
      )
    ) %>%
    dplyr::filter(dbh >= 100) %>%
    select(sp, gx, gy) %>%
    unique() %>%
    as_tibble
  
  if(fdp == 'bci') dat = bci 
  if(fdp == 'lap') dat = lap
  
  
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
          mean(nn12_positive <= distance_threshold) -
          cumulative_null_prob_threshold
        
        return(tibble(sp = char, keep = connected > 0))
      }
    ) %>%
    filter(keep == TRUE) %>%
    pull(sp)
  
  xy = 
    expand_grid(
      sp1 = selected_species,
      sp2 = selected_species
    ) %>%
    filter(sp2 >= sp1)
  
  sp_dtf = 
    future_map2_dfr(
      .x = xy$sp1,
      .y = xy$sp2,
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
        
        connection = mean(nn12_positive <= distance_threshold)
        
        return(tibble(sp1, sp2, connection))
      }
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
  
  if(FALSE){
    seriated_order = 
      1 * (res_matrix > threshold_index) %>% 
      seriate(method = 'BEA_TSP', control = 'rep') %>% 
      get_order
    
    seriated_species = colnames(res_matrix)[seriated_order]
    
    res_seriated = 
      res %>%
      mutate(
        sp1 = factor(sp1, levels = seriated_species),
        sp2 = factor(sp2, levels = seriated_species)
      ) %>%
      arrange(sp1, sp2)
    
  }
  
  graph = 
    graph_from_adjacency_matrix(
      adjacency_matrix, 
      mode = 'undirected',
      diag = TRUE
    )
  
  graph_weighted = 
    graph_from_adjacency_matrix(
      res_matrix, 
      mode = 'undirected',
      weighted = TRUE,
      diag = TRUE
    )
  
  graph_no_self_loops = 
    graph_from_adjacency_matrix(
      adjacency_matrix, 
      mode = 'undirected', 
      diag = FALSE
    )
  
  # communities = cluster_walktrap(graph)
  # communities = cluster_spinglass(graph)
  communities = cluster_louvain(graph)
  
  if(fdp == 'lap'){
    groups = 
      with(
        communities, 
           1 * (membership == 3) + 
           2 * (membership == 2) + 
           3 * (membership == 1)
      )
    
    communities$membership = groups
  }
  
  membership = 
    membership(communities) %>%
    enframe %>%
    rename(group = value)
  
  modularity = modularity(communities)
  
  null_modularity = 
    future_map_dbl(1:1e3, function(k){
      
      null_adjacency_matrix =
        randomize_adjacency_matrix(
          nvertices = nrow(adjacency_matrix),
          nedges = sum(degree(graph))/ 2 - nrow(adjacency_matrix)
        )
      
      null_graph = 
        graph_from_adjacency_matrix(
          null_adjacency_matrix, 
          mode = 'undirected', 
          diag = TRUE
        )
      
      null_communities = cluster_louvain(null_graph)
      
      modularity(null_communities)
  })
  
  plot_null_modularity = 
    tibble(
      null = null_modularity, 
      data = modularity
    ) %>% 
    ggplot(aes(null)) + 
    geom_density(fill = 'grey') + 
    geom_vline(aes(xintercept = data), color = colors$red, size = 2) +
    labs(x = 'modularity') +
    ggtitle('Modularity') +
    theme(aspect.ratio = 1)
  
  adjacency_vs_communities_tibble =
    res %>%
    mutate(
      connected = 1 * (connection > cumulative_null_prob_threshold)
    ) %>%
    left_join(
      tibble(
        sp1 = membership$name,
        sp1_community = membership$group
      ),
      by = 'sp1'
    ) %>%
    left_join(
      tibble(
        sp2 = membership$name,
        sp2_community = membership$group
      ),
      by = 'sp2'
    ) %>%
    mutate(
      same_community = 1 * (sp1_community == sp2_community),
      consensus_community = as.factor(ifelse(sp1_community == sp2_community, sp1_community, 0))
    )
  
  membership_tibble = 
    adjacency_vs_communities_tibble %>%
    select(
      sp1, 
      sp1_community
    ) %>%
    unique %>%
    rename(
      sp = sp1,
      community = sp1_community
    ) %>%
    mutate(
      sp = as.character(sp),
      community = factor(community),
      nullcommunity = sample(community)
    )
  
  adjacency_vs_communities_for_plotting = 
    adjacency_vs_communities_tibble %>%
    mutate(
      sp1 = as.character(sp1), 
      sp2 = as.character(sp2)
    ) %>% 
    arrange(
      sp1_community, 
      sp2_community
    ) %>% 
    mutate(
      sp1 = factor(sp1, levels = unique(sp1)), 
      sp2 = factor(sp2, levels = unique(sp1)), 
      connected = factor(connected)
    ) 
  
  
  #### Plots #####
  if(!hpcc){
    
    theme_set(theme_bw())
    theme_update(panel.grid = element_blank())
    
    plot_all_data =
      bci %>% 
      filter(sp %in% abuns$sp[1:20]) %>%
      ggplot(aes(gx, gy, group = sp, color = sp)) + 
      geom_point() + 
      theme(legend.position = 'none')
    
    if(FALSE){
      plot_seriation_binary = 
        res_seriated %>%
        mutate(
          connected = ifelse(sp1 != sp2, 1, -1) * (connection > 0)
        ) %>%
        mutate(connected = factor(connected)) %>%
        ggplot(aes(sp1, sp2, fill = connected)) +
        geom_raster() +
        theme(aspect.ratio = 1) +
        labs(fill = 'connected') +
        scale_fill_manual(values = c('grey20', 'white', 'grey50')) +
        theme(legend.position = 'none') +
        theme(axis.text.x=element_blank())
      
      plot_seriation = 
        res_seriated %>%
        mutate(qtl = findInterval(connection, quantile(connection, 0:4/4))) %>%
        ggplot(aes(sp1, sp2, fill = qtl)) +
        geom_raster() +
        theme(aspect.ratio = 1)
      
      plot_seriation_vs_communities =
        adjacency_vs_communities_tibble %>%
        ggplot(aes(sp1, sp2, fill = consensus_community)) +
        geom_raster() +
        theme(aspect.ratio = 1) +
        theme(axis.text.x=element_blank()) +
        theme(legend.position = 'none') +
        scale_fill_manual(values = c('white', colors$red, colors$green, colors$blue))
      
    }
    
    plot_communities = 
      as.ggplot(
        ~plot(
          communities, 
          graph_no_self_loops, 
          col = with(colors, c(red, green, blue, yellow))[communities$membership]
        ) 
      ) +
      theme(aspect.ratio = 1) +
      theme(plot.margin = margin(-2, -3, -3, -3, 'cm'))
    
    plot_adjacency_vs_communities =
      adjacency_vs_communities_for_plotting %>% 
      ggplot(aes(sp1, sp2, fill = connected)) + 
      geom_raster() +
      theme(aspect.ratio = 1) +
      theme(legend.position = 'none') +
      scale_fill_manual(values = c('white', 'grey30')) +
      theme(
        axis.text.x = 
          element_text(
            angle = 90,
            vjust = 0.2,
            face = 'bold',
            color = 
              with(colors, c(red, green, blue, yellow))[
                adjacency_vs_communities_for_plotting %>% 
                  select(sp2, sp2_community) %>% 
                  unique() %>%
                  pull(sp2_community)
              ]
          )
      ) +
      theme(
        axis.text.y = 
          element_text(
            face = 'bold',
            color = 
              with(colors, c(red, green, blue, yellow))[
                adjacency_vs_communities_for_plotting %>% 
                select(sp2, sp2_community) %>% 
                unique() %>%
                pull(sp2_community)
              ]
          )
      )
    
    plot_fdp = 
      dat %>%
      mutate(sp = as.character(sp)) %>%
      inner_join(membership_tibble, by = 'sp') %>%
      ggplot(aes(gx, gy, color = community, group = community)) +
      geom_point() +
      labs(color = 'cluster') +
      theme(aspect.ratio = ifelse(fdp == 'bci', .5, 1)) +
      scale_color_manual(values = with(colors,c(red, green, blue, yellow)))
    
    plot_fdp_null = 
      dat %>%
      mutate(sp = as.character(sp)) %>%
      inner_join(membership_tibble, by = 'sp') %>%
      ggplot(aes(gx, gy, color = nullcommunity, group = nullcommunity)) +
      geom_point() +
      labs(color = 'cluster') +
      theme(aspect.ratio = ifelse(fdp == 'bci', .5, 1)) +
      scale_color_manual(values = with(colors,c(red, green, blue, yellow)))
    
    gridExtra::grid.arrange(
      plot_adjacency_vs_communities,
      plot_communities,
      plot_fdp,
      layout_matrix = rbind(c(1, 2), c(3, 3)),
      top = 
        paste(
          'All censuses (1982 - 2010)        d* = ',
          distance_threshold, 
          'm        Nmin = ', 
          abundance_threshold,
          '       No. species = ',
          length(selected_species),
          '       prop. trees > 10 cm dbh analyzed = ',
          round(nrow(dat_filtered) / nrow(bci), 2),
          '       modularity = ',
          round(modularity, 2),
          '       cor[adjacent, same community] = ',
          round(with(adjacency_vs_communities_tibble, cor(connected, same_community)), 2)
        )
    )
    
    
    print(paste('modularity = ', round(modularity, 2)))
    
    if(FALSE){
      gridExtra::grid.arrange(
        plot_bci,
        plot_bci_null,
        nrow = 2,
        top = 
          paste(
            'All censuses        d* = ',
            distance_threshold, 
            'm        Nmin = ', 
            abundance_threshold,
            '       No. species = ',
            length(selected_species),
            '       modularity = ',
            round(modularity, 2)
          )
      )  
    }
    
    
  }
  
  ## Kernel density estimation
  
  if(fdp == 'bci') nx = 50
  if(fdp == 'lap') nx = 25
  ny = 25
  
  evalpoints = 
    expand_grid(
      x = (1:nx - .5) * 20,
      y = (1:ny - .5) * 20
    )
  
  df = 
    dat |> 
    inner_join(
      membership |> 
        rename(sp = name, soiltype = group) |> 
        mutate(soiltype = as.factor(soiltype)), 
      by = 'sp')
  
  results =
    membership %>%
    pull(group) %>%
    unique %>% 
    sort %>%
    map_dfr(
      .f = function(group){
        x = 
          df |> 
          filter(soiltype == group) |> 
          select(gx, gy) |> 
          kde(
            eval.points = evalpoints,
            bgridsize = c(Lx / 20, Ly / 20), 
            xmin = c(0, 0), 
            xmax = c(1000, 500)
          ) 
        
        x$estimate %>% 
          as_tibble %>% 
          bind_cols(evalpoints) %>%
          pivot_wider(names_from = y, values_from = value) %>%
          pivot_longer(-x, names_to = 'y') %>% 
          mutate(
            y = str_remove(y, 'V') %>% 
              as.numeric,
            value = value / sum(value),
            soiltype = group
          ) %>%
          return
      }
    ) %>%
    group_by(x, y) %>%
    mutate(
      soiltype = 
        factor(
          soiltype, 
          levels = 
            membership %>%
            pull(group) %>%
            unique %>% 
            sort
        ),
      prob = value / sum(value)
    ) %>%
    ungroup()
  
  if(fdp == 'bci') asprat = .5
  if(fdp == 'lap') asprat = 1
  
  if(fdp == 'bci' & distance_threshold == 11){
    results %<>%
      mutate(
        soiltype = 
          1 * (soiltype == 1) +
          3 * (soiltype == 2) +
          2 * (soiltype == 3) +
          4 * (soiltype == 4)
      ) %>%
      mutate(soiltype = factor(soiltype, levels = 1:4))
  }
  
  plot_rasters =
    results %>% 
    ggplot(aes(x, y, fill = prob)) + 
    geom_raster() + 
    facet_wrap(~soiltype) +
    theme(aspect.ratio = .5) + 
    scale_fill_gradientn(colors = terrain.colors(100)) +
    theme(strip.background = element_rect(fill = 'orange')) +
    theme(aspect.ratio = asprat)
  
  plot_majority_vote = 
    results %>% 
    group_by(x, y) %>% 
    slice_max(prob) %>% 
    ggplot(aes(x, y, fill = soiltype)) + 
    geom_raster() + 
    theme(aspect.ratio = asprat)
  
  gridExtra::grid.arrange(
    plot_majority_vote,
    plot_rasters,
    ncol = 1
  )
  
}

if(do.network.validation){
  
{
  scenario = make_index(35)
  
  datafile = datafiles[scenario]
  dat = get(load(datafile))
  
  mean_abuns = 
    lapply(dat$evolving_community[94:100], function(sp) as_tibble(table(sp))) %>%
    bind_rows %>% 
    group_by(sp) %>% 
    summarize(n = mean(n), .groups = 'drop') %>% 
    arrange(desc(n))
  
  selected_species = as.numeric(gtools::mixedsort(mean_abuns$sp[1:45]))
  
  rownames(dat$C) = 1:220
  A = dat$C %*% t(dat$C)
  B = A[selected_species, selected_species]
  adjacency_matrix = round((B - min(B[B > 0]) / (max(B) - min(B[B > 0]))))
  
  if(FALSE){
    
    chosen_distance_in_pixels = 4
    
    censuses = 94:100
    
    jobIDs = 100 * (scenario - 1) + censuses
    
    selected_species = 
      lapply(dat$evolving_community[94:100], function(sp) as_tibble(table(sp))) %>%
      bind_rows %>% 
      group_by(sp) %>% 
      summarize(n = mean(n), .groups = 'drop') %>% 
      filter(n > sqrt(1 / cumulative_null_prob[chosen_distance_in_pixels])) %>%
      pull(sp)
    
    
    dtf = 
      tibble(jobID = jobIDs) %>%
      ddply(.(jobID), function(df){
        jobID = df$jobID
        readname = paste0('distance_based_analysis_job_', jobID, '.RData')
        if(file.exists(readname)) res = get(load(readname))
      }) %>% 
      as_tibble
    
    res = 
      dtf %>%
      mutate(
        sp1 = as.character(sp1),
        sp2 = as.character(sp2)
      ) %>%
      select(-census, -jobID) %>%
      group_by(scenario, sp1, sp2) %>%
      summarize_all(list(function(x){ mean(x, na.rm = TRUE)})) %>%
      ungroup
    
    res = 
      res %>%
      bind_rows(
        res %>% 
          filter(sp1 != sp2) %>%
          rename(sp1 = sp2, sp2 = sp1)
      ) %>%
      unique
    
    species_list = unique(res$sp1)
    
    stats = NULL
    # for(chosen_distance_in_pixels in 1:144){
    selected_column_from_res = 4 + chosen_distance_in_pixels
    
    
    adjacency_matrix =
      res %>% 
      mutate(
        sp1 = factor(sp1, levels = as.character(1:220)),
        sp2 = factor(sp2, levels = as.character(1:220))
      ) %>%
      pivot_wider(
        id_cols = sp1, 
        names_from = sp2, 
        values_from = selected_column_from_res,
        values_fill = 0
      ) %>%
      arrange(sp1) %>%
      select(gtools::mixedsort(species_list)) %>%
      as.matrix 
    
    adjacency_matrix = 
      adjacency_matrix[
        colnames(adjacency_matrix) %in% selected_species, 
        colnames(adjacency_matrix) %in% selected_species
      ]
    
    adjacency_matrix = 1 * (adjacency_matrix > 0)  
  }
  
  graph = 
    graph_from_adjacency_matrix(
      adjacency_matrix, 
      mode = 'undirected', 
      diag = TRUE
    )
  
  graph_no_self_loops = 
    graph_from_adjacency_matrix(
      adjacency_matrix, 
      mode = 'undirected', 
      diag = FALSE
    )
  
  communities = cluster_walktrap(graph)
  
  modularity = modularity(communities)
  membership = membership(communities)
  
  true_plotting_order = 
    tibble(
      sp = selected_species, 
      group = 1 + (sp - 1) %% ncol(dat$C)
    ) %>% 
    arrange(group)
  
  inferred_plotting_order = 
    tibble(sp = names(membership), cluster = membership) %>% 
    arrange(cluster)
  
  x = 
    A %>% 
    as_tibble %>% 
    mutate(sp1 = colnames(A)) %>% 
    pivot_longer(cols = -sp1, names_to = 'sp2', values_to = 'index') %>% 
    filter(sp1 %in% as.character(true_plotting_order$sp)) %>% 
    filter(sp2 %in% as.character(true_plotting_order$sp)) %>% 
    mutate(
      sp1 = factor(sp1, levels = true_plotting_order$sp), 
      sp2 = factor(sp2, levels = true_plotting_order$sp)
    ) %>% 
    ggplot(aes(sp1, sp2, fill = index)) + 
    geom_raster() +
    theme(aspect.ratio = 1)
  
  y = 
    A %>% 
    as_tibble %>% 
    mutate(sp1 = colnames(A)) %>% 
    pivot_longer(cols = -sp1, names_to = 'sp2', values_to = 'index') %>% 
    filter(sp1 %in% as.character(inferred_plotting_order$sp)) %>% 
    filter(sp2 %in% as.character(inferred_plotting_order$sp)) %>% 
    mutate(
      sp1 = factor(sp1, levels = inferred_plotting_order$sp), 
      sp2 = factor(sp2, levels = inferred_plotting_order$sp)
    ) %>% 
    ggplot(aes(sp1, sp2, fill = index)) + 
    geom_raster() +
    theme(aspect.ratio = 1)
  
  xx = 
    adjacency_matrix %>% 
    as_tibble %>% 
    mutate(sp1 = colnames(adjacency_matrix)) %>% 
    pivot_longer(cols = -sp1, names_to = 'sp2', values_to = 'index') %>% 
    filter(sp1 %in% as.character(true_plotting_order$sp)) %>% 
    filter(sp2 %in% as.character(true_plotting_order$sp)) %>% 
    mutate(
      sp1 = factor(sp1, levels = true_plotting_order$sp), 
      sp2 = factor(sp2, levels = true_plotting_order$sp)
    ) %>% 
    ggplot(aes(sp1, sp2, fill = index)) + 
    geom_raster() +
    theme(aspect.ratio = 1)
  
  yy = 
    adjacency_matrix %>% 
    as_tibble %>% 
    mutate(sp1 = colnames(adjacency_matrix)) %>% 
    pivot_longer(cols = -sp1, names_to = 'sp2', values_to = 'index') %>% 
    filter(sp1 %in% as.character(inferred_plotting_order$sp)) %>% 
    filter(sp2 %in% as.character(inferred_plotting_order$sp)) %>% 
    mutate(
      sp1 = factor(sp1, levels = inferred_plotting_order$sp), 
      sp2 = factor(sp2, levels = inferred_plotting_order$sp)
    ) %>% 
    ggplot(aes(sp1, sp2, fill = index)) + 
    geom_raster() +
    theme(aspect.ratio = 1)
  
  gridExtra::grid.arrange(x, y, xx, yy, nrow = 2)
  
  plot_communities = 
    as.ggplot(
      ~plot(
        communities, 
        graph_no_self_loops, 
        col = c(colors$red, colors$green, colors$blue)[communities$membership]
      ) 
    ) +
    theme(aspect.ratio = 1) +
    ggtitle(
      paste(
        'No. spp analyzed = ', length(selected_species),
        '    nominal niche index = ', seq(0, 1, l = 21)[1 + (scenario - 1) %% 20],
        '    No. soil types = ', ncol(dat$C),
        '    No. groups found = ', max(membership),
        '    modularity = ', round(modularity, 2)
      )
    )
  
  null_modularity = NULL
  for(null in 1:1e3){
    null_adjacency_matrix =
      randomize_adjacency_matrix(
        nvertices = nrow(adjacency_matrix),
        nedges = (sum(adjacency_matrix) - psych::tr(adjacency_matrix)) / 2
      )
    
    null_graph = 
      graph_from_adjacency_matrix(
        null_adjacency_matrix, 
        mode = 'undirected', 
        diag = TRUE
      )
    
    null_communities = cluster_walktrap(null_graph)
    
    null_modularity = c(null_modularity, modularity(null_communities))
  }
  
  gridExtra::grid.arrange(plot_communities)
  
  print(table(membership))
  
  print(mean(null_modularity > modularity))  
    
}
    
    ##############################
    
    confusion = 
      caret::confusionMatrix(
        data = res_long$empirical_clustering, 
        reference = res_long$reference_clustering
      )
    
    accuracy = confusion$overall['Accuracy']
    
    if(do.plots){
      plot1 = 
        res_long %>% 
        ggplot(aes(sp1, sp2, fill = clustering)) + 
        geom_tile() +
        labs(
          x = 'species 1',
          y = 'species 2',
          fill = 'obs/exp pairs - 1'
        )
      
      plot2 = 
        res_long %>% 
        ggplot(aes(sp1, sp2, fill = empirical_clustering)) + 
        geom_tile() +
        labs(
          x = 'species 1',
          y = 'species 2',
          fill = 'inferred sim'
        )
      
      plot3 = 
        res_long %>%
        select(sp1, sp2, reference_clustering) %>%
        unique %>%
        ggplot(aes(sp1, sp2, fill = reference_clustering)) +
        geom_tile() +
        labs(
          x = 'species 1',
          y = 'species 2',
          fill = 'actual sim'
        )
      
      plot4 = 
        res_long %>%
        ggplot(aes(clustering, true_similarity_binary)) +
        geom_point() +
        geom_smooth(
          method = 'glm',
          formula = 'y ~ x',
          method.args = list(family = 'binomial'),
          se = FALSE
        ) +
        theme(aspect.ratio = 1) +
        labs(
          x = 'proportion of close pairs relative to expectations',
          y = 'actual similarity (rounded to 0 or 1)'
        )
      
      gridExtra::grid.arrange(
        plot1, 
        plot2, 
        plot3,
        plot4,
        nrow = 2,
        top = 
          paste(
            '  N Soil Types = ', c(5, 10, 20, 50)[1 + (scenario - 1) %/% 20], 
            '  Niche connected = ', round(1 - mean_cosine_rows(dat$C), 2),
            '  Sensitivity = ', round(confusion$byClass['Sensitivity'], 2),
            '  Specificity = ', round(confusion$byClass['Specificity'], 2),
            '  Balanced Accuracy = ', round(confusion$byClass['Balanced Accuracy'], 2)
          )
      )  
    }
    
    stats = 
      stats %>%
      rbind(
        tibble(
          scenario = scenario,
          n_soil_types = c(5, 10, 20, 50)[1 + (scenario - 1) %/% 20],
          niche_index = round(1 - mean_cosine_rows(dat$C), 2),
          threshold_distance_in_pixels = chosen_distance_in_pixels,
          sensitivity = confusion$byClass['Sensitivity'],
          specificity = confusion$byClass['Specificity'],
          balanced_accuracy = confusion$byClass['Balanced Accuracy']
        )  
      )  
  }
  
  
 # }
