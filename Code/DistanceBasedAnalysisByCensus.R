library(wordspace) ## for function dist.matrix()
library(seriation) ## for function seriate()
library(caret) ## for function confusionMatrix()
library(ggplotify) ## to convert from base plot to grob
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
                ## a graph then find communities in the graph
library(tidyverse)


run.analysis = 0
wrangle.results = 0
wrangle.confusion = 0
do.bci = 1
do.network.validation = 0

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

## ================ Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc <- !(user <- Sys.info()['user']) %in% c('rafael', 'rdand')

## pdf of the distance d between two random points in a square of side a, from Johan Philip 2007
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

## pdf of the distance d between two random points in a square of sides a and b > a, from Johan Philip 2007
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
  lines(x, F, col = red)  
}

distances = 0:ceiling(sqrt(2) * 144)/144
pdexp = sapply(distances, function(d) euclidean2d_distance_probability(d, 1))

cumulative_null_prob = cumsum(pdexp * c(0, diff(distances)))

setwd('~/SpatialNiche/Data/20200709-80scen-A20kS220-H-betanoise-biotime2e5/')
datafiles = paste0('scenario_', 1:80,'.RData')

scenario = make_index(19)

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

if(do.bci){
  
  distance_threshold = 10 ## in meters
  abundance_threshold = 200
  
  distances_bci = seq(0, distance_threshold, by = 1e-5)
  
  pdexp_bci = 
    sapply(
      distances_bci, function(d) 
        euclidean2d_distance_probability_ab(d, a = 500, b = 1000)
    )
  
  cumulative_null_prob_threshold = sum(pdexp_bci) * 1e-5
  
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
  
  seed = make_index(1)
  set.seed(seed)
  
  # bci = 
  #   get(load('~/SpatialNiche/Data/BCI/bci_coords_allcensuses.RData')) %>%
  #   select(sp, gx, gy) %>%
  #   unique
  
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
  
  abuns = 
    bci %>%
    group_by(sp) %>%
    count %>%
    arrange(desc(n)) %>%
    ungroup
  
  bci_filtered = 
    bci %>%
    inner_join(abuns, by = 'sp') %>%
    filter(n >= abundance_threshold)
    
  selected_species = sort(unique(as.character(bci_filtered$sp)))
  
  selected_species =
    {
      bci_filtered %>%
        ddply(
          .(sp),
          function(df){
            d = dist(cbind(df$gx, df$gy))
            nn12_positive = d[d > 0]
            connected =
              mean(nn12_positive <= distance_threshold) -
              cumulative_null_prob_threshold
          }
        ) %>%
        filter(V1 > 0) %>%
        mutate(sp = as.character(sp)) %>%
        pull(sp)
    }

  sp_dtf =
    {
      tidyr::crossing(
        sp1 = selected_species,
        sp2 = selected_species
      ) %>%
        filter(sp2 >= sp1) %>%
        ddply(
          .(sp1, sp2), 
          function(df){
            sp1 = df$sp1
            sp2 = df$sp2
            
            # if(!hpcc) writeLines(paste(sp1, sp2))
            
            df1 = 
              bci_filtered %>%
              filter(sp == sp1)
            
            df2 = 
              bci_filtered %>%
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
              connected = 
                mean(nn12_positive <= distance_threshold) - 
                cumulative_null_prob_threshold
          }
        ) %>%
        as_tibble %>%
        rename(pairs_within_radius = V1)
      
    }
  
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
      values_from = pairs_within_radius,
      values_fill = NA
    ) %>%
    select(-1) %>% 
    as.matrix 
  
  threshold_index = 0
  
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
  
  adjacency_matrix = 1 * (res_matrix > 0)
  
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

  adjacency_vs_communities_tibble =
    res %>%
    mutate(
      connected = 1 * (pairs_within_radius > 0)
    ) %>%
    left_join(
      tibble(
        sp1 = names(membership),
        sp1_community = membership
      ),
      by = 'sp1'
    ) %>%
    left_join(
      tibble(
        sp2 = names(membership),
        sp2_community = membership
      ),
      by = 'sp2'
    ) %>%
    mutate(
      same_community = 1 * (sp1_community == sp2_community),
      consensus_community = as.factor(ifelse(sp1_community == sp2_community, sp1_community, 0))
    ) # %>%
    # mutate(
    #   sp1 = factor(sp1, levels = seriated_species),
    #   sp2 = factor(sp2, levels = seriated_species)
    # )
  
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
          connected = ifelse(sp1 != sp2, 1, -1) * (pairs_within_radius > 0)
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
        mutate(qtl = findInterval(pairs_within_radius, quantile(pairs_within_radius, 0:4/4))) %>%
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
        scale_fill_manual(values = c('white', red, green, blue))
      
    }
    
    plot_communities = 
      as.ggplot(
        ~plot(
          communities, 
          graph_no_self_loops, 
          col = c(red, green, blue)[communities$membership]
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
              c(red, green, blue)[
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
              c(red, green, blue)[
                adjacency_vs_communities_for_plotting %>% 
                select(sp2, sp2_community) %>% 
                unique() %>%
                pull(sp2_community)
              ]
          )
      )
    
    plot_bci = 
      bci %>%
      mutate(sp = as.character(sp)) %>%
      inner_join(membership_tibble, by = 'sp') %>%
      ggplot(aes(gx, gy, color = community, group = community)) +
      geom_point() +
      labs(color = 'cluster') +
      theme(aspect.ratio = .5) +
      scale_color_manual(values = c(red, green, blue))
    
    plot_bci_null = 
      bci %>%
      mutate(sp = as.character(sp)) %>%
      inner_join(membership_tibble, by = 'sp') %>%
      ggplot(aes(gx, gy, color = nullcommunity, group = nullcommunity)) +
      geom_point() +
      labs(color = 'cluster') +
      theme(aspect.ratio = .5) +
      scale_color_manual(values = c(red, green, blue))
    
    gridExtra::grid.arrange(
      plot_adjacency_vs_communities,
      plot_communities,
      plot_bci,
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
          round(nrow(bci_filtered) / nrow(bci), 2),
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
        col = c(red, green, blue)[communities$membership]
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
