library(tidyverse)
library(wordspace) ## for function dist.matrix()
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
                ## a graph then find communities in the graph
library(furrr)  ## for parallel computing
library(parallel)  

plan(multisession, workers = detectCores() - 1)

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


calculate_cutoffs = 
  function(distance_threshold, Lx, Ly, distance_step = 1e-5){
    
    distances_fdp = seq(0, distance_threshold, by = distance_step) 
    
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
    
    return(
      list(
        nstar = abundance_threshold, 
        pstar = cumulative_null_prob_threshold
      )
    )
  }


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

network_analysis =
  function(fdp, distance_threshold, seed){
    
    # fdp = scenario$fdp
    # distance_threshold = scenario$distance_threshold
    # seed = scenario$seed
    
    set.seed(seed)
    
    if(fdp == 'bci'){
      dat = bci 
      Lx = 1000
      Ly = 500
    } 
    
    if(fdp == 'lap'){
      dat = lap
      Lx = 500
      Ly = 500
    }
    
    cutoffs = calculate_cutoffs(distance_threshold, Lx, Ly)
    abundance_threshold = cutoffs$nstar
    cumulative_null_prob_threshold = cutoffs$pstar
    
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
    
    communities_list = 
      list(
        communities_walktrap = cluster_walktrap(graph),
        communities_spinglass = cluster_spinglass(graph),
        communities_louvain = cluster_louvain(graph)
      )
    
    number_of_clusters = 
      sapply(
        communities_list, 
        length
      )
    
    modularity = 
      sapply(
        communities_list, 
        modularity
      )
    
    return(
      tibble(
        fdp, 
        distance_threshold, 
        abundance_threshold, 
        algorithm = c('walktrap', 'spinglass', 'louvain'),
        number_of_clusters, 
        modularity
      )
    )
  }

scenarios = 
  expand_grid(
    fdp = c('bci', 'lap'),
    distance_threshold = 5:25,
    seed = 1
  )

res =
  scenarios %>%
  future_pmap_dfr(
    .f = network_analysis
  )

saveRDS(res, '~/SpatialNiche/Data/network_analysis_results.rds')
