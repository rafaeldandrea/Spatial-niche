library(wordspace) ## for function dist.matrix()
library(seriation) ## for function seriate()
library(distances) ## for functions distances(), nearest_neighbor_search()
library(janitor)   ## for function get_dupes()

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




# distfun = function(id1, id2){
#   xymap %>% 
#     filter(site %in% c(id1, id2)) %>% 
#     dist %>% 
#     as.numeric
# } 

distfun = function(id1, id2){
  x1 = xymap$x[id1]
  x2 = xymap$x[id2]
  y1 = xymap$y[id1]
  y2 = xymap$y[id2]
  
  sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2)
} 


number_of_distance_breaks = 50
number_of_dominant_species = 20
number_of_nearest_neighbors = 1

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

euclidean2d_distance_probability <- function(d, a){
    if(d > a * sqrt(2)) return(0)
    f =
      ifelse(
        d < a,
        -4 * d / a ^ 3 + pi / a ^ 2 + d ^ 2 / a ^ 4,
        -2 / a ^ 2 + 4 / a ^ 2 * asin(a / d) + 4 / a ^ 3 * sqrt(d ^ 2 - a ^ 2) - pi / a ^ 2 - d ^ 2 / a ^ 4
      )
    return(2 * d * f)
} 

distances = seq(0, sqrt(2), length = number_of_distance_breaks)
pdexp = sapply(distances, function(d) euclidean2d_distance_probability(d, 1))

ind = make_index(19)
setwd('~/SpatialNiche/Data/20200709-80scen-A20kS220-H-betanoise-biotime2e5/')
datafiles = paste0('scenario_', 1:80,'.RData')
datafile = datafiles[ind]
dat = get(load(datafile))

censuses = dat$evolving_community
linear_community_dimension = sqrt(length(censuses[[1]]))

xymap = 
  crossing(
    x = seq(linear_community_dimension) / linear_community_dimension, 
    y = seq(linear_community_dimension) / linear_community_dimension
  ) %>% 
  arrange(y) %>% 
  rowid_to_column(var = 'site')

dists = 
  xymap %>% 
  as.matrix %>%
  distances()

community_mode =
  crossing(
    x = seq(linear_community_dimension) / linear_community_dimension,
    y = seq(linear_community_dimension) / linear_community_dimension
  ) %>%
  arrange(y) %>%
  mutate(
    species =
      do.call(cbind, censuses) %>%
      apply(., 1, getmode)
  ) %>%
  mutate(dupe_count = 1) %>%
  rowid_to_column(var = 'site')

# community_mode =
#   crossing(
#     x = seq(linear_community_dimension) / linear_community_dimension,
#     y = seq(linear_community_dimension) / linear_community_dimension
#   ) %>%
#   arrange(y) %>%
#   rowid_to_column(var = 'site') %>%
#   bind_cols(
#     as_tibble(do.call(cbind, censuses[1:100]))
#   ) %>%
#   pivot_longer(-c(x, y, site), names_to = 'census', values_to = 'species') %>%
#   select(-census) %>% 
#   get_dupes(site, species) %>%
#   unique

abundance_restrictions = 
  community_mode %>% 
  group_by(species) %>% 
  count %>%
  arrange(desc(n)) %>%
  filter(n > number_of_nearest_neighbors)

# abundances = 
#   community_mode %>%
#   group_by(species) %>%
#   summarize(n = sum(dupe_count), .groups = 'drop') %>%
#   arrange(desc(n))

community_mode = 
  community_mode %>%
  filter(species %in% abundance_restrictions$species[1:20])


# community_mode = 
#   community_mode %>%
#   filter(species %in% abundances$species[1:20])

# community_mode =
#   crossing(
#     x = seq(linear_community_dimension) / linear_community_dimension,
#     y = seq(linear_community_dimension) / linear_community_dimension
#   ) %>%
#   arrange(y) %>%
#   mutate(
#     species = censuses[[length(censuses)]]
#   )

species_list = sort(unique(community_mode$species))
  
dominant_species = 
  plyr::count(community_mode$species) %>% 
  arrange(desc(freq)) %>% 
  as_tibble %>%
  pull(x)

number_of_dominant_species = min(number_of_dominant_species, length(dominant_species))
number_of_dominant_species = length(dominant_species)


res = 
  crossing(
    sp1 = species_list, 
    sp2 = species_list
  ) %>% 
  filter(sp2 >= sp1) %>%
  arrange(sp1) %>%
  ddply(
    .(sp1, sp2), 
    function(df){
      sp1 = df$sp1
      sp2 = df$sp2
      
      writeLines(paste(sp1, sp2))
      
      df1 = 
        community_mode %>%
        filter(species == sp1)
      
      df2 = 
        community_mode %>%
        filter(species == sp2)
      
      nn12 = 
        nearest_neighbor_search(
          dists, 
          k = number_of_nearest_neighbors, 
          query_indices = df1$site, 
          search_indices = df2$site
        )
      
      mean_nn_dist = 
        nn12 %>%
        as_tibble %>%
        mutate(k = 1:number_of_nearest_neighbors) %>% 
        pivot_longer(-k, names_to = 'treeID', values_to = 'neighbor') %>% 
        mutate(
          sp1 = sp1,
          sp2 = sp2
        ) %>%
        mutate(treeID = as.integer(treeID)) %>% 
        mutate(dist = distfun(treeID, neighbor)) %>%
        inner_join(
          df1 %>%
            rename(
              treeID = site,
              n1 = dupe_count
            ) %>%
            select(treeID, n1),
          by = 'treeID'
        ) %>%
        inner_join(
          df2 %>%
            rename(
              neighbor = site,
              n2 = dupe_count
            ) %>%
            select(neighbor, n2),
          by = 'neighbor'
        ) %>%
        group_by(sp1, sp2) %>%
        summarize(weighted_dist = weighted.mean(dist, n1 * n2), .groups= 'drop')
    }
  ) %>%
  as_tibble

res = 
  res %>%
  bind_rows(
    res %>%
      rename(tsp1 = sp2, tsp2 = sp1) %>% 
      rename(sp1 = tsp1, sp2 = tsp2)
  ) %>% 
  unique

res_matrix =
  res %>% 
  pivot_wider(
    id_cols = sp1, 
    names_from = sp2, 
    values_from = weighted_dist
  ) %>%
  select(-1) %>% 
  as.matrix 

seriated_order = 
  res_matrix %>% 
  seriate %>% 
  get_order
 

res_long = 
  res %>%
  mutate(
    sp1 = factor(sp1, levels = species_list[seriated_order]),
    sp2 = factor(sp2, levels = species_list[seriated_order])
  ) %>%
  arrange(sp1, sp2) %>%
  mutate(
    true_similarity = as.numeric((dat$C %*% t(dat$C))[seriated_order, seriated_order])
  )

plot0 = 
  crossing(
    x = seq(linear_community_dimension) / linear_community_dimension,
    y = seq(linear_community_dimension) / linear_community_dimension
  ) %>%
  arrange(y) %>%
  rowid_to_column(var = 'site') %>%
  bind_cols(tibble(species = censuses[[100]])) %>%
  ggplot(aes(x, y, fill = species)) +
  geom_tile()

plot1 = 
  res_long %>%
  ggplot(aes(sp1, sp2, fill = weighted_dist)) +
  geom_tile()

plot2 = 
  res_long %>%
  ggplot(aes(sp1, sp2, fill = true_similarity)) +
  geom_tile()

plot3 = 
  res_long %>%
  ggplot(aes(weighted_dist, true_similarity)) +
  geom_point() +
  geom_smooth(method = 'lm')

gridExtra::grid.arrange(
  plot0,
  plot1,
  plot2,
  plot3,
  nrow = 2
)


if(FALSE){
df1 = 
  community_mode %>% 
  filter(species == sp1) %>% 
  select(site, x, y) %>% 
  get_dupes(site) %>%
  unique

df2 = 
  community_mode %>%
  filter(species == sp2) %>% 
  select(site, x, y) %>% 
  get_dupes(site) %>%
  unique

multiplicity = outer(df1$dupe_count, df2$dupe_count)

nearest_neighbor_counts = 
  ddply(df1, .(x, y), function(v){
    x1 = v$x
    y1 = v$y
    n1 = v$dupe_count
    
    n2 = 
      df2 %>%
      mutate(ed = sqrt((x1 - x) ^ 2 + (y1 - y) ^ 2)) %>%
      filter(ed > 0 & ed <= 2 / linear_community_dimension) %>%
      pull(dupe_count) %>%
      sum
    
    return(n1 * n2)
  }) %>%
  as_tibble

observed_distance_dbn =
  dist.matrix(
    df1 %>% select(x, y) %>% unique %>% as.matrix,
    df2 %>% select(x, y) %>% unique %>% as.matrix,
    method = 'euclidean'
  ) %>%
  hist(., breaks = distances, probability = TRUE, plot = FALSE)
}

if(FALSE){
  similarity = NULL
  for(sp1 in dominant_species[1:number_of_dominant_species]){ 
    for(sp2 in dominant_species[1:number_of_dominant_species]){
      observed_distance_dbn = 
        community_mode %>% 
        filter(species == sp1) %>% 
        select(x, y) %>% 
        as.matrix() %>%
        dist.matrix(
          community_mode %>%
            filter(species == sp2) %>% 
            select(x, y) %>% 
            as.matrix(),
          method = 'euclidean'
        ) %>%
        hist(., breaks = distances, probability = TRUE, plot = FALSE)
      
      pdobs = c(0, observed_distance_dbn$density)
      
      similarity = c(similarity, pdobs[2] / pdexp[2])
      
      show.plots = FALSE
      if(show.plots){
        df = 
          tibble(
            distance = distances,
            expected = pdexp,
            observed = pdobs
          ) %>%
          filter(distance > 0)
        
        plot1 = 
          community_mode %>%
          filter(species %in% c(sp1, sp2)) %>%
          mutate(species = factor(species)) %>%
          ggplot(aes(x, y, fill = species)) +
          geom_tile()
        
        plot2 = 
          df %>%
          pivot_longer(-distance, names_to = 'distribution') %>%
          ggplot() +
          geom_line(aes(distance, value, group = distribution, color = distribution))
        
        plot3 =
          df %>%
          filter(distance <= 1) %>%
          ggplot() +
          geom_line(aes(distance, observed / expected - 1)) +
          geom_point(aes(distance, observed / expected - 1)) +
          geom_hline(aes(yintercept = 0), color = red)
        
        gridExtra::grid.arrange(
          plot1,
          plot2,
          plot3,
          nrow = 1
        )
      }
    }
  }
  
  similarity_matrix = matrix(similarity, number_of_dominant_species)
  seriated_order = get_order(seriate(similarity_matrix))
  blockdiag_similarity_matrix = similarity_matrix[seriated_order, seriated_order]
  
  scaled_community_matrix = blockdiag_similarity_matrix / max(blockdiag_similarity_matrix)
  seriated_species = dominant_species[1: number_of_dominant_species][seriated_order]
  colnames(scaled_community_matrix) = seriated_species
  rownames(scaled_community_matrix) = seriated_species
  
  scaled_community_matrix_dtf = 
    tibble(
      sp1 = rep(seriated_species, each = number_of_dominant_species), 
      sp2 = rep(seriated_species, times = number_of_dominant_species),
      estimated_similarity = as.numeric(scaled_community_matrix),
      true_similarity = as.numeric((dat$C %*% t(dat$C))[seriated_species, seriated_species])
    ) %>%
    mutate(true_similarity = true_similarity / max(true_similarity))
  
  plot_community_mode = 
    community_mode %>%
    filter(species %in% dominant_species) %>%
    ggplot(aes(x, y, fill = species)) +
    geom_tile()
  
  plot_scaled_community_matrix = 
    scaled_community_matrix_dtf %>%
    mutate(
      sp1 = factor(sp1, levels = unique(sp1)),
      sp2 = factor(sp2, levels = unique(sp2))
    ) %>%
    ggplot(aes(sp1, sp2, fill = estimated_similarity)) +
    geom_tile() +
    labs(
      x = 'species 1',
      y = 'species 2'
    )
  
  plot_true_similarity = 
    scaled_community_matrix_dtf %>%
    mutate(
      sp1 = factor(sp1, levels = unique(sp1)),
      sp2 = factor(sp2, levels = unique(sp2))
    ) %>%
    ggplot(aes(sp1, sp2, fill = true_similarity)) +
    geom_tile() +
    labs(
      x = 'species 1',
      y = 'species 2'
    )
  
  plot_correlation =
    scaled_community_matrix_dtf %>%
    ggplot(aes(estimated_similarity, true_similarity)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0)
  
  gridExtra::grid.arrange(
    plot_community_mode,
    plot_scaled_community_matrix,
    plot_true_similarity,
    plot_correlation,
    nrow = 2
  )
  
}
