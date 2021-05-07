#### ================= My OWN IMPLEMENTATION OF SDBW ======================
data = x
numclusters = 10
cluster_vector = sample(numclusters, size = nrow(data), replace = TRUE)
thedistance = 'euclidean'

SDbw <- function(data, numclusters, cluster_vector, metric = 'euclidean'){
  df <-
    data %>%
    as_tibble %>%
    rowid_to_column(var = 'tuple') %>%
    mutate(cluster_id = cluster_vector) %>%
    pivot_longer(cols = -c(tuple, cluster_id), names_to = 'feature')
  
  Var = function(x) sum(x ^ 2 - mean(x) ^ 2) / length(x)
  
  sigma_X <- apply(data, 2, Var)
  
  norm_sigma_X <- as.numeric(sqrt(sigma_X %*% sigma_X))
  
  sigma_Ci <-
    df %>%
    group_by(cluster_id, feature) %>%
    summarize(sigma = Var(value))
  
  norm_sigma_Ci <-
    sigma_Ci %>%
    group_by(cluster_id) %>%
    summarize(norm = sqrt(sigma %*% sigma))
  
  stdev <- 1/ numclusters * sqrt(sum(norm_sigma_Ci$norm))
  
  Scat <- 1 / numclusters * sum(norm_sigma_Ci$norm) / norm_sigma_X
  
  centroids <-
    df %>%
    group_by(cluster_id, feature) %>%
    summarize(coord = mean(value)) %>%
    ungroup %>%
    pivot_wider(names_from = feature, values_from = coord)
  
  medoids <- 
    df %>% 
    filter(tuple %in% GDAtools::medoids(dist(data), cl = cluster_vector)) %>%
    pivot_wider(names_from = feature, values_from = value) %>%
    select(-tuple)
  
  distance_to_centroid <- ## dims: centroids x tuples
    t(apply(centroids[, -1], 1, function(vector) 
      apply(data, 1, function(x) 
        dist(rbind(vector, x))
      )
    )) 
  
  distance_to_medoid <- ## dims: centroids x tuples
    t(apply(medoids[, -1], 1, function(vector) 
      apply(data, 1, function(x) 
        dist(rbind(vector, x))
      )
    )) 
  
  centroid_density = apply(distance_to_centroid, 1, function(d) sum(d <= stdev))
  
  medoid_density = apply(distance_to_medoid, 1, function(d) sum(d <= stdev))
  
  maxdensity_centroidpairs = outer(centroid_density, centroid_density, function(x, y) pmax(x, y))
  
  maxdensity_medoidpairs = outer(medoid_density, medoid_density, function(x, y) pmax(x, y))
  
  midpoints = NULL 
  for(i in seq(numclusters)) 
    for(j in seq(numclusters)) 
      midpoints = rbind(midpoints, 1 / 2 * (centroids[i, -1] + centroids[j, -1]))
  
  distance_to_midpoint <- ## dims: midpoints x tuples
    t(apply(midpoints, 1, function(vector) 
      apply(data, 1, function(x) 
        dist(rbind(vector, x))
      )
    ))
  
  belongs = matrix(nrow = numclusters ^ 2, ncol = nrow(data))
  for(clusterpair in seq(numclusters ^ 2)) 
    for(tuple in seq(nrow(data))){
      coords = pos2coord(pos = clusterpair, dim.mat = c(numclusters , numclusters))
      belongs[clusterpair, tuple] = 1 * (df$cluster_id[tuple] %in% coords)
    } 
  
  density_midpoints = matrix(apply(distance_to_midpoint * belongs, 1, function(d) sum(d > 0 & d <= stdev)), numclusters, numclusters)
  
  # dens_ratio <- 
  #   if(all(maxdensity_centroidpairs > 0)){
  #     density_midpoints / maxdensity_centroidpairs 
  #   } else{
  #     density_midpoints / maxdensity_medoidpairs 
  #   }
  
  dens_ratio = density_midpoints / maxdensity_medoidpairs 
  
  Dens_bw = 1 / numclusters / (numclusters - 1) * sum(dens_ratio * upper.tri(dens_ratio))
  
  S_Dbw = Scat + Dens_bw
  
  return(S_Dbw)
}


## ============= TEST ON ARTIFICIAL DATA ===============

set.seed(0)

minNST = 2
maxNST = 20

x <-
  rbind(
    matrix(rnorm(100, 0, 1), ncol = 2), 
    matrix(rnorm(100, 8, 1), ncol = 2), 
    matrix(rnorm(100, 13, .5), ncol = 2), 
    matrix(rnorm(100, 17, 1), ncol = 2),
    matrix(rnorm(100, 21, .5), ncol = 2)
  );

# x = rbind(x, cbind(runif(20, min(x[, 1]), max(x[, 1])), runif(20, min(x[, 2]), max(x[, 2]))))

dis_sites <- vegdist(x, method = "euclidean")

clus_sites <- hclust(dis_sites, method = 'ward.D2')

SD <-
  sapply(minNST:maxNST, function(numclust){
    cluster_index <- cutree(clus_sites, numclust)
    scatt <- clv.Scatt(x, cluster_index)
    dis <- clv.Dis(scatt$cluster.center)
    #dens.bw <- clv.DensBw(x, cluster_index, scatt)
    
    # compute SD and S_Dbw indicies
    clv.SD(scatt$Scatt, dis, alfa = numclust) # alfa is equal to number of clusters
  })


mySDbw <-
  sapply(minNST:maxNST, function(numclust){
    cluster_index <- cutree(clus_sites, numclust)
    SDbw(data = x, numclusters = numclust, cluster_vector = cluster_index)
  })


result <-
  NbClust::NbClust(
    data = x,
    diss = NULL,
    distance = 'euclidean',
    min.nc = minNST,
    max.nc = maxNST,
    method = 'ward.D2',
    index = 'sdbw'
  )

df = tibble(clv.sd = SD, nbclust.sd = result$All.index, NST = minNST:maxNST)

plot <-
  df %>%
  ggplot() +
  geom_line(aes(NST, clv.sd), color = red) +
  geom_point(aes(NST, clv.sd), color = red) +
  geom_line(aes(NST, nbclust.sd), color = blue) +
  geom_point(aes(NST, nbclust.sd), color = blue)

gridExtra::grid.arrange(plot)

plot(x)











