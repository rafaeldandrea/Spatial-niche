library(tidyverse)
library(magrittr)
library(furrr)
library(parallel)
library(readxl)
library(parallelDist) ## for function parDist()
library(igraph) ## to convert the adjacency matrix constructed from the pairs connected into 
## a graph then find communities in the graph
library(tidyverse)
library(furrr)  ## for parallel computing
library(parallel)  
library(magrittr)
library(RANN)  ## for neighbor-finding function nn2()
library(sparr)

yas_link="https://github.com/rafaeldandrea/Spatial-niche/tree/main/Data/YasuniData"


#Read soil data

nutrients=read.csv(
  url("https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/resoilnutrientdatarequest/yas_20x20_soil.csv?raw=true"
  )
)%>%as_tibble()


#Read tree data

yas=read.csv(
  url('https://github.com/rafaeldandrea/Spatial-niche/tree/main/Data/YasuniData/yasuni.csv?raw=true'
      ))%>%
  as_tibble()

Lx=500
Ly=500

if(do.clustering.analysis){
  
      #Baldeck cutoffs
      baldeck_cutoff = 
        yas %>%
        group_by(sp) %>%
        summarize(
          baldeck = quantile(dbh, .56), 
          .groups ='drop'
        )

  dat = 
  yas %>%
  inner_join(baldeck_cutoff) %>%
  filter(dbh > baldeck) %>%
  mutate(fdp = 'yas ')
    
    parameters = 
     expand_grid(
    algorithm = 'louvain',
    d_cutoff = c(20),
    Lx = 500,
    Ly = 500,
    weighted = TRUE,
    seed = 0
    ) 

    


    
    fdp_analyzed = 
      parameters %>%
      future_pmap_dfr(
        .f = function(
          algorithm,
          d_cutoff,
          Lx,
          Ly,
          weighted,
          seed
        ){
          
          if(seed > 0){
            data %<>%
              mutate(sp = sample(sp))
          }
          
          result = 
            dat %>%
            adjacency_matrix(
              d_cutoff = d_cutoff, 
              Lx = Lx, 
              Ly = Ly,
              method='ann'
            ) %>%
            cluster_analysis(
              algorithm = algorithm,
              weighted = weighted
            ) %>%
            pluck('result') %>%
            mutate(
              d_cutoff = d_cutoff,
              seed = seed
            )%>%
            return()
        },
        .options = furrr_options(seed = TRUE)
      ) %>%
      rename(sp = name)
    
    dir.create(save_date, showWarnings = FALSE)
    
    
    saveRDS(fdp_analyzed, file = filename)
    
    cluster_data = readRDS(filename)
    
 }
    

if(do.kde.analysis){
  cluster_data=readRDS("20220114yas_clustering_analysis.rds")
  
  
  #Baldeck cutoffs
  baldeck_cutoff = 
    yas %>%drop_na()%>%
    group_by(sp) %>%
    summarize(
      baldeck = quantile(dbh, .56), 
      .groups ='drop'
    )
  
  
  dat <-yas%>%
    inner_join(baldeck_cutoff) %>%
    filter(dbh > baldeck) %>%
    mutate(fdp = 'yas')%>%select(-baldeck)
  
  
  census_data = dat
  

  data = 
    census_data %>% 
    select(sp, gx, gy) %>%
    full_join(
      cluster_data %>%
        filter(seed == 0) %>%
        select(sp, d_cutoff, group)
    ) %>%drop_na()%>%
    select( d_cutoff, gx, gy, group)%>%mutate(census=7)
  
  
  
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
  
  kde_full<-kde_full%>%
    inner_join(
      kde_full %>%
        group_by(d_cutoff, x, y) %>%
        slice_max(density, n = 1, with_ties = FALSE) %>%
        rename(soiltype = group) %>%
        ungroup() %>%
        select(-density)%>%
        mutate(fdp='yas')
    )
  
  kde_full$soiltype<-factor(kde_full$soiltype,levels=c(1,2,3,4))
  saveRDS(kde_full,"20220115yas_kde_full.RDS")
  
}



if(do.kde.AlP.analysis){
  
  ## The story here is that our plant groups segregate in the 3-d space 
  ## whose axes are correlated nutrients + Al + P. Group 3 is associated
  ## with high nutrient concentration, while the other 3 groups are associated
  ## with low nutrient concentration. Groups 2 and 3 distinguish themselves
  ## on the P vs Al plane. Unclear how Group 1 differs from 2 and 3.
  
  plants =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/20211205/20211205_kde_full.rds?raw=true'
      )
    ) %>%
    filter(
      census == 8,
      d_cutoff == 20
    ) %>%
    select(
      x, y, group, density
    ) %>%
    pivot_wider(names_from = group, values_from = density)
  
  names(plants) = c('x', 'y', 'g1', 'g2', 'g3', 'g4')
  
  water =
    read.table(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
      ),
      header = TRUE
    ) %>%
    as_tibble() %>%
    mutate(
      x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
    ) %>%
    replace_na(list(x = 10, y = 10)) %>%
    group_by(x, y) %>%
    summarize(water = mean(swp), .groups = 'drop')
  
  nutrients =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/Manuscript-Data/bci_nutrient_data.rds?raw=true'
      )
    ) %>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )
  
  data = 
    plants %>%
    inner_join(nutrients) 
  
  data_simplified =
    data %>%
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density') %>%
    group_by(x, y) %>%
    slice_max(density) %>%
    ungroup() %>%
    select(-density)
  
  data_simplified %>%
    ggplot(aes(Al, P)) +
    geom_point() +
    facet_wrap(~group)
  
  BD = 
    function(plantgroup, h0 = .3, data = data_simplified){
      if(plantgroup %in% unique(data_simplified$group)){
        foo = with(
          data %>% 
            filter(group == plantgroup), 
          ppp(Al, P, xrange = range(Al), yrange = range(P))
        )  
      } else{
        foo = with(
          data,
          ppp(Al, P, xrange = range(Al), yrange = range(P))
        )
      }
      
      bar = as.numeric(bivariate.density(foo, h0 = h0)$z$v)
      
      return(bar)
    }
  
  Al_min = min(data$Al)
  Al_max = max(data$Al)
  P_min = min(data$P)
  P_max = max(data$P)
  
  res = 
    expand_grid(
      Al = seq(Al_min, Al_max, length = 128),
      P = seq(P_min, P_max, length = 128)
    ) %>%
    arrange(Al, P) %>%
    bind_cols(
      tibble(
        g0 = BD('all'),
        g1 = BD('g1'),
        g2 = BD('g2'),
        g3 = BD('g3'),
        g4 = BD('g4')
      )
    )
  
  plot_densities = 
    res %>%
    filter(g0 > 0.01) %>%
    select(-g0) %>%
    pivot_longer(g1:g4, names_to = 'group') %>%
    group_by(group) %>%
    mutate(value = scale(value), .groups = 'drop') %>%
    ggplot() +
    geom_raster(aes(Al, P, fill = value)) +
    geom_point(aes(Al, P), data = data_simplified) +
    facet_wrap(~group) +
    scale_fill_gradientn(colors = terrain.colors(2000)) +
    theme(aspect.ratio = 1) +
    labs(fill = 'scaled\ndensity')
  
  set.seed(0)
  pca = 
    nutrients %>%
    filter(x<500)%>%
    select(Al:pH) %>% 
    pcaMethods::pca(nPcs = 2, scale = 'uv', center = TRUE)
  
  
  df_scores = 
    pca@scores %>% 
    as_tibble() %>% 
    bind_cols(data_simplified) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  df_loadings = 
    pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(pca@loadings))) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  plot = 
    ggplot() + 
    geom_point(aes(PC1, PC2, color = group), data = df_scores) + 
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings)), 
        y = rep(0, nrow(df_loadings)), 
        xend = 10*PC1, 
        yend = 10*PC2
      ), 
      data = df_loadings, 
      arrow = arrow(length = unit(.5, "cm"))
    ) +
    geom_text(
      aes(10 * PC1, 10 * PC2, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC1,
      nudge_y = .5 * df_loadings$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85))
  
  BD_PC = 
    function(plantgroup, h0 = .3, data = df_scores){
      if(plantgroup %in% unique(data_simplified$group)){
        foo = with(
          data %>% 
            filter(group == plantgroup), 
          ppp(PC1, PC2, xrange = range(PC1), yrange = range(PC2))
        )  
      } else{
        foo = with(
          data,
          ppp(PC1, PC2, xrange = range(PC1), yrange = range(PC2))
        )
      }
      
      bar = 
        as.numeric(
          bivariate.density(
            foo, 
            h0 = h0, 
            xy = 
              list(
                x = seq(min(df_scores$PC1), max(df_scores$PC1), length = 128),
                y = seq(min(df_scores$PC2), max(df_scores$PC2), length = 128)
              ))$z$v
        )
      
      return(bar)
    }
  
  res = 
    expand_grid(
      PC1 = seq(min(df_scores$PC1), max(df_scores$PC1), length = 128),
      PC2 = seq(min(df_scores$PC2), max(df_scores$PC2), length = 128)
    ) %>%
    arrange(PC1, PC2) %>%
    bind_cols(
      tibble(
        g0 = BD_PC('all'),
        g1 = BD_PC('g1'),
        g2 = BD_PC('g2'),
        g3 = BD_PC('g3'),
        g4 = BD_PC('g4')
      )
    ) %>%
    replace_na(list(g1 = 0, g2 = 0, g3 = 0, g4 = 0)) %>%
    mutate(across(g1:g4, function(x) (x - min(x)) / (max(x) - min(x))))
  
  
  plot_densities = 
    res %>%
    filter(g0 > 1e-6) %>%
    select(-g0) %>%
    pivot_longer(g1:g4, names_to = 'group') %>%
    ggplot() +
    geom_tile(aes(PC1, PC2, fill = value)) +
    geom_point(aes(PC1, PC2), data = df_scores, alpha = .2) +
    facet_wrap(~group) +
    scale_fill_gradientn(colors = terrain.colors(2000)) +
    theme(aspect.ratio = 1) +
    labs(fill = 'scaled\ndensity') +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings)), 
        y = rep(0, nrow(df_loadings)), 
        xend = 10*PC1, 
        yend = 10*PC2
      ), 
      data = df_loadings, 
      arrow = arrow(length = unit(.5, "cm"))
    ) +
    geom_text(
      aes(10 * PC1, 10 * PC2, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC1,
      nudge_y = .5 * df_loadings$PC2,
      color = 'red'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) 
}

if(do.pca.analysis){
  

  plants=readRDS(
    url("https://github.com/rafaeldandrea/Spatial-niche/tree/main/Data/YasuniData/20220115yas_kde_full.RDS?raw=true")
    )%>%
    filter(d_cutoff==20)%>%
    select(x,y,group,density)%>%
    pivot_wider(names_from = group, values_from = density)
  
  names(plants) = c('x', 'y', 'g1', 'g2', 'g3', 'g4')
  
  water =
    read.table(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/BCI_SWP_map_mid_dry_season_regular.txt'
      ),
      header = TRUE
    ) %>%
    as_tibble() %>%
    mutate(
      x = seq(10, 990, 20)[cut(x, breaks = seq(0, 1000, 20), labels = FALSE)],
      y = seq(10, 490, 20)[cut(y, breaks = seq(0, 500, 20), labels = FALSE)]
    ) %>%
    replace_na(list(x = 10, y = 10)) %>%
    group_by(x, y) %>%
    summarize(water = mean(swp), .groups = 'drop')
  
  
  nutrients=read.csv(
    url("https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/resoilnutrientdatarequest/yas_20x20_soil.csv?raw=true"
    )
  )%>%as_tibble()%>%
    mutate(
      across(
        Al:pH, 
        function(x) scale(x)[, 1]
      )
    )
  
nutrients%<>%
  mutate(across(
    c(x,y),
  function(x) 20*x-10
  ))

  
  data = 
    plants %>%
    inner_join(nutrients) 
  
  data_simplified =
    data %>%
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density') %>%
    group_by(x, y) %>%
    slice_max(density) %>%
    ungroup() %>%
    select(-density)
  
  plot_BCI = 
    df_scores %>% 
    ggplot(aes(x, y, fill = group))+
    geom_tile()+
    theme(aspect.ratio = .5) +
    scale_fill_manual(
      values = 
        c(
          colors$red, 
          colors$green, 
          colors$blue, 
          colors$yellow)
    )
  
  set.seed(0)
  pca = 
    nutrients %>% 
    inner_join(water) %>%
    select(Al:water) %>% 
    pcaMethods::pca(nPcs = 4, scale = 'uv', center = TRUE)
  
  df_scores = 
    pca@scores %>% 
    as_tibble() %>% 
    bind_cols(data_simplified) %>%
    mutate(PC1 = -PC1, PC2 = -PC2) 
  
  if(do.kmeans.analysis){
    k = 
      future_pmap_dfr(
        expand_grid(
          groups = 2:10, 
          seed = 0:100
        ),
        function(groups, seed){
          foo = 
            df_scores %>% 
            select(Al:water)
          if(seed > 0){
            set.seed(seed)
            foo = apply(foo, 2, sample)
          }
          wss = 
            kmeans(
              foo, 
              centers = groups, 
              nstart = 100
            )$tot.withinss
          return(tibble(groups = groups, seed = seed, wss = wss))
        } 
      )
    
    k4 = 
      kmeans(
        df_scores %>% 
          select(Al:water),
        centers = 4,
        nstart = 100
      )
    
    df_scores = 
      df_scores %>%
      mutate(kmeans = k4$cluster)
    
    plot_mapping = 
      df_scores %>% 
      count(group, kmeans) %>%
      ggplot(aes(kmeans, n)) +
      geom_col() +
      facet_wrap(~group) +
      labs(
        x = 'soil group',
        y = 'quadrat count'
      ) +
      theme(aspect.ratio = 1)
    
    plot_gap = 
      k %>% 
      filter(seed == 0) %>% 
      inner_join(
        k %>% 
          filter(seed > 0) %>% 
          group_by(groups) %>% 
          summarize(
            null = mean(wss), 
            .groups = 'drop')
      ) %>% 
      mutate(gap = null - wss) %>%
      ggplot(aes(groups, gap)) +
      geom_line() + 
      geom_point() +
      xlab('number of groups') +
      theme(aspect.ratio = 1)
    
    grid.arrange(plot_gap, plot_mapping, ncol = 2)
    
  }
  
  df_loadings = 
    pca@loadings %>% 
    as_tibble() %>% 
    bind_cols(tibble(feature = rownames(pca@loadings))) %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  
  plot_PCA = 
    ggplot() + 
    geom_point(aes(PC1, PC2, color = group), data = df_scores) + 
    theme(aspect.ratio = 1) +
    geom_segment(
      aes(
        x = rep(0, nrow(df_loadings)), 
        y = rep(0, nrow(df_loadings)), 
        xend = 10*PC1, 
        yend = 10*PC2
      ), 
      data = df_loadings, 
      arrow = arrow(length = unit(.25, "cm"))
    ) +
    geom_text(
      aes(10 * PC1, 10 * PC2, label = feature), 
      data = df_loadings, 
      nudge_x = .5 * df_loadings$PC1,
      nudge_y = .5 * df_loadings$PC2,
      color = 'darkgreen'
    ) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    labs(
      x = 'PC1',
      y = 'PC2'
    ) +
    theme(legend.position = c(.85, .85)) +
    scale_color_manual(
      values = 
        c(
          colors$red, 
          colors$green, 
          colors$blue, 
          colors$yellow)
    )
  
  plot_group = function(gr){
    plot = 
      df_scores %>%
      ggplot() + 
      geom_point(
        aes(PC1, PC2), 
        data = df_scores %>% filter(group != gr),
        color = 'grey',
        alpha = .3
      ) + 
      geom_point(
        aes(PC1, PC2), 
        data = df_scores %>% filter(group == gr),
        color = 'blue'
      ) +
      theme(aspect.ratio = 1) +
      geom_hline(yintercept = 0, color = 'gray') +
      geom_vline(xintercept = 0, color = 'gray') +
      labs(
        x = 'PC1',
        y = 'PC2'
      ) +
      coord_cartesian(
        xlim = range(c(10 * df_loadings$PC1, df_scores$PC1)), 
        ylim = range(c(10 * df_loadings$PC2, df_scores$PC2))
      )
    
    return(plot)
  }
  
  p1 = plot_group('g1')
  p2 = plot_group('g2')
  p3 = plot_group('g3')
  p4 = plot_group('g4')
  
  plot_grid(plot, plot_grid(p1, p2, p3, p4), axis = 'tb')  
  
  s2d = 
    plot(
      df_scores[,1:2], 
      pch = 20, 
      col = (1:4)[match(df_scores$group, c('g1', 'g2', 'g3', 'g4'))]
    )
  s3d = 
    scatterplot3d(
      df_scores[,1:3], 
      pch = 20, 
      color = (1:4)[match(df_scores$group, c('g1', 'g2', 'g3', 'g4'))]
    )
  
  plot_loadings_angle = 
    df_loadings %>%
    mutate(angle = 360 / (2 * pi) * atan(PC2/ PC1)) %>%
    mutate(feature = factor(feature, levels = feature[order(angle)])) %>%
    ggplot(aes(feature, angle, label = feature)) +
    geom_text() +
    theme(legend.position = 'none')
  
  scatterplot3d(
    df_scores$`N(min)`, 
    df_scores$P, 
    df_scores$Fe, 
    pch = 20,
    c(2, 3, 4, 6)[match(df_scores$group, c('g1', 'x', 'x','g4'))]
  )
  
  data %>% 
    pivot_longer(g1:g4, names_to = 'group', values_to = 'density') %>% 
    pivot_longer(Al:water, names_to = 'feature') %>% 
    group_by(feature, group) %>% 
    summarize(
      cor = cor(density, value), 
      .groups = 'drop'
    ) %>% 
    ggplot(aes(feature, cor)) + 
    facet_wrap(~group) + 
    geom_col()
  
  # feature groups
  # group 1: B K Ca Zn Mg N(min) Cu
  # group 2: Mn Fe pH N
  # group 3: Al water P
  
  p124 = 
    df_scores %>% 
    filter(group != 'g3') %>% 
    ggplot(aes(pH, P, color = group)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    scale_color_manual(values = c(colors$red, colors$green, colors$yellow)) +
    theme(aspect.ratio = 1)
  
  pAl = 
    df_scores %>% 
    ggplot(aes(pH, Al, color = group)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    scale_color_manual(values = c(colors$red, colors$green, colors$blue, colors$yellow)) +
    theme(aspect.ratio = 1)
  
  p14 = 
    df_scores %>% 
    filter(group %in% c('g1', 'g4')) %>% 
    ggplot(aes(N, Fe, color = group)) + 
    geom_point() +
    ggtitle('G1 vs G4') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray') +
    scale_color_manual(values = c(colors$red, colors$yellow))
  
  p13 = 
    df_scores %>% 
    filter(group %in% c('g1', 'g3')) %>% 
    ggplot(aes(PC1, Al, color = group)) + 
    geom_point() +
    ggtitle('G1 vs G3') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')
  
  p12 = 
    df_scores %>% 
    filter(group %in% c('g1', 'g2')) %>% 
    ggplot(aes(N, P, color = group)) + 
    geom_point() +
    ggtitle('G1 vs G2') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')
  
  p24 = 
    df_scores %>% 
    filter(group %in% c('g2', 'g4')) %>% 
    ggplot(aes(N, P, color = group)) + 
    geom_point() +
    ggtitle('G2 vs G4') +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_vline(xintercept = 0, color = 'gray')
  
  plot_grid(p12, p13, p14, p24)
  
}


