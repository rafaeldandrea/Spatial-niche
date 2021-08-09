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
# load('/Users/wangdz/Downloads/Spatial-niche-main-2/Data/colors.rdata')


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


distances = 0:ceiling(sqrt(2) * 144)/144
pdexp = sapply(distances, function(d) euclidean2d_distance_probability(d, 1))

cumulative_null_prob = cumsum(pdexp * c(0, diff(distances)))

#setwd('~/SpatialNiche/Data/20200709-80scen-A20kS220-H-betanoise-biotime2e5/')
#datafiles = paste0('scenario_', 1:80,'.RData')



if(do.fdp){
  if(TRUE){
    bci_all = NULL
    for(i in 1:7){
      bci_all = 
        bci_all %>%
        rbind(
          get(
            load(
              url(
                paste0("https://github.com/rafaeldandrea/BCI/blob/master/bci.full",i,".rdata?raw=true")
              )
            )
          ) %>%
            as_tibble %>%
            filter(!is.na(gx)) %>%
            filter(!is.na(gy)) %>%
            mutate(census = i)
        )
    }
    bci_recruit = 
      bci_all %>%
      filter(dbh > 100) %>% # keeping only the adult trees
      group_by(treeID)%>%
      slice_min(census) %>%
      filter(census!=1) %>%
      select("sp","gx","gy",'census') %>%
      ungroup()
  }

  #these are the final results
  verify_result=NULL
  hist_results = NULL
  ##using for loop to go through census 1-7
  for(b in 1:7){  
    
  distance_threshold = 10 ## d* in meters
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
  

  set.seed(0)
  
 
  bci = 
    get(
      load(
        url(
          paste0("https://github.com/rafaeldandrea/BCI/blob/master/bci.full",b,".rdata?raw=true")
        )
      )
    )  %>%
    dplyr::filter(dbh >= 100) %>%
    select(sp, gx, gy) %>%
    unique() %>%
    as_tibble
  
  if(fdp == 'bci') dat = bci 
  if(fdp == 'lap'){
    dat = lap
    lap = 
      read.table(
        '~/SpatialNiche/Data/relaplanadacensusdata/Censo 2_La Planada.txt', 
        header = TRUE
      ) %>%
      as_tibble %>%
      filter(dbh >= 100) %>%
      select(sp, gx, gy)
  }
  
  
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
      theme(aspect.ratio = ifelse(fdp == 'bci', .5, 1)) #+
      # scale_color_manual(values = with(colors,c(red, green, blue, yellow)))
    
    plot_fdp_null = 
      dat %>%
      mutate(sp = as.character(sp)) %>%
      inner_join(membership_tibble, by = 'sp') %>%
      ggplot(aes(gx, gy, color = nullcommunity, group = nullcommunity)) +
      geom_point() +
      labs(color = 'cluster') +
      theme(aspect.ratio = ifelse(fdp == 'bci', .5, 1)) #+
      # scale_color_manual(values = with(colors,c(red, green, blue, yellow)))
    
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
  
"1.label every species by the inferred sniche according to data from census 1. 
  Also label every quadrat by the inferred soil type based on the kde, also with labels from census 1.
 2. With the labels from Step 1, look at every recruitment event that took place in any of the censuses (2-7),
    and write down whether the tree recruited in its own soil type or a different soil type.
 3. take the average of same type recruiting of censuses 2-7, calculate std across censuses 2-7
 4. Repeat steps 1 and 2 and 3 for labels from census 2, 3, 4, 5, 6, and 7."
  
  #Below are the code for running recruitment verification
  same_type_recruiting   = 1   
  if (same_type_recruiting){
    
    ##choose new recruit trees 
    ##  has dbh greater than 100 in census T but less than 100 in census T-1. 
    results_by_majority = results %>% 
      group_by(x, y) %>% 
      slice_max(prob)  %>% 
      ungroup()
  
   ##parameters for gridding
    if(fdp == 'bci') {
      ncx = 50
      ncy  = 25
      gxmin = 0
      gxmax = 1000
      gymin = 0
      gymax = 500
    }
    
    #assign membership to young tress
    bci_recruit_member =  bci_recruit %>% inner_join(membership_tibble,by='sp')
    bci_single_census = dat_filtered %>% inner_join(membership_tibble,by='sp')
    ##grid the landscape into cells
    gxbreaks=seq(gxmin,gxmax,length=ncx+1)
    cellindex_x=findInterval(bci_recruit_member$gx,gxbreaks,rightmost.closed=TRUE)
    gybreaks=seq(gymin,gymax,length=ncy+1)
    cellindex_y=findInterval(bci_recruit_member$gy,gybreaks,rightmost.closed=TRUE)
    
    #assign cell numbers to young trees
    bci_recruit_member=
      bci_recruit_member %>%
      mutate(
        cx=cellindex_x,
        cy=cellindex_y,
        cellindex=cx+(cy-1)*ncx
      )
    #assign cell numbers to inferred soil type
    results_by_majority =
      results_by_majority %>%
      mutate(
        cx= (x-10)/20+1,
        cy = (y-10)/20+1,
        cellindex=cx+(cy-1)*ncx
      )
    ##assign each tree with soil type, according to cell number
    bci_recruit_soiltype = bci_recruit_member %>% 
      inner_join(results_by_majority , by='cellindex') %>%
      select('sp','gx','gy','census','community','nullcommunity','soiltype','cellindex') %>%
      mutate(recruit = community == soiltype ) 
    
    ##bci_recruit_soiltype contains all new recruitsin every census and 
    ##whether they recruit in the same soiltype
    
    ##calculate the percentage of same type recruitment happens on  census 2-7 for each soiltype
    ##percentage_recruits contains the chance for same type recruitment 
    percentage_recruits = bci_recruit_soiltype %>% 
      group_by(census,community) %>% 
      summarise(correct_recruit = sum(recruit==TRUE),
                wrong_recruit = sum(recruit==FALSE),
                percentage = correct_recruit/(correct_recruit+wrong_recruit))%>%
                ungroup()
    ##calculate mean and std, across census 2-7, with census 1-7 labeling
    percent_tibble =  percentage_recruits%>%
      group_by(community) %>% 
      summarise(mean=mean(percentage),std=sd(percentage))%>%
      ungroup()%>%
      mutate(census=b)
    
    percentage_recruits$community = as.numeric(percentage_recruits$community)
    ## preparing data for plotting histogram,  exp = 1/cluster number
    hist  =  percentage_recruits %>% group_by(census) %>% 
      mutate(expect = 1/max(community)) %>% 
      ungroup()%>% 
      mutate(delta = (percentage - expect) /expect )
    
    ##bind results using the for loop to go through census 1-7
    hist_results = hist_results %>%
      rbind(hist)
    
    verify_result = 
      verify_result %>%
      rbind(percent_tibble)
    
  
    
    
  }
  
  }}

#bar plots with error bars
verify_result%>%
  ggplot(aes(census, mean,fill = community)) +
  geom_bar(position = "dodge",
           stat = "identity")+
  geom_errorbar(aes(x=census, ymin=mean-std, ymax=mean+std), width=.1,
                position=position_dodge(.9)) +
  scale_y_continuous(name="Pencentage" ,breaks = c(0.1,0.2,0.3,0.4,0.5,0.6), limits=c(0,0.6)) +
  scale_x_discrete(name="Census" ,limits=c("1","2","3","4","5","6","7"))+
  labs(title=paste0("Same type recruiting, distance thresh = ", distance_threshold, "m"))+
  theme_classic() 
#hist plots
hist_results %>%
ggplot(aes(delta))+
  geom_histogram(binwidth=0.05)+
  labs(title=paste0("Same type recruiting, distance thresh = ", distance_threshold, "m"))+
  #scale_y_continuous(name="count",breaks = c(2,4,6,8,10,12,14),limits=c(0,15))+
  #xlab("(obs-exp)/exp")+
  theme_classic() 
closeAllConnections()
