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

do.data = 1
do.plots = 0
fdp = 'lap'

if(do.clustering.analysis){
  
  if(do.data){
    mypc = (Sys.info()['sysname'] == 'Windows') || (Sys.info()['sysname'] =='Darwin')
    
    cores = if(mypc) 4 else detectCores() - 10
    plan(multisession, workers = cores)
    save_directory = paste0('~/SpatialNiche/', save_date, '/')
    
    source('~/SpatialNiche/Code/clustering_functions.R')
  
    Ly = 500
    if(fdp == 'bci'){
      Lx = 1000
      bci =
        readRDS(
          url(
            'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
          )
        )
      filename = paste0(save_directory, 'bci_clustering_analysis.rds')
      
    }
    if(fdp == 'lap') 
    {Lx = 500
    bci  = read_all_laplanada()
    filename = paste0(save_directory, 'lap_clustering_analysis.rds')
    }
    
    parameters = 
      expand_grid(
        thecensus = 1:2,
         algorithm = c('louvain', 'walktrap'),
         #algorithm = c('louvain'),
        d_cutoff = seq(10, 30, by = 2),
        self_loops = FALSE,
        d_step = 1e-5,
        Lx = Lx,
        Ly = 500,
        autolinked_species_only = TRUE,
        #weighted = c(TRUE),
        weighted = c(TRUE, FALSE),
        seed = 0:10
      ) %>%
      filter(algorithm=='walktrap'|seed==0)
    
    chunk_size = cores
    
    for(piece in seq(nrow(parameters) / chunk_size)){
      
      indices = chunk_size * (piece - 1) + seq(chunk_size)
      
      parms = parameters[indices, ]
      
      bci_analyzed = 
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
              bci #%>%
              #filter(census == thecensus)
            
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
        )
      
      save_date = gsub('-', '', Sys.Date())
      save_directory = paste0('~/SpatialNiche/', save_date, '/')
      dir.create(save_directory, showWarnings = FALSE)
      

      if(file.exists(filename)){
        bci_analyzed = 
          readRDS(filename) %>%
          bind_rows(bci_analyzed)
      }
      
      saveRDS(bci_analyzed, file = filename)
      
    }
  }
  
  if(do.plots){
    data = readRDS('~/SpatialNiche/20210820/lap_clustering_analysis.rds')
    
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

if(do.recruitment.analysis){
  
  library(sparr) # for function bivariate.density()
  
  ## Determine whether working on SeaWulf (SBU hpc) or personal computer
  seawulf = as.logical(Sys.info()['user'] == 'rdrocha')
  cores = if(seawulf) detectCores() else 6
  plan(multisession, workers = cores)
  
  if(!exists('combined_data')){
    census_data = 
      readRDS('c:/users/rdand/Google Drive/GitHub/Spatial-niche/Data/all_data.rds')
    
    cluster_data = 
      readRDS('~/SpatialNiche/Data/20210722/bci_clustering_analysis.rds') %>%
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
  
  if(!file.exists('~/SpatialNiche/Data/bci_inferred_soiltypes.rds')){
    KernelDensityEstimation = 
      function(gx, gy, Lx = 1000, Ly = 500, quadrat_length = 20, ...){
        
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
      function(Census, Algorithm, Seed, D_cutoff, Group){
        df = 
          combined_data %>% 
          filter(
            census == Census,
            algorithm == Algorithm,
            seed == Seed,
            d_cutoff == D_cutoff,
            group == Group
          ) %>%
          select(gx, gy)
        
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
        .options = furrr_options(seed = TRUE)
      )
    
    # kde_full = 
    #   combined_data %>%
    #   filter(
    #     algorithm == 'louvain',
    #     seed == 0
    #   ) %>%
    #   group_by(
    #     algorithm,
    #     seed,
    #     d_cutoff,
    #     group,
    #     census
    #   ) %>%
    #   future_pmap_dfr(
    #     .f = KernelDensityEstimation,
    #     .options = furrr_options(seed = TRUE)
    #   )
    
    
    soiltype = 
      kde_full %>%
      group_by(algorithm, seed, d_cutoff, census, x, y) %>%
      slice_max(density, n = 1) %>%
      ungroup %>%
      rename(soiltype = group) %>%
      select(-density)
    
    kde_full %<>%
      left_join(
        soiltype,
        by = c('x', 'y', 'census', 'algorithm', 'seed', 'd_cutoff')
      ) %>% 
      mutate(fdp = 'bci')
    
    saveRDS(kde_full, file = '~/SpatialNiche/Data/bci_inferred_soiltypes.rds')
    
  }else{
    kde_full = readRDS('~/SpatialNiche/Data/bci_inferred_soiltypes.rds')
  }
  
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
  
  res =
    rec_df %>%
    group_by(
      labelling_census, 
      algorithm, 
      seed, 
      d_cutoff, 
      fdp, 
      weighted, 
      autolinked_species_only,
      self_loops,
      d_step,
      number_of_groups,
      group
    ) %>%
    summarize(
      recruits = n(),
      matches = sum(group == soiltype),
      nsoiltypes = unique(number_of_groups),
      theta = (nsoiltypes - 1) / (recruits / matches - 1),
      .groups = 'drop'
    )
  
  res_summary = 
    res %>%
    group_by(
      d_cutoff, 
      algorithm, 
      seed, 
      fdp, 
      weighted, 
      autolinked_species_only,
      self_loops,
      d_step,
    ) %>%
    summarize(
      theta_mean = mean(theta),
      theta_se = sd(theta) / sqrt(n()),
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
    nutrients = 
      read_excel(
        '~/SpatialNiche/Data/bci.block20.data-original.xls', 
        sheet = 2
      ) %>%
      pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
      group_by(nutrient) %>%
      mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
      ungroup
  }
    nutrients_wide = 
      nutrients %>%
      select(-value) %>%
      pivot_wider(names_from = nutrient, values_from = standardized)
    
    kde_full = readRDS('~/SpatialNiche/Data/bci_inferred_soiltypes.rds')
    
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
    
    MLres = 
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
            select(Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, `N(min)`, pH, soiltype) %>%
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
    
    saveRDS(MLres, file = '~/SpatialNiche/Data/bci_C5_soiltype_vs_nutrients.rds')
    
  }  
  
  if(do.plots){
    
    MLres = readRDS('~/SpatialNiche/Data/bci_C5_soiltype_vs_nutrients.rds')
    
    MLres_summary = 
      MLres %>% 
      group_by(d_cutoff) %>% 
      summarize(
        mean = mean(Kappa), 
        se = sd(Kappa) / sqrt(n()), 
        .groups = 'drop'
      )
    
    plot_bars = 
      MLres_summary %>% 
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
    
    
    nutrients = 
      read_excel(
        '~/SpatialNiche/Data/bci.block20.data-original.xls', 
        sheet = 2
      ) %>%
      pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
      group_by(nutrient) %>%
      mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
      ungroup
    
    kde_full = readRDS('~/SpatialNiche/Data/bci_inferred_soiltypes.rds')
    
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
  
  
  trait_data_raw = read_excel('~/BCI/BCITRAITS_20101220.xlsx')
  
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
  cluster_data = 
    readRDS('~/SpatialNiche/Data/20210722/bci_clustering_analysis.rds') %>%
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
  
  nutrients = 
    read_excel(
      '~/SpatialNiche/Data/bci.block20.data-original.xls', 
      sheet = 2
    ) %>%
    pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
    group_by(nutrient) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup
  
  kde_full = readRDS('~/SpatialNiche/Data/bci_inferred_soiltypes.rds')
  
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
    inner_join(ranked_groups) %>%
    inner_join(pca_data)
  
  dcut = 10
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
      fdp,
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
    census_data = 
      readRDS('c:/users/rdand/Google Drive/GitHub/Spatial-niche/Data/all_data.rds')
    
    cluster_data = 
      readRDS('~/SpatialNiche/Data/20210722/bci_clustering_analysis.rds') %>%
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
