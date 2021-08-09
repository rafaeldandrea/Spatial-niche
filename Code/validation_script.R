library(tidyverse)
library(magrittr)
library(RandomFields)
library(furrr)
library(parallel)

## Determine whether working on SeaWulf (SBU hpc) or personal computer
seawulf = as.logical(Sys.info()['user'] == 'rdrocha')

do.data = 1
do.plots = 0
do.plots2 = 0
do.envelope.test = 0

if(do.data){
  cores = if(seawulf) detectCores() else 7
  plan(multisession, workers = cores)
  
  source('~/SpatialNiche/R_Scripts/clustering_functions.R')
  
  generate_landscape =
    function(
      Lx,
      Ly,
      quadrat_length,
      rangepar,
      sillpar,
      nuggetpar,
      num_soiltypes,
      seed
    ) {
      stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
      RFoptions(seed = seed)
      
      if(num_soiltypes == 1){
        soiltype = 1
      }else{
        stress = 
          RFsimulate(
            RMgauss(
              scale = rangepar + 1e-16,
              var = 2 * sillpar * (1 - nuggetpar)
            ) + 
              RMtrend(mean = 0) + 
              RMnugget(var = 2 * sillpar * nuggetpar),
            x = seq(Lx / quadrat_length),
            y = seq(Ly / quadrat_length)
          )@data$variable1
        
        soiltype = 
          stress %>%
          cut(
            num_soiltypes, 
            labels = FALSE
          )
      }
      
      landscape = 
        expand_grid(
          y = seq(Ly / quadrat_length),
          x = seq(Lx / quadrat_length)
        ) %>%
        mutate(
          soiltype = as.factor(soiltype)
        )
      
      return(landscape)
    }
  
  recruitment = 
    function(
      soiltype, 
      species_group,
      species_abundance,
      theta
    ){
      prob_mat = 
        outer(
          soiltype, 
          species_group, 
          function(s, g) ifelse(s == g, theta, 1)
        )
      
      species = 
        apply(
          prob_mat, 
          1, 
          function(probs){
            sample(
              seq_along(probs), 
              size = 1, 
              prob = probs * species_abundance
            )
          } 
        )
      
      return(species)
    }
  
  build_recruits = 
    function(
      coords, 
      landscape, 
      species_group,
      species_abundance,
      Lx,
      Ly,
      quadrat_length,
      theta
    ){
      df =
        coords %>%
        mutate(
          x = 
            cut(
              gx, 
              breaks = seq(0, Lx, quadrat_length), 
              labels = FALSE
            ),
          y = 
            cut(
              gy, 
              breaks = seq(0, Ly, quadrat_length), 
              labels = FALSE
            )
        ) %>%
        left_join(
          landscape, 
          by = c('x', 'y')
        ) %>%
        mutate(
          species = 
            soiltype %>%
            recruitment(
              species_group = species_group,
              species_abundance = species_abundance,
              theta = theta
            )
        ) %>%
        return()
    }
  
  ## average intercensus mortality and births are each 
  ## typically ~= 10% of the number of living trees
  bci_all = readRDS('~/BCI/all_data.rds')
  
  censuses = 
    expand_grid(
      census1 = 1:7,
      census2 = 1:7
    ) %>%
    filter(census2 == census1 + 1)
  
  ndeaths = 
    censuses %>%
    pmap_dbl(
      .f = function(census1, census2){
        a = 
          bci_all %>% 
          filter(
            census == census1, 
            status == 'A',
            dbh >= 100
          )
        
        d = 
          bci_all %>% 
          filter(
            census == census2, 
            status == 'D'
          )
        
        length(intersect(a$treeID, d$treeID)) %>%
          return
      }
    )
  
  nbirths = 
    censuses %>%
    pmap_dbl(
      .f = function(census1, census2){
        a1 = 
          bci_all %>% 
          filter(
            census == census1, 
            status == 'A',
            dbh < 100
          )
        
        a2 = 
          bci_all %>% 
          filter(
            census == census2, 
            status == 'A',
            dbh >= 100
          )
        
        length(intersect(a1$treeID, a2$treeID)) %>%
          return
      }
    )
  
  mean_births = mean(c(nbirths, ndeaths))
  sd_births = sd(nbirths)
  
  mean_deaths = mean(c(nbirths, ndeaths))
  sd_deaths = sd(ndeaths)
  
  
  bci =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_all_censuses.rds?raw=true'
      )
    )
  
  foo = 
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/all_bci_censuses.rds?raw=true'
      )
    )
  
  bci %<>% 
    filter(sp %in% intersect(sp, foo$name))
  
  abuns = 
    bci %>%
    count(census, sp)
  
  number_of_species = 
    abuns %>% 
    filter(n >= 40) %>%
    pull(sp) %>% 
    unique() %>% 
    length()
  
  average_community_size = 
    bci %>%
    count(census) %>%
    pull(n) %>%
    mean()
  
  parameters =
    expand_grid(
      average_community_size = average_community_size,
      nspecies = number_of_species,
      nsoiltypes = c(2, 3, 4, 5, 6, 10, 15),
      ncensuses = 100,
      d_cutoff = c(10, 20), 
      d_step = 1e-5,
      Lx = 1000,
      Ly = 500,
      quadrat_length = 20,
      rangepar = 20, 
      sillpar = 1, 
      nuggetpar = .001, 
      seed = 0, 
      theta = c(1, 2, 5, 10, 1e5),
      clustering_algorithm = c('louvain', 'walktrap'),
      autolinked_species_only = TRUE,
      weighted = c(FALSE, TRUE),
      self_loops = FALSE
    )
  
  if(!seawulf){
    parameters %<>%
      filter(
        nsoiltypes %in% c(1, 10, 15)
      )
  }
  
  if(seawulf){
    parameters %<>% 
      filter(
        clustering_algorithm == 'louvain', 
        d_cutoff == 20, 
        rangepar == 20, 
        quadrat_length == 20, 
        weighted == TRUE,
        self_loops == FALSE
      )
  }
  
  simulation = 
    function(
      average_community_size,
      nspecies,
      nsoiltypes,
      ncensuses,
      d_cutoff, 
      d_step,
      Lx,
      Ly,
      quadrat_length,
      rangepar, 
      sillpar, 
      nuggetpar, 
      seed, 
      theta,
      clustering_algorithm,
      autolinked_species_only,
      weighted,
      self_loops
    ){
      census = 0
      
      species_list = seq(nspecies)
      
      species_group = 
        rep(seq(nsoiltypes), each = nspecies / nsoiltypes)
      
      delta = nspecies - length(species_group)
      if(delta > 0){
        species_group = c(rep(1, delta), species_group)
      }
      
      landscape = 
        generate_landscape(
          Lx = Lx, 
          Ly = Ly, 
          quadrat_length = quadrat_length, 
          rangepar = rangepar, 
          sillpar = sillpar, 
          nuggetpar = nuggetpar, 
          seed = seed, 
          num_soiltypes = nsoiltypes
        )
      
      community = 
        tibble(
          gx = runif(average_community_size, min = 0, max = Lx),
          gy = runif(average_community_size, min = 0, max = Ly)
        ) %>%
        build_recruits(
          landscape = landscape,
          species_group = species_group,
          species_abundance = 1,
          Lx = Lx,
          Ly = Ly,
          quadrat_length = quadrat_length,
          theta = theta
        )
      
      data = 
        tibble(
          census = census
        ) %>%
        bind_cols(community)
      
      census_deaths = 
        rnorm(
          ncensuses, 
          mean = mean_deaths,
          sd = sd_deaths
        )
      
      census_births = 
        rnorm(
          ncensuses, 
          mean = mean_births,
          sd = sd_births
        )
      
      while(census < ncensuses){
        census = census + 1
        ndeaths = census_deaths[census]
        nbirths = census_births[census]
        
        abundances = 
          community %>%
          count(species) %>%
          right_join(
            tibble(
              species = seq(nspecies)
            ),
            by = 'species'
          ) %>%
          replace_na(list(n = 0))
        
        births = 
          tibble(
            gx = runif(nbirths, min = 0, max = Lx),
            gy = runif(nbirths, min = 0, max = Ly)
          ) %>%
          build_recruits(
            landscape = landscape,
            species_group = species_group,
            species_abundance = abundances$n,
            Lx = Lx,
            Ly = Ly,
            quadrat_length = quadrat_length,
            theta = theta
          )
        
        community %<>%
          sample_n(size = nrow(community) - ndeaths) %>%
          bind_rows(births)
        
        if(census %% 10 == 0){
          data %<>%
            bind_rows(
              tibble(
                census = census
              ) %>%
                bind_cols(community)  
            )   
        }
        
      }
      
      return(data)
    }
  
  analyze = 
    function(
      data,
      average_community_size,
      nspecies,
      nsoiltypes,
      ncensuses,
      d_cutoff, 
      d_step,
      Lx,
      Ly,
      quadrat_length,
      rangepar, 
      sillpar, 
      nuggetpar, 
      seed, 
      theta,
      clustering_algorithm,
      autolinked_species_only,
      weighted,
      self_loops
    ){
      result = 
        data %>%
        pull(census) %>%
        unique() %>%
        future_map_dfr(
          .options = furrr_options(seed = TRUE),
          .f = function(cen){
            dat = 
              data %>%
              filter(census == cen) %>%
              rename(sp = species)
            
            A = 
              adjacency_matrix(
                dat, 
                autolinked_species_only = autolinked_species_only, 
                d_cutoff = d_cutoff, 
                d_step = d_step, 
                Lx = Lx, 
                Ly = Ly
              )
            
            clusters = 
              cluster_analysis(
                A = A, 
                algorithm = clustering_algorithm, 
                weighted = weighted, 
                self_loops = self_loops
              )
            
            clusters$result %>%
              mutate(census = cen) %>%
              return
          }
        ) %>%
        bind_cols(
          tibble(
            average_community_size,
            nspecies,
            nsoiltypes,
            ncensuses,
            d_cutoff, 
            d_step,
            Lx,
            Ly,
            quadrat_length,
            rangepar, 
            sillpar, 
            nuggetpar, 
            seed, 
            theta,
            autolinked_species_only
          )
        )
      
      return(result)
      
    }
  
  wrapper = 
    function(
      average_community_size,
      nspecies,
      nsoiltypes,
      ncensuses,
      d_cutoff, 
      d_step,
      Lx,
      Ly,
      quadrat_length,
      rangepar, 
      sillpar, 
      nuggetpar, 
      seed, 
      theta,
      clustering_algorithm,
      autolinked_species_only,
      weighted,
      self_loops
    ){
      data = 
        simulation(
          average_community_size,
          nspecies,
          nsoiltypes,
          ncensuses,
          d_cutoff, 
          d_step,
          Lx,
          Ly,
          quadrat_length,
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed, 
          theta,
          clustering_algorithm,
          autolinked_species_only,
          weighted,
          self_loops
        )
      
      result = 
        analyze(
          data,
          average_community_size,
          nspecies,
          nsoiltypes,
          ncensuses,
          d_cutoff, 
          d_step,
          Lx,
          Ly,
          quadrat_length,
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed, 
          theta,
          clustering_algorithm,
          autolinked_species_only,
          weighted,
          self_loops
        )
      
      return(result)
    }
  
  wrapper_simulation = 
    function(
      average_community_size,
      nspecies,
      nsoiltypes,
      ncensuses,
      d_cutoff, 
      d_step,
      Lx,
      Ly,
      quadrat_length,
      rangepar, 
      sillpar, 
      nuggetpar, 
      seed, 
      theta,
      clustering_algorithm,
      autolinked_species_only,
      weighted,
      self_loops
    ){
      data = 
        simulation(
          average_community_size,
          nspecies,
          nsoiltypes,
          ncensuses,
          d_cutoff, 
          d_step,
          Lx,
          Ly,
          quadrat_length,
          rangepar, 
          sillpar, 
          nuggetpar, 
          seed, 
          theta,
          clustering_algorithm,
          autolinked_species_only,
          weighted,
          self_loops
        )
      
      return(
        data %>% 
          bind_cols(
            tibble(
              average_community_size,
              nspecies,
              nsoiltypes,
              ncensuses,
              Lx,
              Ly,
              quadrat_length,
              rangepar, 
              sillpar, 
              nuggetpar, 
              seed, 
              theta
            )
          )
      )
    }
  
  if(seawulf){
    
    chunk_size = nrow(parameters)
    
    for(piece in seq(nrow(parameters) / chunk_size)){
      
      indices = chunk_size * (piece - 1) + seq(chunk_size)
      
      parms = parameters[indices, ]
      
      results = 
        parms %>%
        future_pmap_dfr(
          .f = wrapper_simulation,
          .options = furrr_options(seed = TRUE)
        )    
      
      save_date = gsub('-', '', Sys.Date())
      save_directory = paste0('~/SpatialNiche/Data/', save_date, '/')
      dir.create(save_directory, showWarnings = FALSE)
      
      filename = paste0(save_directory, 'simulation_data.rds')
      
      if(file.exists(filename)){
        results = 
          readRDS(filename) %>%
          bind_rows(results)
      }
      
      saveRDS(results, file = filename)
      
    }
    
  }
  

}

if(do.plots){
  
  theme_set(theme_bw())
  theme_update(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  result = 
    readRDS('~/SpatialNiche/Data/20210720/simulation_analysis.rds')
  
  plot_ngroups = 
    result %>%
    mutate(number_of_groups = factor(number_of_groups)) %>%
    group_by(nsoiltypes, number_of_groups) %>%
    count() %>%
    ungroup() %>%
    group_by(nsoiltypes) %>%
    mutate(p = n / sum(n) * 100) %>%
    ggplot(aes(number_of_groups, p, fill = number_of_groups)) +
    geom_col(color = 'black') +
    theme(legend.position = 'none') +
    facet_wrap(~nsoiltypes, labeller = label_both) +
    labs(x = 'number of groups', y = 'probability (%)')
  
  plot_modularity =
    result %>% 
    mutate(ngroups = factor(number_of_groups)) %>%
    ggplot(aes(ngroups, modularity, fill = ngroups)) + 
    geom_boxplot() +
    facet_wrap(~nsoiltypes, labeller = label_both) +
    theme(legend.position = 'none') +
    labs(x = 'number of groups')
  
  gridExtra::grid.arrange(plot_ngroups, plot_modularity, nrow = 1 )
  
}

if(do.plots2){
  
  theme_set(theme_bw())
  theme_update(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  data = 
    readRDS('~/SpatialNiche/Data/20210720/simulation_analysis.rds')
  
  null =
    data %>%
    select(-c(name, group)) %>% 
    unique() %>% 
    filter(
      theta == 1,
      nsoiltypes == 2
    ) %>%
    group_by(
      weighted, 
      algorithm, 
      d_cutoff,
      rangepar,
      quadrat_length,
    ) %>%
    summarize(
      mu = mean(modularity),
      sigma = sd(modularity),
      .groups = 'drop'
    )
  
  summary = 
    data %>%
    select(-c(name, group)) %>% 
    unique() %>%
    filter(modularity > 0) %>%
    group_by(algorithm, nsoiltypes, d_cutoff, rangepar, quadrat_length, theta, weighted) %>%
    summarize(
      ngroups = mean(number_of_groups), 
      sengroups = sd(number_of_groups) / sqrt(n()), 
      mod = mean(modularity), 
      senmod = sd(modularity) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    left_join(
      null,
      by = c('algorithm', 'weighted', 'quadrat_length', 'rangepar', 'd_cutoff')
    ) %>%
    mutate(
      effect_size = (mod - mu) / sigma,
      seffect_size = senmod / sigma
    ) %>%
    mutate(
      theta = factor(theta), 
      nsoiltypes = factor(nsoiltypes)
    ) %>%
    filter(
      nsoiltypes != 1,
      d_cutoff == 20,
      theta == 2
    )
  
  plot_ngroups = 
    summary %>%
    ggplot() +
    geom_col(
      aes(
        nsoiltypes, 
        ngroups, 
        group = theta, 
        fill = theta
      ),
      position = position_dodge(width = .9)
    ) +
    geom_errorbar(
      aes(
        x = nsoiltypes,
        ymin = ngroups - 0 * sengroups, 
        ymax = ngroups + 2 * sengroups,
        group = theta
      ), 
      position = position_dodge(width = .9),
      width = .2,
      color = 'darkgrey'
    ) +
    facet_grid(weighted ~ algorithm, labeller = label_both) +
    labs(x = 'true number of soil types', y = 'inferred number of sniches') +
    ggtitle('Inferred number of sniches')
  
  plot_modularity = 
    summary %>%
    ggplot() +
    geom_col(
      aes(
        nsoiltypes, 
        mod, 
        group = theta, 
        fill = theta
      ),
      position = position_dodge(width = .9)
    ) +
    geom_errorbar(
      aes(
        x = nsoiltypes,
        ymin = mod - 0 * senmod, 
        ymax = mod + 2 * senmod,
        group = theta
      ), 
      position = position_dodge(width = .9),
      width = .25,
      color = 'darkgrey'
    ) +
    facet_grid(weighted ~ algorithm, labeller = label_both) +
    labs(x = 'true number of soil types', y = 'modularity') +
    ggtitle('Modularity')
  
  plot_effect = 
    summary %>%
    ggplot() +
    geom_col(
      aes(
        nsoiltypes, 
        effect_size, 
        group = theta, 
        fill = theta
      ),
      position = position_dodge(width = .9)
    ) +
    geom_errorbar(
      aes(
        x = nsoiltypes,
        ymin = effect_size - 0 * seffect_size, 
        ymax = effect_size + 2 * seffect_size,
        group = theta
      ), 
      position = position_dodge(width = .9),
      width = .25,
      color = 'darkgrey'
    ) +
    facet_grid(weighted ~ algorithm, labeller = label_both) +
    labs(x = 'true number of soil types', y = '(modularity - mean(null)) / sd(null)') +
    ggtitle('Modularity - Effect size')
  
  gridExtra::grid.arrange(plot_ngroups, plot_modularity, plot_effect, nrow = 1)
}

if(do.envelope.test){
  theme_set(theme_bw())
  theme_update(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  data = 
    readRDS('~/SpatialNiche/Data/20210720/simulation_analysis.rds')
  
  parms = 
    expand_grid(
      nspecies = 105,
      cen = seq(0, 100, 10),
      nst = sort(setdiff(unique(data$nsoiltypes), 1)),
      specialization = setdiff(sort(unique(data$theta)), 1)
    )
  
  res =
    parms %>%
    pmap_dfr(
      .f = function(nspecies, cen, nst, specialization){
        subdat = 
          data %>%
          filter(
            theta == specialization,
            algorithm == 'louvain',
            census == cen,
            weighted == TRUE,
            d_cutoff == 10,
            nsoiltypes == nst
          ) %>%
          mutate(
            sp = factor(name, levels = 1:nspecies),
            group = as.numeric(group)
          )
        
        species_group = 
          rep(seq(nst), each = nspecies / nst)
        
        delta = nspecies - length(species_group)
        if(delta > 0){
          species_group = c(rep(1, delta), species_group)
        }
        
        dtf =
          subdat %>%
          left_join(
            tibble(
              sp = as.factor(seq(nspecies)),
              true_group = species_group
            ), 
            by = 'sp'
          )
        
        true_nst = 
          dtf %>% 
          pull(true_group) %>% 
          unique %>% 
          length
        
        inferred_nst = 
          dtf %>% 
          pull(group) %>% 
          unique %>% 
          length
        
        dtf %>% 
          mutate(
            true_nst = true_nst,
            inferred_nst = inferred_nst
          ) %>%
          return()
      }
  )
  
  theta_plot = 5
  census_plot = 0
  
  plot_nst = 
    res %>%
    group_by(
      algorithm, 
      weighted, 
      d_cutoff, 
      nsoiltypes, 
      theta, 
      census
    ) %>%
    summarize(
      true_nst = length(unique(true_group)),
      inferred_nst = length(unique(group)),
      .groups = 'drop'
    ) %>%
    group_by(
      algorithm, 
      weighted, 
      d_cutoff, 
      nsoiltypes, 
      theta
    ) %>%
    summarize(
      mean = mean(inferred_nst),
      sd = sd(inferred_nst) / sqrt(n()),
      mean_true = mean(true_nst),
      .groups = 'drop'
    ) %>%
    mutate(true = factor(round(mean_true)), inferred = mean) %>%
    ggplot(aes(true, inferred)) +
    geom_errorbar(aes(x = true, ymin = mean - 2 * sd, ymax = mean + 2 * sd)) +
    geom_col(fill = 'plum4') +
    facet_wrap(~ theta, nrow = 1)
  
  plot = 
    res %>%
    filter(theta == theta_plot, nsoiltypes >= 10, census < 100) %>%
    # filter(census == census_plot) %>% 
    mutate(
      group = factor(group), 
      true_group = factor(true_group)
    ) %>%
    ggplot(aes(group, true_group)) +
    geom_violin(aes(group, true_group, group = group)) + 
    theme(legend.position = 'none') +
    # geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_jitter(height = .1, width = .1) +
    labs(x = 'inferred', y = 'true') +
    facet_wrap(nsoiltypes ~ census, scales = 'free', nrow = 4, labeller = label_both) +
    # facet_grid(nsoiltypes ~ census, scales = 'free') +
    ggtitle(paste('theta = ', theta_plot))
    # ggtitle(paste('census =', census_plot))
  
  plot %>% show
  
  plot_summary = 
    res %>% 
    group_by(algorithm, weighted, d_cutoff, nsoiltypes, theta, sp) %>%
    count(group, true_group) %>% 
    ungroup %>%
    filter(theta > 1) %>%
    mutate(
      true_group = factor(true_group),
      group = factor(group)
    ) %>%
    ggplot(aes(group, true_group, size = n)) +
    geom_point()+
    theme(legend.position = 'none') +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    labs(x = 'inferred', y = 'true') +
    facet_grid(theta ~ nsoiltypes, scales = 'free', labeller = label_both)
  
  obs = 
    res %>% 
    group_by(
      algorithm, 
      weighted, 
      d_cutoff, 
      nsoiltypes, 
      theta, 
      group,
      census
    ) %>%
    summarize(
      mean_diff = mean(dist(true_group)),
      sd_diff = sd(dist(true_group)),
      .groups = 'drop'
    ) %>%
    mutate(null = 0)
  
  null = 
    1:100 %>%
    map_dfr(
    .f = function(k){
      res %>% 
        group_by(
          algorithm, 
          weighted, 
          d_cutoff, 
          nsoiltypes, 
          theta,
          census
        ) %>%
        mutate(
          null_true_group = sample(true_group)
        ) %>% 
        ungroup %>%
        group_by(
          algorithm, 
          weighted, 
          d_cutoff, 
          nsoiltypes, 
          theta,
          group,
          census
        ) %>%
        summarize(
          mean_diff = mean(dist(null_true_group)),
          sd_diff = sd(dist(null_true_group)),
          .groups = 'drop'
        ) %>%
        mutate(null = k) %>%
        return()
    }
  )
  
  bar = 
    obs %>%
    bind_rows(null)
  
  foobar = 
    bar %>%
    group_by(
      algorithm, 
      weighted, 
      d_cutoff, 
      nsoiltypes, 
      theta
    ) %>%
    summarize(
      pval = 
        t.test(
          mean_diff[null == 0], 
          mean_diff[null > 0],
          alternative = 'less'
        )$p.value,
      .groups = 'drop'
    )
  
  plot_obs = 
    obs %>%
    left_join(foobar) %>%
    mutate(
      theta = factor(theta),
      nsoiltypes = factor(nsoiltypes)
    ) %>%
    ggplot(aes(nsoiltypes, mean_diff)) +
    geom_boxplot() +
    facet_wrap(~theta) +
    geom_text(
      aes(
        nsoiltypes, 
        rep(7.8, length(nsoiltypes)), 
        label = ifelse(pval <= .05, '*', '')
      ),
      color = 'red',
      size = 10
    ) +
    coord_cartesian(ylim = c(0, 8))
  
  plot_null = 
    null %>%
    mutate(
      theta = factor(theta),
      nsoiltypes = factor(nsoiltypes)
    ) %>%
    ggplot(aes(nsoiltypes, mean_diff)) +
    geom_boxplot() +
    facet_wrap(~theta) +
    coord_cartesian(ylim = c(0, 8))
    
  gridExtra::grid.arrange(plot_obs, plot_null, nrow = 1)
  
  plot1 = 
    dtf %>%
    ggplot(aes(group, sp, color = difference)) +
    geom_point() +
    ggtitle(
      paste('true #clus: ', true_nst, '   inferred #clus: ', inferred_nst)
    ) +
    labs(color = 'diff') +
    scale_x_continuous(breaks = seq(inferred_nst))
  
  plot2 = 
    dtf %>%
    count(group, true_group) %>%
    ggplot(aes(group, true_group, fill = n)) +
    geom_raster() +
    ggtitle(
      paste('true #clus: ', true_nst, '   inferred #clus: ', inferred_nst)
    ) +
    scale_fill_gradientn(colors = terrain.colors(15))  +
    scale_x_continuous(breaks = seq(inferred_nst)) +
    scale_y_continuous(breaks = seq(true_nst))
  
  # gridExtra::grid.arrange(plot1, plot2, nrow = 1)
  
  
}
