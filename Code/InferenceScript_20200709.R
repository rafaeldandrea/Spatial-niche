## ========== Load libraries ==============
library(tidyverse)

## ============= Directive parameters ===========
perform_inference_on_consensus_data = 0
perform_inference_on_averaged_cells = 0
perform_gap_inference_on_averaged_cells = 0
perform_kkmeans_inference_on_averaged_cells = 0
perform_inference_on_averaged_cells_save_all_NST_results = 0
do_plot_SD_by_numclusters = 0
do_plot_consensus_communities = 0
do_plot_true_matrices = 0
do_plot_variance_over_time = 0
calculate_likelihood = 0
do_plot_inferred_landscapes = 0
do_gap_analysis = 0
do_kkmeans_analysis = 0
analyze_inference_on_averaged_cells_save_all_NST_results = 1



## ================ Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc <- !(user <- Sys.info()["user"]) %in% c("rafael", "wangdz", "rdandrea-dw9")

## ========== Load functions =========
loadfuns <- {
  lapply(
    list.files(
      path = switch(
        match(user, c("rafael", "wangdz", "rdandrea", "dw9", "rdandrea-dw9")),
        "~/R_Functions",
        "/Users/wangdz/Documents/UIUC/Odywer/R_Functions",
        "~/R_Functions",
        "/home/dw9/R_Functions",
        "/home/rdandrea/R_Functions"
      ),
      pattern = "[.][R]$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ),
    source
  )
}

## ========= Set working directory ===========
datadir = 
  switch (
    5,
    '~/SpatialNiche/Data/20200618-20scen-NST300S300-H-betanoise-biotime1e5/',
    '~/SpatialNiche/Data/20200620-20scen-NST20S200-H-betanoise-biotime1e5',
    '~/SpatialNiche/Data/20200622-20scen-NST20S200-H-betanoise-biotime1e5',
    '~/SpatialNiche/Data/20200623-20scen-A70kNST20S220-H-betanoise-biotime2e5',
    '~/SpatialNiche/Data/20200709-80scen-A20kS220-H-betanoise-biotime2e5'
  )

setwd(datadir)

## ========== Functions =====================
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Square_Grid <- {
  function(
    landscape_length = NULL, 
    landscape_area = NULL, 
    cell_length = NULL,
    cell_area = NULL
  ) {
    
    if(!is.null(landscape_area) & is.null(landscape_length)){
      A = landscape_area
    } else if(is.null(landscape_area) & !is.null(landscape_length)){
      A = landscape_length ^ 2
    } else{ 
      stop('Must provide either landscape area or landscape length')
    }
    
    if(!is.null(cell_area) & is.null(cell_length)){
      NC = A / cell_area
    } else if(is.null(cell_area) & !is.null(cell_length)){
      NC = A / cell_length ^ 2
    } else{ 
      stop('Must provide either cell area or cell length')
    } 
    
    cellgrid <-
      as_tibble(pos2coord(pos = seq(A), dim.mat = c(sqrt(A), sqrt(A)))) %>%
      mutate(
        cx = 1 + (V1 - 1) %/% (max(V1) / sqrt(NC)),
        cy = 1 + (V2 - 1) %/% (max(V2) / sqrt(NC)),
        cellname = paste0(
          formatC(cx, width = nchar(sqrt(NC)), flag = "0"),
          formatC(cy, width = nchar(sqrt(NC)), flag = "0")
        ),
        cell = match(cellname, unique(cellname))
      ) %>%
      pull(cell)
    return(cellgrid)
  } 
}

Infer_Affinity_Matrix <- {
  function(
    species_map, 
    cell_length, 
    min_num_clusters, 
    max_num_clusters,
    num_clusters = NULL,
    index = 'sdindex'
  ){
    if(!is.null(num_clusters)){
      min_num_clusters = max_num_clusters = num_clusters
    }
    
    grid = Square_Grid(landscape_area = length(species_map), cell_length = cell_length)
    
    cell_by_species = table(grid, species_map)
    
    dis_sites <- vegan::vegdist(cell_by_species, method = "jaccard")
    
    clustresult <- 
      NbClust::NbClust(
        data = cell_by_species, 
        diss = dis_sites, 
        distance = NULL, 
        min.nc = min_num_clusters, 
        max.nc = min(max_num_clusters, length(species_map)/cell_length ^ 2 - 1),
        method = 'ward.D2', 
        index = index
      )
    
    inferred_landscape = as.numeric(clustresult$Best.partition[grid])
    
    abundance_by_soiltype = table(species_map, inferred_landscape, dnn = c('species', 'soil_type'))
    
    inferred_matrix = 
      as_tibble(abundance_by_soiltype / rowSums(abundance_by_soiltype)) %>%
      rename(affinity = n) %>%
      mutate(species = as.numeric(species))
    
    soiltype_order = NULL
    for(sp in sort(unique(species_map))){
      soiltypes_sorted =
        inferred_matrix %>%
        filter(species == sp) %>%
        arrange(desc(affinity)) %>%
        pull(soil_type) %>%
        unique
      soiltype_order = c(soiltype_order, soiltypes_sorted[!soiltypes_sorted %in% soiltype_order][1])
    }
    soiltype_order = c(soiltype_order, setdiff(seq(clustresult$Best.nc['Number_clusters']), soiltype_order))
    
    clustresult$grid = 
      grid
    
    clustresult$inferred_landscape = 
      match(inferred_landscape, as.numeric(soiltype_order))
    
    clustresult$abundance_by_soiltype = 
      table(species_map, clustresult$inferred_landscape, dnn = c('species', 'soil_type'))
    
    clustresult$inferred_matrix =
      inferred_matrix %>%
      mutate(soil_type = factor(match(soil_type, soiltype_order))) %>%
      arrange(species, soil_type)
    
    return(clustresult) 
  }
  
} 

Infer_Affinity_Matrix_from_CellTable <- {
  function(
    celltable,
    min_num_clusters, 
    max_num_clusters,
    num_clusters = NULL,
    index = 'sdindex'
  ){
    if(!is.null(num_clusters)){
      min_num_clusters = max_num_clusters = num_clusters
    }
    
    dis_sites <- vegan::vegdist(celltable, method = "jaccard")
    
    clustresult <- 
      NbClust::NbClust(
        data = celltable, 
        diss = dis_sites, 
        distance = NULL, 
        min.nc = min_num_clusters, 
        max.nc = min(max_num_clusters, nrow(celltable) - 1),
        method = 'ward.D2', 
        index = index
      )
    
    inferred_landscape = 
      tibble(
        gridcell = seq(nrow(celltable)),
        soil_type = clustresult$Best.partition
      )
    
    foo = as_tibble(celltable)
    names(foo) = c('gridcell', 'species', 'abundance')
    foo = 
      foo %>%
      mutate(
        gridcell = as.numeric(gridcell),
        species = as.numeric(species)
      )
    
    abundance_by_soiltype = 
      foo %>%
      left_join(inferred_landscape, by = 'gridcell') %>%
      group_by(soil_type, species) %>%
      summarize(abundance = sum(abundance)) %>%
      ungroup() %>%
      pivot_wider(names_from = 'species', values_from = 'abundance') %>%
      select(-soil_type) %>% 
      t()
    
    inferred_matrix = 
      as_tibble(abundance_by_soiltype / rowSums(abundance_by_soiltype)) %>%
      bind_cols(species = seq((ncol(celltable)))) %>%
      pivot_longer(
        cols = -species,
        names_to = 'soil_type',
        values_to = 'affinity'
      )
    
    soiltype_order = NULL
    for(sp in seq(ncol(celltable))){
      soiltypes_sorted =
        inferred_matrix %>%
        filter(species == sp) %>%
        arrange(desc(affinity)) %>%
        pull(soil_type) %>%
        unique
      soiltype_order = c(soiltype_order, soiltypes_sorted[!soiltypes_sorted %in% soiltype_order][1])
    }
    soiltype_order = c(soiltype_order, setdiff(paste0('V', seq(clustresult$Best.nc['Number_clusters'])), soiltype_order))
    soiltype_order = soiltype_order[!is.na(soiltype_order)]
    
    clustresult$inferred_matrix =
      inferred_matrix %>%
      mutate(soil_type = factor(match(soil_type, soiltype_order))) %>%
      arrange(species, soil_type)
    
    soiltype_order_numeric = sapply(strsplit(soiltype_order, split = 'V'), function(l) as.numeric(l[2]))
    
    clustresult$inferred_landscape = 
      inferred_landscape %>%
      mutate(soil_type = match(soil_type, soiltype_order_numeric))
    
    clustresult$abundance_by_soiltype = 
      abundance_by_soiltype[, soiltype_order_numeric]
    
    return(clustresult) 
  }
  
} 

## ========== Scenarios =====================
scenarios = 
  crossing(
    cell_length = 
      switch(
        3,
        c(3, 4, 6, 8, 9, 12, 16, 18),
        c(3, 4, 6, 8, 11, 12),
        c(2, 3, 4, 6, 8, 9, 12)
      ),
    scenario = 1:80
  ) %>%
  rowid_to_column(var = 'jobid')

## filter to niche index > .5
trimmed_scenarios = 
  scenarios %>% 
  select(-jobid) %>%
  filter(1 + (scenario - 1) %% 20 > 10) %>%
  rowid_to_column(var = 'jobid')
## ==========================================



if(perform_inference_on_consensus_data){
  ## ========= Read job id, data ==================
  ind = make_index(1)
  scen = scenarios[ind, ]
  list2env(scen, envir = environment())
  
  filename = paste0('scenario_', scenario, '.RData')
  
  dat = get(load(filename))
  community = matrix(unlist(dat$evolving_community), 144 ^ 2)
  
  
  ## ========== Set censuses to use in analysis ==========
  threshold_census = ifelse(scenario < 15, 24, 3)
  if(scenario == 20) threshold_census = ncol(community)
  column_index = threshold_census:ncol(community)
  
  
  ## =========== Perform inference ======================
  results = list()
  
  the_num_clusters = 
    switch(
      2,
      c(2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 150, 200, 250),
      seq(2, 30, by = 2)
    )
  
  for(num_clusters in the_num_clusters){
    
    for(census_block in c(0, 1:14)) {
      
      if(!hpcc) 
        writeLines(
          paste('number of clusters =', num_clusters,
                'census block =', census_block)
        )
      
      if(census_block > 0){
        selected_columns = column_index[1 + (column_index - threshold_census) %/% 7 == census_block]
      } else{
        selected_columns = column_index
      }
      
      if(length(selected_columns) > 0){
        
        if(length(selected_columns) > 1){
          consensus_community = 
              apply(community[, selected_columns], 1, function(speciesid) getmode(speciesid))
        }
       
        if(length(selected_columns) == 1){
          consensus_community = community[, selected_columns]
        } 
        
        res =
          Infer_Affinity_Matrix(
            species_map = consensus_community,
            cell_length = cell_length,
            num_clusters = num_clusters
          )
        
        res$input_parameters = 
          list(
            scenario = scenario,
            cell_length = cell_length,
            num_clusters = num_clusters,
            census_block = census_block,
            jobid = jobid,
            selected_censuses = selected_columns
          )
        
        results[[length(results) + 1]] = res
      }
    }
    
    if(hpcc) save(results, file = paste0('analysis_jobID_', jobid, '.RData'))
    
  } 
  
  
}

if(perform_inference_on_averaged_cells){
  ## ========= Read job id, data ==================
  ind = make_index(19)
  scen = trimmed_scenarios[ind, ]
  list2env(scen, envir = environment())
  
  filename = paste0('scenario_', scenario, '.RData')
  
  dat = get(load(filename))
  community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
  grid = Square_Grid(landscape_area = nrow(community), cell_length = cell_length)
  
  ## ========== Set censuses to use in analysis ==========
  # threshold_census = ifelse(scenario < 15, 17, 3)
  # if(scenario == 20) threshold_census = ncol(community)
  threshold_census = 30
  column_index = threshold_census:ncol(community)
  
  
  ## =========== Perform inference ======================
  results = list()
  
  # the_num_clusters = 
  #   switch(
  #     2,
  #     c(2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 150, 200, 250),
  #     seq(2, 30, by = 2),
  #     seq(3, 1.5 * ncol(dat$C), by = 2)
  #   )
  # 
  # for(num_clusters in the_num_clusters){
    
    for(census_block in 0) {
      
      # if(!hpcc) 
      #   writeLines(
      #     paste('number of clusters =', num_clusters,
      #           'census block =', census_block)
      #   )
      
      if(census_block > 0){
        selected_columns = column_index[1 + (column_index - threshold_census) %/% 7 == census_block]
      } else{
        selected_columns = column_index
      }
      
      if(length(selected_columns) > 0){
        
        averaged_community = 
          cbind(
            grid = grid,
            community[, selected_columns]
          )
        colnames(averaged_community) = c('grid', paste0('census_', formatC(selected_columns, digits = 1, flag = '0')))
        
        averaged_community = as_tibble(averaged_community)
        
        avg_community = 
          averaged_community %>%
          pivot_longer(cols = -grid, names_to = 'census', values_to = 'species')
        
        celltable = with(avg_community, table(grid, species))
        
        
        # res =
        #   Infer_Affinity_Matrix_from_CellTable(
        #     celltable = celltable,
        #     num_clusters = num_clusters,
        #     index = 'sdindex'
        #   )
        
        res =
          Infer_Affinity_Matrix_from_CellTable(
            celltable = celltable,
            min_num_clusters = 2, 
            max_num_clusters = max(10, 1.5 * ncol(dat$C)),
            index = 'sdindex'
          )
        
        res$input_parameters = 
          list(
            scenario = scenario,
            cell_length = cell_length,
            # num_clusters = num_clusters,
            census_block = census_block,
            jobid = jobid,
            selected_censuses = selected_columns
          )
        
        results[[length(results) + 1]] = res
      }
    }
    
    if(hpcc) 
      save(
        results, 
        file = paste0('analysis_avg_community_scneario_', scenario, '_celllength_', cell_length,'.RData')
      )
    
  } 

if(perform_inference_on_averaged_cells_save_all_NST_results){
  ## ========= Read job id, data ==================
  ind = make_index(19)
  scen = trimmed_scenarios[ind, ]
  list2env(scen, envir = environment())
  
  filename = paste0('scenario_', scenario, '.RData')
  
  dat = get(load(filename))
  community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
  grid = Square_Grid(landscape_area = nrow(community), cell_length = cell_length)
  
  ## ========== Set censuses to use in analysis ==========
  # threshold_census = ifelse(scenario < 15, 17, 3)
  # if(scenario == 20) threshold_census = ncol(community)
  threshold_census = 30
  column_index = threshold_census:ncol(community)
  
  
  ## =========== Perform inference ======================
  results = list()
  
  the_num_clusters =
    switch(
      3,
      c(2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 150, 200, 250),
      seq(2, 30, by = 2),
      seq(3, 1.5 * ncol(dat$C), by = 2)
    )

  for(num_clusters in the_num_clusters){
  
    for(census_block in 0) {
      
      if(!hpcc)
        writeLines(
          paste('number of clusters =', num_clusters,
                'census block =', census_block)
        )
      
      if(census_block > 0){
        selected_columns = column_index[1 + (column_index - threshold_census) %/% 7 == census_block]
      } else{
        selected_columns = column_index
      }
      
      if(length(selected_columns) > 0){
        
        averaged_community = 
          cbind(
            grid = grid,
            community[, selected_columns]
          )
        colnames(averaged_community) = c('grid', paste0('census_', formatC(selected_columns, digits = 1, flag = '0')))
        
        averaged_community = as_tibble(averaged_community)
        
        avg_community = 
          averaged_community %>%
          pivot_longer(cols = -grid, names_to = 'census', values_to = 'species')
        
        celltable = with(avg_community, table(grid, species))
        
        
        res =
          Infer_Affinity_Matrix_from_CellTable(
            celltable = celltable,
            num_clusters = num_clusters,
            index = 'sdindex'
          )
        
        res$true_C = dat$C
        
        res$input_parameters = 
          list(
            scenario = scenario,
            cell_length = cell_length,
            num_clusters = num_clusters,
            census_block = census_block,
            jobid = jobid,
            selected_censuses = selected_columns
          )
        
        results[[length(results) + 1]] = res
        
        if(hpcc) 
          save(
            results, 
            file = paste0('analysis_allcensuses_scenario_', scenario, '_celllength_', cell_length,'.RData')
          )
      }
    }
  }
  
  
  
} 

if(perform_gap_inference_on_averaged_cells){
  
  ## ========= Read job id, data ==================
  ind = make_index(19)
  scen = scenarios[ind, ]
  list2env(scen, envir = environment())
  
  dat = get(load(paste0('scenario_', scenario, '.RData')))
  community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
  
  
  ## ========= file name to save results =========
  filename = paste0('gap_analysis_scenario_', scenario, '_cell_length_', cell_length,'.RData')
  
  
  ## ========= fixed parameters ============
  num_nulls = 10
  max_nc = 30
  
  ## ========== create grid according to cell size ===============================
  grid = Square_Grid(landscape_area = nrow(community), cell_length = cell_length)
  
  ## ========== Set censuses to use in analysis ==========
  threshold_census = 17
  column_index = threshold_census:ncol(community)
  
  
  ## =========== Perform inference ======================
  for(census_block in 0:12) {
      
    if(!hpcc)  writeLines(paste('census block =', census_block))
    
    if(census_block > 0){
      selected_columns = column_index[1 + (column_index - threshold_census) %/% 7 == census_block]
    } else{
      selected_columns = column_index
    }
      
    averaged_community = 
      cbind(
        grid = grid,
        community[, selected_columns]
      )
    colnames(averaged_community) = c('grid', paste0('census_', formatC(selected_columns, digits = 1, flag = '0')))
    
    averaged_community = as_tibble(averaged_community)
    
    avg_community = 
      averaged_community %>%
      pivot_longer(cols = -grid, names_to = 'census', values_to = 'species')
    
    indices = {
      crossing(
        null = seq(0, num_nulls), 
        num_clusters = seq(2, max_nc, by = 2)
      ) %>%
      ddply(.(null, num_clusters), function(df){
        
        seed = df$null
        num_clusters = df$num_clusters
        
        writeLines(paste('null = ', seed,'   nc = ', num_clusters))
        
        set.seed(seed)
        
        if(seed == 0){
          celltable = with(avg_community, table(grid, species)) / length(selected_columns)
        } 
        
        if(seed != 0){
          null_avg_community = 
            avg_community %>% 
            group_by(census) %>% 
            mutate(null_species_map = sample(species)) %>% 
            ungroup
          
          celltable = with(null_avg_community, table(grid, null_species_map)) / length(selected_columns)
        } 
        
        dis_cells = vegan::vegdist(celltable, method = "jaccard")
        
        sdindex = 
          NbClust::NbClust(
            data = celltable, 
            diss = dis_cells, 
            distance = NULL, 
            min.nc = num_clusters, 
            max.nc = num_clusters,
            method = 'ward.D2', 
            index = 'sdindex'
          )$All.index
        
        return(sdindex)
      }) %>%
      as_tibble %>%
      rename(sdindex = V1) %>%
      mutate(
        niche_index_scenario = scenario,
        cell_length = cell_length,
        census_block = census_block
      )
    }
    
    if(file.exists(filename)){
      foo = get(load(filename))
      results = 
        bind_rows(
          foo,
          indices
        )
    } else{
      results = indices
    }
    
    if(hpcc) save(results, file = filename)
  }
}

if(perform_kkmeans_inference_on_averaged_cells){
  
  ## ========= Read job id, data ==================
  ind = make_index(79)
  scen = scenarios[ind, ]
  list2env(scen, envir = environment())
  
  dat = get(load(paste0('scenario_', scenario, '.RData')))
  community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
  
  
  ## ========= file name to save results =========
  filename = paste0('kkmeans_inference_scenario_', scenario, '_cell_length_', cell_length,'.RData')
  
  
  ## ========= fixed parameters ============
  num_nulls = 10
  max_nc = 30
  
  ## ========== create grid according to cell size ===============================
  grid = Square_Grid(landscape_area = nrow(community), cell_length = cell_length)
  
  ## ========== Set censuses to use in analysis ==========
  threshold_census = 17
  column_index = threshold_census:ncol(community)
  
  
  ## =========== Perform inference ======================
  for(census_block in 0:12) {
    
    if(!hpcc)  writeLines(paste('census block =', census_block))
    
    if(census_block > 0){
      selected_columns = column_index[1 + (column_index - threshold_census) %/% 7 == census_block]
    } else{
      selected_columns = column_index
    }
    
    averaged_community = 
      cbind(
        grid = grid,
        community[, selected_columns]
      )
    colnames(averaged_community) = c('grid', paste0('census_', formatC(selected_columns, digits = 1, flag = '0')))
    
    averaged_community = as_tibble(averaged_community)
    
    avg_community = 
      averaged_community %>%
      pivot_longer(cols = -grid, names_to = 'census', values_to = 'species')
    
    for(null in seq(0, num_nulls)){
      
      if(null == 0){
        celltable = with(avg_community, table(grid, species)) / length(selected_columns)
      } 
      
      if(null != 0){
        set.seed(null)
        
        null_avg_community = 
          avg_community %>% 
          group_by(census) %>% 
          mutate(null_species_map = sample(species)) %>% 
          ungroup
        
        celltable = with(null_avg_community, table(grid, null_species_map)) / length(selected_columns)
      } 
      
      for(num_clusters in seq(2, max_nc, by = 2)){
        
        writeLines(paste('null = ', null,'   nc = ', num_clusters))
        
        dispersion = NULL
        while(is.null(dispersion) | class(dispersion) == 'try-error'){
          dispersion = try(kernlab::withinss(kernlab::kkmeans(unclass(celltable), centers = num_clusters)), silent = TRUE)
        }
        total_withinss = sum(dispersion)
        
        indices = 
          tibble(
            niche_index_scenario = scenario,
            cell_length = cell_length,
            census_block = census_block,
            num_clusters = num_clusters,
            null = null,
            total_withinss = total_withinss
          )
        
        if(file.exists(filename)){
          foo = get(load(filename))
          results = 
            bind_rows(
              foo,
              indices
            )
        } else{
          results = indices
        }
        
        if(hpcc) save(results, file = filename)
        
      }
    }
  }
}

if(calculate_likelihood){
  
  dat = get(load('analysis_avg_community_jobID_10.RData'))
  
  grid = Square_Grid(
    
  )
  
  loglik = 
    calculateLikelihood(
      species_map = community,
      inferred_landscape,
      affinity_matrix = estC
    )
}

if(do_plot_SD_by_numclusters){
  dtf = NULL
  for(niscen in 1:20){
    jobids = 
      scenarios %>% 
      filter(scenario == niscen) %>%
      pull(jobid)
    
    foo = list()
    for(j in jobids){
      b = 
        try(
          # get(load(paste0('analysis_jobID_', j, '.RData'))),
          get(load(paste0('analysis_avg_community_jobID_', j, '.RData'))),
          silent = TRUE
        )
      if(class(b) != 'try-error')
        foo = c(foo, b)
    }
    
    bar = 
      t(
        sapply(foo, function(l){
          pars = l$input_parameters
          return(
            c(
              scenario = niscen,
              jobid = pars$jobid,
              cell_length = pars$cell_length,
              num_clusters = pars$num_clusters,
              census_block = pars$census_block,
              sdindex = as.numeric(l$Best.nc['Value_Index'])
            )
          )
        })
      ) %>%
      as_tibble
    
    dtf = rbind(dtf, bar)
  }
  
  res = 
    dtf %>%
    filter(census_block != 0) %>%
    group_by(scenario, cell_length, num_clusters) %>%
    summarize(sd_mean = mean(sdindex), .groups = 'drop')
  
  res_allcensuses = 
    dtf %>%
    filter(census_block == 0) %>%
    group_by(scenario, cell_length, num_clusters) %>%
    summarize(sd_allcensuses = mean(sdindex), .groups = 'drop')
  
  res = 
    res %>%
    left_join(res_allcensuses)
  
  plot_SD_by_numclusters = 
    res %>% 
    filter(cell_length == 6) %>% 
    mutate(cell_length = factor(cell_length)) %>% 
    ggplot(aes(num_clusters, sd_mean, group = cell_length, color = cell_length)) + 
    geom_line() + 
    geom_point() + 
    facet_wrap(~ scenario, scale = 'free') + 
   # scale_x_log10() +
    scale_y_log10()
  
  
  gridExtra::grid.arrange(plot_SD_by_numclusters)  
}

if(do_plot_consensus_communities){
  
  dtf = NULL
  for(scenario in 1:20){
    filename = paste0('scenario_', scenario, '.RData')
    
    dat = get(load(filename))
    community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
    
    
    ## ========== Set censuses to use in analysis ==========
    # threshold_census = ifelse(scenario < 15, 17, 3)
    # if(scenario == 20) threshold_census = ncol(community)
    threshold_census = 17
    column_index = threshold_census:ncol(community)
    
   # if(scenario < 20){
      consensus_community = 
        apply(community[, column_index], 1, function(speciesid) getmode(speciesid))
    #} else {
     # consensus_community = community[, column_index]
    #}
    
    dtf = cbind(dtf, consensus_community)
    
  }
  
  dtf = as_tibble(dtf)
  names(dtf) = paste0('scenario_', formatC(1:20, digits = 1, flag = '0'))
  
  coords = 
    as_tibble(
      pos2coord(pos = 1:length(dat$Landscape), dim.mat = c(sqrt(length(dat$Landscape)), sqrt(length(dat$Landscape))))
    )
  names(coords) = c('x', 'y')
  
  foo = 
    dtf %>%
    rowid_to_column(var = 'site') %>%
    bind_cols(coords) %>%
    pivot_longer(cols = -c(site, x, y), names_to = 'scenario', values_to = 'species') 
  
  plot_consensus_communities = 
    foo %>%
    ggplot(aes(x, y, fill = species)) +
    geom_tile() +
    facet_wrap(~ scenario) +
    theme(aspect.ratio = 1)
  
  gridExtra::grid.arrange(plot_consensus_communities)
}

if(do_plot_inferred_landscapes){
  
  landscape_length = sqrt(length(get(load('scenario_1.rdata'))$Landscape))
  
  foo = list()
  for(j in 1:20){
    b = 
      try(
        # get(load(paste0('analysis_jobID_', j, '.RData'))),
        get(load(paste0('analysis_avg_community_jobID_', j, '.RData'))),
        silent = TRUE
      )
    if(class(b) != 'try-error')
      foo = c(foo, b)
  }
  
  bar = 
    t(
      sapply(foo, function(l){
        pars = l$input_parameters
        return(
          c(
            jobid = pars$jobid,
            cell_length = pars$cell_length,
            num_clusters = pars$num_clusters,
            census_block = pars$census_block,
            sdindex = as.numeric(l$Best.nc['Value_Index'])
          )
        )
      })
    ) %>%
    as_tibble %>%
    dplyr::rename(scenario = jobid) %>%
    rowid_to_column(var = 'foo_index') %>%
    dplyr::filter(census_block == 0) %>%
    group_by(scenario) %>%
    slice(which.min(sdindex)) %>%
    ungroup
  
  
  foobar = foo[bar$foo_index]
  
  res_inferred_matrix = NULL
  for(l in 1:20){
    res_inferred_matrix = 
      rbind(
        res_inferred_matrix,
          foobar[[l]]$inferred_matrix %>%
          mutate(scenario = l)
        )
  }
  
  plot_inferred_matrix =
    res_inferred_matrix %>%
    ggplot(aes(species, soil_type, fill = affinity)) +
    geom_tile() +
    #theme(aspect.ratio = 1) +
    facet_wrap( ~scenario, scales = 'free')
  
  site_info = 
    as_tibble(
      pos2coord(pos = 1:(landscape_length ^ 2), dim.mat = c(landscape_length, landscape_length))
    ) %>%
    mutate(
      gridcell = 
        Square_Grid(
          landscape_length = landscape_length, 
          cell_length = 3
        )
    )
  names(site_info) = c('x', 'y', 'gridcell')
  
  
  
  res_landscape = NULL
  for(l in 1:20){
    
    site_info =
      site_info %>%
      mutate(
        soil_type = foobar[[l]]$Best.partition[gridcell],
        scenario = l
      ) 
    
    res_landscape = 
      rbind(
        res_landscape,
        site_info
      )
  }
  
  plot_inferred_landscape =
    res_landscape %>%
    ggplot(aes(x, y, fill = soil_type)) +
    geom_tile() +
    theme(aspect.ratio = 1) +
    facet_wrap( ~scenario, scales = 'free')
  
  
  true_C = NULL
  for(ni_scenario in 1:20){
    true_C = 
      rbind(
        true_C,
        as_tibble(get(load(paste0('scenario_', ni_scenario, '.rdata')))$C) %>%
          rowid_to_column(var = 'species') %>%
          pivot_longer(
            cols = -species, 
            names_to = 'soil_type', 
            values_to = 'affinity'
          ) %>%
          mutate(scenario = ni_scenario)
      )
  }
  
  RVs = 
    sapply(1:20, function(ni_scenario){
      
      tC =
        true_C %>% 
        filter(scenario == ni_scenario) %>% 
        select(-scenario) %>%
        pivot_wider(
          names_from = 'soil_type', 
          values_from = 'affinity'
        ) %>% 
        select(-species) %>% 
        as.matrix()
      
      iC =
        res_inferred_matrix %>% 
        filter(scenario == ni_scenario) %>% 
        select(-scenario) %>%
        pivot_wider(
          names_from = 'soil_type', 
          values_from = 'affinity'
        ) %>% 
        select(-species) %>% 
        as.matrix()
      
      RV_coefficient(tC, iC)$RV
      
    })
  
  plot_RV = 
    tibble(
    `niche index` = seq(.05, 1, by = .05),
      `RV coefficient` = RVs
    ) %>%
    ggplot(aes(`niche index`, `RV coefficient`)) +
    geom_line() +
    geom_point() +
    ylim(c(0,1))
  
    
}

if(do_plot_true_matrices){
  
  dtf = NULL
  for(scenario in 1:20){
    filename = paste0('scenario_', scenario, '.RData')
    
    dat = get(load(filename))
    dtf = cbind(dtf, as.numeric(dat$C))
  }
  
  dtf = as_tibble(dtf)
  names(dtf) = paste0('scenario_', formatC(1:20, digits = 1, flag = '0'))
  
  coords = 
    as_tibble(
      pos2coord(pos = 1:4000, dim.mat = c(200, 20))
    )
  names(coords) = c('species', 'soil type')
  
  foo = 
    dtf %>%
    rowid_to_column(var = 'site') %>%
    bind_cols(coords) %>%
    pivot_longer(
      cols = -c(site, species, `soil type`), 
      names_to = 'scenario', 
      values_to = 'affinity'
    ) 
  
  plot_true_matrices = 
    foo %>%
    ggplot(aes(species, `soil type`, fill = affinity)) +
    geom_tile() +
    facet_wrap(~ scenario) +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  gridExtra::grid.arrange(plot_true_matrices)
  
}

if(do_plot_variance_over_time){
  
  dtf = NULL
  for(scenario in 1:20){
    filename = paste0('scenario_', scenario, '.RData')
    
    dat = get(load(filename))
    community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
    var_over_time = apply(community, 2, var)
    
    dtf = cbind(dtf, var_over_time)
    
  }
  
  dtf = as_tibble(dtf)
  names(dtf) = paste0('scenario_', formatC(1:20, digits = 1, flag = '0'))
  dtf =
    dtf %>%
    rowid_to_column(var = 'census') %>%
    pivot_longer(
      cols = -census, 
      names_to = 'scenario', 
      values_to = 'variance'
    )
  
  plot_variance_over_time = 
    dtf %>%
    ggplot(aes(census, variance)) +
    geom_line() +
    facet_wrap(~ scenario) +
    theme(aspect.ratio = 1)
  
  gridExtra::grid.arrange(plot_variance_over_time)
}

if(do_gap_analysis){
  
  chosen_cell_length = 3
  chosen_census_block = 0
  min_numclust = 3
  
  files = list.files(pattern = 'gap_analysis_*')
  
  foo = get(load(files[1]))
  for(char in files[-1]){
    foo = 
      foo %>%
      bind_rows(bar <- get(load(char)))
  }
  
  dat = 
    foo %>%
    group_by(
      niche_index_scenario,
      null, 
      cell_length, 
      census_block
    ) %>%
    mutate(
      standardized_sdindex = -(sdindex - mean(sdindex)) / sd(sdindex)
    ) %>%
    ungroup
  
  results = 
    dat %>%
    filter(null > 0) %>%
    group_by(
      niche_index_scenario, 
      cell_length, 
      census_block, 
      num_clusters
    ) %>%
    summarize(
      mean_null_sdindex = mean(standardized_sdindex), 
      sd_null_sdindex = sd(standardized_sdindex), 
      .groups = 'drop'
    ) %>%
    left_join(
      dat %>%
        filter(null == 0) %>%
        select(-null)
    ) %>%
    mutate(
      gap = standardized_sdindex - mean_null_sdindex,
      sgap = sd_null_sdindex / sqrt(10)
    )
  
  results_for_plot = 
    results %>%
    mutate(
      census_block = factor(census_block),
      cell_length = factor(cell_length)
    ) %>%
    filter(cell_length == chosen_cell_length) %>%
    filter(census_block == chosen_census_block)
  
  plot_sdindex =
    results_for_plot %>%
    ggplot(aes(num_clusters, standardized_sdindex)) +
    geom_errorbar(
      aes(
        x = num_clusters, 
        ymin = mean_null_sdindex - sd_null_sdindex, 
        ymax = mean_null_sdindex + sd_null_sdindex
      ), 
      width = .5,
      color = red
    ) +
    geom_line(aes(num_clusters, mean_null_sdindex), color = red) +
    geom_point(aes(num_clusters, mean_null_sdindex), color = red) +
    geom_line() +
    geom_point() +
    facet_wrap(~ niche_index_scenario, scales = 'free') +
    ggtitle(
      'Standardized SD index of simulated communities = -(index - mean(index)) / sd(index)                 black = data,   red = nulls'
    )
  
  plot_null_sdindex =
    results_for_plot %>%
    ggplot(aes(num_clusters, mean_null_sdindex)) +
    geom_errorbar(
      aes(
        x = num_clusters, 
        ymin = mean_null_sdindex - sd_null_sdindex, 
        ymax = mean_null_sdindex + sd_null_sdindex
      ), 
      width = .5,
      color = red
    ) +
    geom_line(color = red) +
    geom_point(color = red) +
    facet_wrap(~ niche_index_scenario, scales = 'free') +
    ggtitle('Mean standardized SD index of null communities') + 
    ylim(c(-3, 3))
  
  plot_gap =
    results_for_plot %>%
    filter(num_clusters >= min_numclust) %>%
    ggplot(aes(num_clusters, gap)) +
    geom_errorbar(aes(x = num_clusters, ymin = gap - sgap, ymax = gap + sgap), width = .5, color = blue) +
    geom_line(color = blue) +
    geom_point(color = blue) +
    facet_wrap(~ niche_index_scenario, scales = 'free') +
    ggtitle('Gap index of simulated communities: gap = obs - mean(nulls)')
  
  gridExtra::grid.arrange(plot_sdindex, plot_gap, nrow = 1)
  
  best_numclust = 
    results %>%
    filter(cell_length == chosen_cell_length) %>%
    filter(census_block == chosen_census_block) %>%
    filter(num_clusters >= min_numclust) %>% 
    group_by(niche_index_scenario) %>% 
    slice(which.max(gap)) %>% 
    select(cell_length, niche_index_scenario, num_clusters_by_gap = num_clusters) %>%
    left_join(
      results %>%
        filter(cell_length == chosen_cell_length) %>%
        filter(census_block == chosen_census_block) %>%
        filter(num_clusters >= min_numclust) %>% 
        group_by(niche_index_scenario) %>% 
        slice(which.max(standardized_sdindex)) %>% 
        select(cell_length, niche_index_scenario, num_clusters_by_sdindex = num_clusters)
      
    )
  
  plot_best_numclust = 
    best_numclust %>%
    ggplot() +
    geom_line(aes(niche_index_scenario, num_clusters_by_gap), color = blue) +
    geom_point(aes(niche_index_scenario, num_clusters_by_gap), color = blue) +
    geom_line(aes(niche_index_scenario, num_clusters_by_sdindex)) +
    geom_point(aes(niche_index_scenario, num_clusters_by_sdindex)) +
    geom_hline(aes(yintercept = 20), color = green) +
    labs(
      x = 'Niche index scenario',
      y = 'Best number of clusters'
    )
    ggtitle('Best number of clusters.         black = SD index,  blue = gap index.         Range of search = 3 to 30 clusters')
  
  gridExtra::grid.arrange(plot_best_numclust)
  
  result = list()
  for(i in seq(nrow(best_numclust))){
    
    print(i)
    
    scen = best_numclust[i, ]
    scenario = scen$niche_index_scenario
    cell_length = as.numeric(scen$cell_length)
    num_clusters = scen$num_clusters_by_gap
    
    
    dat = get(load(paste0('scenario_', scenario, '.RData')))
    community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
   
    grid = Square_Grid(landscape_area = nrow(community), cell_length = cell_length)
    
    threshold_census = 17
    selected_columns = threshold_census:ncol(community)
    
    averaged_community = 
      cbind(
        grid = grid,
        community[, selected_columns]
      )
    colnames(averaged_community) = c('grid', paste0('census_', formatC(selected_columns, digits = 1, flag = '0')))
    
    averaged_community = as_tibble(averaged_community)
    
    avg_community = 
      averaged_community %>%
      pivot_longer(cols = -grid, names_to = 'census', values_to = 'species')
    
    celltable = with(avg_community, table(grid, species)) / length(selected_columns)
    
    res =
      Infer_Affinity_Matrix_from_CellTable(
        celltable = celltable,
        num_clusters = num_clusters,
        index = 'sdindex'
      )
    
    result = c(result, res)
    
  }
  
  
  
  ## ========== Set censuses to use in analysis ==========
  # threshold_census = ifelse(scenario < 15, 17, 3)
  # if(scenario == 20) threshold_census = ncol(community)
  threshold_census = 17
  column_index = threshold_census:ncol(community)
  
  
}

if(do_kkmeans_analysis){
  
  chosen_cell_length = 3
  chosen_census_block = 0
  min_numclust = 3
  
  files = list.files(pattern = 'kkmeans_inference_*')
  
  foo = get(load(files[1]))
  for(char in files[-1]){
    foo = 
      foo %>%
      bind_rows(bar <- get(load(char)))
  }
  
  dat = 
    foo %>%
    group_by(
      niche_index_scenario,
      null, 
      cell_length, 
      census_block
    ) %>%
    mutate(
      # standardized_index = -(total_withinss - mean(total_withinss)) / sd(total_withinss)
      standardized_index = log(total_withinss)
    ) %>%
    ungroup
  
  results = 
    dat %>%
    filter(null > 0) %>%
    group_by(
      niche_index_scenario, 
      cell_length, 
      census_block, 
      num_clusters
    ) %>%
    summarize(
      mean_null_index = mean(standardized_index), 
      sd_null_index = sd(standardized_index), 
      .groups = 'drop'
    ) %>%
    left_join(
      dat %>%
        filter(null == 0) %>%
        select(-null)
    ) %>%
    mutate(
      gap = standardized_index - mean_null_index,
      sgap = sd_null_index / sqrt(10)
    )
  
  results_for_plot = 
    results %>%
    mutate(
      census_block = factor(census_block),
      cell_length = factor(cell_length)
    ) %>%
    filter(cell_length == chosen_cell_length) %>%
    filter(census_block == chosen_census_block)
  
  plot_index =
    results_for_plot %>%
    ggplot(aes(num_clusters, standardized_index)) +
    geom_errorbar(
      aes(
        x = num_clusters, 
        ymin = mean_null_index - sd_null_index, 
        ymax = mean_null_index + sd_null_index
      ), 
      width = .5,
      color = red
    ) +
    geom_line(aes(num_clusters, mean_null_index), color = red) +
    geom_point(aes(num_clusters, mean_null_index), color = red) +
    geom_line() +
    geom_point() +
    facet_wrap(~ niche_index_scenario, scales = 'free') +
    ggtitle(
      'kmeans index of simulated communities                 black = data,   red = nulls'
    )
  
  plot_null_sdindex =
    results_for_plot %>%
    ggplot(aes(num_clusters, mean_null_sdindex)) +
    geom_errorbar(
      aes(
        x = num_clusters, 
        ymin = mean_null_sdindex - sd_null_sdindex, 
        ymax = mean_null_sdindex + sd_null_sdindex
      ), 
      width = .5,
      color = red
    ) +
    geom_line(color = red) +
    geom_point(color = red) +
    facet_wrap(~ niche_index_scenario, scales = 'free') +
    ggtitle('Mean standardized SD index of null communities') + 
    ylim(c(-3, 3))
  
  plot_gap =
    results_for_plot %>%
    ggplot(aes(num_clusters, gap)) +
    geom_errorbar(aes(x = num_clusters, ymin = gap - sgap, ymax = gap + sgap), width = .5, color = blue) +
    geom_line(color = blue) +
    geom_point(color = blue) +
    facet_wrap(~ niche_index_scenario, scales = 'free') +
    ggtitle('Gap index of simulated communities: gap = obs - mean(nulls)')
  
  gridExtra::grid.arrange(plot_index, plot_gap, nrow = 1)
  
  plot_null_index = 
    dat %>% 
    filter(null > 0) %>% 
    group_by(niche_index_scenario, census_block, cell_length, num_clusters) %>% 
    mutate(mean_index = mean(total_withinss)) %>% 
    ungroup %>% 
    filter(cell_length == 3, census_block == 0) %>% 
    ggplot(aes(num_clusters, mean_index)) + 
    geom_line(color = red) + 
    geom_point(color = red) + 
    geom_vline(aes(xintercept = 20), color = green) + 
    facet_wrap(~ niche_index_scenario, scales = 'free_y')
  
  plot_kmeans_index = 
    dat %>% 
    mutate(null = factor(null)) %>%
    group_by(niche_index_scenario, cell_length, null, census_block) %>%
    mutate(
      mean_withinss = total_withinss / num_clusters,
      scaled_index = (total_withinss - mean(total_withinss)) / sd(total_withinss)
    ) %>%
    ungroup %>% 
    filter(cell_length == 3, census_block == 0) %>% 
    ggplot(aes(num_clusters, mean_withinss, group = null, color = null)) + 
    geom_line() + 
    geom_point() + 
    geom_vline(aes(xintercept = 20), color = green) + 
    facet_wrap(~ niche_index_scenario, scales = 'free_y')
  
  best_numclust = 
    results %>%
    filter(cell_length == chosen_cell_length) %>%
    filter(census_block == chosen_census_block) %>%
    filter(num_clusters >= min_numclust) %>% 
    group_by(niche_index_scenario) %>% 
    slice(which.max(gap)) %>% 
    select(cell_length, niche_index_scenario, num_clusters_by_gap = num_clusters) %>%
    left_join(
      results %>%
        filter(cell_length == chosen_cell_length) %>%
        filter(census_block == chosen_census_block) %>%
        filter(num_clusters >= min_numclust) %>% 
        group_by(niche_index_scenario) %>% 
        slice(which.max(standardized_sdindex)) %>% 
        select(cell_length, niche_index_scenario, num_clusters_by_sdindex = num_clusters)
      
    )
  
  plot_best_numclust = 
    best_numclust %>%
    ggplot() +
    geom_line(aes(niche_index_scenario, num_clusters_by_gap), color = blue) +
    geom_point(aes(niche_index_scenario, num_clusters_by_gap), color = blue) +
    geom_line(aes(niche_index_scenario, num_clusters_by_sdindex)) +
    geom_point(aes(niche_index_scenario, num_clusters_by_sdindex)) +
    geom_hline(aes(yintercept = 20), color = green) +
    labs(
      x = 'Niche index scenario',
      y = 'Best number of clusters'
    )
  ggtitle('Best number of clusters.         black = SD index,  blue = gap index.         Range of search = 3 to 30 clusters')
  
  gridExtra::grid.arrange(plot_best_numclust)
  
  result = list()
  for(i in seq(nrow(best_numclust))){
    
    print(i)
    
    scen = best_numclust[i, ]
    scenario = scen$niche_index_scenario
    cell_length = as.numeric(scen$cell_length)
    num_clusters = scen$num_clusters_by_gap
    
    
    dat = get(load(paste0('scenario_', scenario, '.RData')))
    community = matrix(unlist(dat$evolving_community), length(dat$Landscape))
    
    grid = Square_Grid(landscape_area = nrow(community), cell_length = cell_length)
    
    threshold_census = 17
    selected_columns = threshold_census:ncol(community)
    
    averaged_community = 
      cbind(
        grid = grid,
        community[, selected_columns]
      )
    colnames(averaged_community) = c('grid', paste0('census_', formatC(selected_columns, digits = 1, flag = '0')))
    
    averaged_community = as_tibble(averaged_community)
    
    avg_community = 
      averaged_community %>%
      pivot_longer(cols = -grid, names_to = 'census', values_to = 'species')
    
    celltable = with(avg_community, table(grid, species)) / length(selected_columns)
    
    res =
      Infer_Affinity_Matrix_from_CellTable(
        celltable = celltable,
        num_clusters = num_clusters,
        index = 'sdindex'
      )
    
    result = c(result, res)
    
  }
  
  
  
  ## ========== Set censuses to use in analysis ==========
  # threshold_census = ifelse(scenario < 15, 17, 3)
  # if(scenario == 20) threshold_census = ncol(community)
  threshold_census = 17
  column_index = threshold_census:ncol(community)
  
  
}

if(analyze_inference_on_averaged_cells_save_all_NST_results){
  theresults = NULL
  for(scenario in unique(trimmed_scenarios$scenario)){
    for(cell_length in unique(trimmed_scenarios$cell_length)){
      dat = get(load(paste0('analysis_allcensuses_scenario_', scenario, '_celllength_', cell_length, '.RData')))
      dtf = 
        as_tibble(t(sapply(dat, function(l) l$Best.nc))) %>%
        bind_cols(
          sapply(
            dat,
            function(l){
              C = l$true_C
              eC = 
                l$inferred_matrix %>%
                pivot_wider(
                  names_from = 'soil_type', 
                  values_from = 'affinity'
                ) %>% 
                select(-species) %>%
                as.matrix()
              true_ni = 1 - mean_cosine_rows(C)
              estimated_ni = 1 - mean_cosine_rows(eC)
              RV_index = RV_coefficient(C, eC)$RV
              return(
                c(
                  true_niche_index = true_ni,
                  estimated_niche_index = estimated_ni,
                  RV = RV_index
                )
              )
            }
          ) %>%
            t() %>%
            as_tibble
        ) %>%
        mutate(
          scenario = scenario,
          cell_length = cell_length
        )
      
      theresults = 
        theresults %>%
        rbind(dtf)
    }
  }
  
  res = theresults
  
  plot1 = 
    res %>%
    mutate(cell_length = factor(cell_length)) %>%
    ggplot(aes(Number_clusters, Value_Index, group = cell_length, color = cell_length)) +
    geom_line() +
    ggtitle(paste0('true niche index = ', true_niche_index)) +
    facet_wrap(~scenario, scales = 'free', ncol = 10)
}
