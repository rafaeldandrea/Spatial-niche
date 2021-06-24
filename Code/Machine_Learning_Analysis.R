## Trait and Soil by Cluster Machine Learning

## Determine whether working on SeaWulf (SBU hpc) or personal computer
seawulf = as.logical(Sys.info()['user'] == 'rdrocha')

## Read array argument from command line
index = as.numeric(commandArgs(TRUE))[1]
if(is.na(index)) index = 1

## operation specifiers
trait_based_analysis = 0
soil_based_analysis = 0
do.naive.bayes = 0
do.replicates = seawulf          ## Repeatedly perform the C5.0 classifier on the raw trait data
do.pca = 0                       ## Perform the C5.0 classifier on pca data
do.partial.randomization = 0     ## Test C5.0 on partially randomized nutrient features (package fgpt)
do.randomfield = 0               ## Test C5.0 on random fields with same autocorrelation scale, number, 
                                 ##      and prevalence of soiltypes as our data
do.matrix.analysis = 1           ## Estimate the affinity matrix

## libraries
packages = 
  c(
    'tidyverse', 
    'readxl', 
    'magrittr', 
    'gmodels', 
    'caret', 
    'irr', 
    'C50', 
    'fgpt', 
    'pcaMethods',
    'gstat',        ## for function variogram()
    'RandomFields', ## for function RandomFields() defined below
    'cowplot',       ## for function plot_grid() 
    #'slurmR'         ## for parallel computing function Slurm_sapply()
    'furrr',
    'rio',           ## for function import() to read xls file directly from github
    'e1071'          ## for C5.0 algorithm
  )

install.packages(setdiff(packages, rownames(installed.packages())))  
sapply(packages, library, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)

plan(multisession, workers = 40 * seawulf + 7 * !seawulf)

## plotting theme
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

## parameters
parms = 
  expand_grid(
    scale = c(0, 30, 50, 100, 200, 500, 10000, 1000, 5000, 50000), 
    run = 1:100
  )
scenario = parms[index, ]

if(!seawulf){
  
  if(!file.exists('~/SpatialNiche/Data/species_by_soiltype.rdata')){
    
    source('~/SpatialNiche/R Scripts/DistanceBasedAnalysisByCensus.R')
    
    groups = 
      tibble(
        sp = communities %>% membership %>% names,
        group = communities %>% membership %>% unclass %>% as.factor
      )
    
    save(groups, file = '~/SpatialNiche/Data/species_by_soiltype.rdata')
  }  
    
}

# groups = 
#   readRDS(
#     url(
#       'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/species_by_soiltype.rds?raw=true'
#     )
#   )

groups = 
  readRDS(
    url(
      'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_groups_cluster_louvain.rds?raw=true'
    )
  ) %>%
  mutate(group = community)

## Soil-based analysis

if(soil_based_analysis){
  # dtf = NULL
  # for(i in 1:3){
  #   dtf = 
  #     dtf %>%
  #     bind_rows(
  #       read.csv(
  #         paste0(
  #           'https://raw.githubusercontent.com/rafaeldandrea/Spatial-niche/main/Data/soiltype',i,'_dist_10_abd_50.csv'
  #         )
  #       ) %>%
  #       as_tibble %>%
  #       mutate(soiltype = i)
  #     )
  # }
  
  # soiltype = 
  #   dtf %>%
  #   pivot_longer(-c(X, soiltype), names_to = 'Y') %>%
  #   mutate(Y = as.numeric(str_remove(Y, 'X'))) %>%
  #   mutate(x = (1 + X) * 20 - 10, y = (1 + Y) * 20 - 10) %>%
  #   select(X, Y, x, y, soiltype, value) %>%
  #   group_by(soiltype) %>%
  #   mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
  #   ungroup
  
  soiltype =
    readRDS(
      url(
        'https://github.com/rafaeldandrea/Spatial-niche/blob/main/Data/bci_20by20grid_soiltype_louvain.rds?raw=true'
      )
    )
  
  nutrients =
    rio::import(
      'https://github.com/rafaeldandrea/Spatial-niche/raw/main/Data/bci.block20.data-original.xls', 
      which = 2
    ) %>% 
    pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
    group_by(nutrient) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup
  
  unified = 
    soiltype %>%
    mutate(mark = as.factor(soiltype)) %>%
    select(x, y, mark, value = prob) %>%
    bind_rows(
      nutrients %>%
        mutate(mark = nutrient) %>%
        select(x, y, mark, value = standardized)
    ) 
  
  dtf_wide = 
    unified %>%
    pivot_wider(names_from = mark, values_from = value) %>%
    group_by(x, y) %>% 
    mutate(soiltype = which.max(c(`1`, `2`, `3`))) %>% 
    ungroup %>%
    mutate(soiltype = factor(soiltype, levels = c(1, 2, 3)))
  
  pca_model =
    dtf_wide %>%
    select(6:18) %>%
    pca(method = 'ppca', scale = 'none', center = FALSE)

  dtf_wide %<>%
    bind_cols(
      pc1 = pca_model@scores[, 1],
      pc2 = pca_model@scores[, 2]
    )
  
  if(do.naive.bayes){
    percentage_train = .8
    set.seed(0)
    train_index = sort(sample(nrow(dtf_wide), size = percentage_train * nrow(dtf_wide)))
    test_index = setdiff(seq(nrow(dtf_wide)), train_index)
    
    train = 
      dtf_wide[train_index, ] %>%
      select(6:18)
    
    train_class = 
      dtf_wide[train_index, ] %>%
      select(soiltype) %$%
      soiltype
    
    test =
      dtf_wide[test_index, ] %>%
      select(6:18)
    
    test_class = 
      dtf_wide[test_index, ] %>%
      select(soiltype) %$%
      soiltype
    
    NB_model = naiveBayes(train, train_class, laplace = 1)
    NB_predict = predict(NB_model, test)
    
    NB_cross_table = 
      CrossTable(
        NB_predict, 
        test_class, 
        prop.chisq = FALSE,
        prop.c = FALSE,
        prop.r = FALSE,
        dnn = c('predicted', 'actual')
      )
    
    NB_confusion_matrix =
      confusionMatrix(
        as.factor(NB_predict),
        as.factor(test_class)
      )
    
    set.seed(0)
    features = dtf_wide %>% select(6:18)
    soiltypes = dtf_wide %$% soiltype
    
    repeated_kfold_crossvalidation =
      sapply(1:10, function(i){
        folds = createFolds(soiltypes, k = 10)
        cv_results = 
          lapply(folds, function(x){
            train = features[-x, ]
            test = features[x, ]
            nb_model = naiveBayes(train, soiltypes[-x], laplace = 1)
            nb_predict = predict(nb_model, test)
            kappa = kappa2(tibble(soiltypes[x], nb_predict))$value
            return(kappa)
          }) %>%
          unlist %>%
          mean
      })
    
    
    mean_kappa = mean(repeated_kfold_crossvalidation)
    boxplot_cv_results = 
      repeated_kfold_crossvalidation %>%
      as_tibble %>%
      mutate(dummy = 'a', kappa = value) %>%
      ggplot(aes(dummy, kappa)) +
      geom_boxplot(fill = 'orange') +
      coord_cartesian(ylim = c(0, 1)) +
      theme(aspect.ratio = 5) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
      ggtitle(paste('Naive Bayes mean kappa =', round(mean_kappa, 2)))
    
    boxplot_cv_results %>%
      show
  }
  
  soil_dtf = 
    dtf_wide %>%
    select(6:19)
  
  ## simple test - shorter runtime
  if(FALSE){
    C5_model = train(soiltype ~ ., data = soil_dtf, method = 'C5.0')
    C5_predict = predict(C5_model, soil_dtf)
  }
  
  ## compare C5.0 on nutrients vs coordinates
  if(FALSE){
    C5_model_soil =
      train(
        soiltype ~ .,
        data = soil_dtf,
        method = 'C5.0',
        na.action = na.pass,
        trControl = trainControl(method = 'repeatedcv', repeats = 10),
        metric = 'Kappa'
     )

    C5_model_coords =
      train(
        soiltype ~ .,
        data = dtf_wide %>% select(x, y, soiltype),
        method = 'C5.0',
        na.action = na.pass,
        trControl = trainControl(method = 'repeatedcv', repeats = 10),
        metric = 'Kappa'
      )
  }
  
  ## test C5.0 on Gaussian random field
  if(do.randomfield){
    seed = scenario$run
    
    RandomField = function(
      Lx = 50, 
      Ly = Lx, 
      rangepar, 
      sillpar, 
      nuggetpar, 
      seed = seed, 
      plot = FALSE
    ) {
        stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
        
        RFoptions(seed = seed)
        
        stress = 
          RFsimulate(
            RMgauss(
              scale = rangepar + 1e-16,
              var = 2 * sillpar * (1 - nuggetpar)
            ) + RMtrend(mean = 0) + RMnugget(var = 2 * sillpar * nuggetpar),
            x = 1:Lx,
            y = 1:Ly
          )@data$variable1
        
        if (plot)
          plot(
            raster::raster(matrix(stress, Lx, Ly)),
            las = 1,
            xlab = 'x-coordinate',
            ylab = 'y-coordinate'
          )
        
        return(stress)
        
    }
    
    Gaussian_vgm_optim = function(parms, sample.vgm){
      range = parms[1]
      sill = parms[2]
      nugget = parms[3]
      dist = c(0, sample.vgm$dist)
      observed = c(0, sample.vgm$gamma)
      predicted = (sill - nugget) * (1 - exp(-dist ^ 2 / range ^ 2)) + nugget * (dist > 0)
      sum_squared_errors = sum((observed - predicted) ^ 2)
      return(sum_squared_errors)
    } 
    
    sample.vgm = 
      variogram(
        soiltype ~ 1, 
        data = dtf_wide, 
        locations = ~ x + y, 
        width = 1
      )
    
    fitted.vgm = 
      optim(
        par = c(range = 20, sill = .6, nugget = 20), 
        fn = Gaussian_vgm_optim,
        sample.vgm = sample.vgm,
        method = "L-BFGS-B",
        lower = c(1, 1, 0)
      )$par
    
    sill = as.numeric(fitted.vgm[2])		## semivariance at the landscape level --> regional heterogeneity
    nugget = as.numeric(fitted.vgm[3])	## semivariance at distance = 1	--> local uniformity
    range = as.numeric(fitted.vgm[1])		## distance at which the semivariance reaches 63% of the sill
    range95 = sqrt(3) * range   ## distance at which the semivariance reaches 95% of the sill
    
    randomfield_test = function(
      seed, 
      .data, 
      .range95, 
      .sill, 
      .func, 
      save.result = seawulf
    ){
      set.seed(seed)
      
      test_data = 
        expand_grid(
          y = seq(10, 490, by = 20),
          x = seq(10, 990, by = 20)
        ) %>%
        mutate(
          field = 
            .func(
              Lx = 1000 / 20, 
              Ly = 500 / 20, 
              rangepar = range95, 
              sillpar = sill, 
              nuggetpar = .015, 
              seed = seed
            ),
          soiltype = 
            cut(
              field, 
              breaks = quantile(
                field, 
                p = c(
                  0, 
                  dtf_wide %>% 
                    count(soiltype) %>% 
                    mutate(p = cumsum(n) / sum(n)) %>% 
                    pull(p)
                )
              ), ## these quantiles match the relative prevalence 
              ## of each soiltype in the real data
              labels = FALSE,
              include.lowest = TRUE
            ) %>%
            as.factor
        )
      
      C5_model =
        train(
          soiltype ~ .,
          data = .data %>% mutate(soiltype = test_data$soiltype),
          method = 'C5.0',
          trControl = trainControl(method = 'repeatedcv', repeats = 10),
          metric = 'Kappa'
        )
      
      if(save.result){
        res = 
          bind_cols(
            seed = seed,
            C5_model$results
          ) %>%
          as_tibble %>%
          slice_max(Kappa, n = 1, with_ties = FALSE)
        
        saveRDS(
          res, 
          file = paste0('~/SpatialNiche/Data/C50_nutrients_RandomField_20210622/randomfield_null_', seed, '.rds')
        )
        
        return(res)
        
      }
    }
    
      # if(seawulf){
        # Slurm_sapply(
        #   X = as.list(1:100), 
        #   FUN = randomfield_test, 
        #   .data = soil_dtf,
        #   .sill = sill,
        #   .range95 = range95,
        #   .func = RandomField,
        #   njobs = 10L,
        #   sbatch_opt = list(partition = 'short-28core')
        # )
      # } 
      
    if(seawulf) ncores = 40 else ncores = 7
    plan(multisession, workers = ncores)
    
    if(seawulf){
      number_of_nulls = 100
      
      results = 
        1:number_of_nulls %>%
        future_map_dfr(
          .f = randomfield_test,
          .data = soil_dtf,
          .sill = sill,
          .range95 = range95,
          .func = RandomField,
          .options = furrr_options(seed = NULL),
          save.result = TRUE
        )
      
      saveRDS(results, file = '~/SpatialNiche/Data/C50_nutrients_RandomField_20210622/randomfield_collated.rds')
      
    }
    
    ## plots
    if(!seawulf){
      
      test_data = 
        expand_grid(
          y = seq(10, 490, by = 20),
          x = seq(10, 990, by = 20)
        ) %>%
        mutate(
          field = 
            RandomField(
              Lx = 1000 / 20, 
              Ly = 500 / 20, 
              rangepar = range95, 
              sillpar = sill, 
              nuggetpar = 0, 
              seed = 10
            ),
          soiltype = 
            cut(
              field, 
              breaks = quantile(
                field, 
                p = c(
                  0, 
                  dtf_wide %>% 
                    count(soiltype) %>% 
                    mutate(p = cumsum(n) / sum(n)) %>% 
                    pull(p)
                  )
                ), ## these quantiles match the relative prevalence 
                   ## of each soiltype in the real data
              labels = FALSE,
              include.lowest = TRUE
            ) %>%
            as.factor
        )
      
      # data = 
      #   1:100 %>%
      #   future_map_dfr(
      #     function(char){ 
      #       foo = 
      #         try(
      #           get(load(paste0('~/SpatialNiche/Data/C50_nutrients/nutrient_null_scale_0_run_',char,'.rdata'))), 
      #           silent = TRUE
      #         )
      #       
      #       if(class(foo) != 'try-error') 
      #         tibble(run = char) %>%
      #         bind_cols(
      #           foo$model$results %>% 
      #             slice_max(Kappa, with_ties = FALSE)
      #         )
      #     }
      #   ) %>%
      #   mutate(type = 'data') %>%
      #   bind_rows(
      #     readRDS('~/SpatialNiche/Data/C50_nutrients_RandomField/randomfield_collated.rds') %>%
      #       mutate(type = 'null') 
      #   )
      
      data = 
        readRDS('~/SpatialNiche/Data/C50_nutrients_RandomField_20210622/randomfield_collated.rds')
      
      plot_test_data = 
        test_data %>% 
        ggplot(aes(x, y, fill = soiltype)) + 
        geom_raster() +
        theme(aspect.ratio = .5) +
        ggtitle('Null')
      
      plot_real_data = 
        dtf_wide %>%
        ggplot(aes(x, y, fill = soiltype)) + 
        geom_raster() +
        theme(aspect.ratio = .5) +
        ggtitle('Data')
      
      plot_kappa =
        data %>%
        filter(seed > 0) %>%
        ggplot(aes(Kappa)) +
        geom_density(fill = 'grey') +
        geom_vline(
          aes(
            xintercept = 
              data %>% 
              filter(is.na(seed)) %>% 
              pull(Kappa)
          ),
          color = 'red',
          size = 2
        ) +
        coord_cartesian(xlim = c(0, 1)) +
        ggtitle('C5.0 Cohen\'s Kappa')
      
      plot_grid(plot_real_data, plot_test_data, plot_kappa, ncol = 1)
      
    } 
  }
  
  ## test C5.0 on partially randomized nutrient features (package fgpt)
  if(do.partial.randomization){
    nutrient_null_c5 = function(scale = 100, run){
      set.seed(run)
      
      if(scale > 0){
        data =
          apply(dtf_wide[, 6:18], 2, function(vector)
            unlist(fgperm(xy = cbind(dtf_wide$x, dtf_wide$y), z = vector, scale = scale, iter = 1))) %>%
          as_tibble %>%
          mutate(soiltype = dtf_wide$soiltype)
      } else{
        data = dtf_wide[6:19]
      }
      
      C5_model_null =
        train(
          soiltype ~ .,
          data = data,
          method = 'C5.0',
          na.action = na.pass,
          trControl = trainControl(method = 'repeatedcv', repeats = 10),
          metric = 'Kappa'
        )
      
      return(C5_model_null)
      
    }
    
    # future_sapply(1:10, nutrient_null_c5, future.seed = NULL)
    
    if(seawulf) model = nutrient_null_c5(scenario$scale, scenario$run)
    
    if(seawulf){
      path = paste0('~/SpatialNiche/Data/C50_nutrients/nutrient_null_scale_', scenario$scale, '_run_', scenario$run, '.rdata')
      data = list(scenario = scenario, model = model)
      save(data, file = path)
    }
    
    collate = FALSE
    if(collate){
      
      if(seawulf){
        library(tidyverse)
        
        parms = expand_grid(scale = c(0, 30, 50, 100, 200, 500, 10000, 1000, 5000, 50000), run = 1:100)
        
        filenames = 
          intersect(
            parms %>%
              mutate(filename = paste0('nutrient_null_scale_', scale, '_run_', run,'.rdata')) %>%
              pull(filename),
            list.files(path = '~/SpatialNiche/Data/C50_nutrients/')
          )
        
        filenames = paste0('~/SpatialNiche/Data/C50_nutrients/', filenames)
        
        res = NULL
        for(name in filenames){
          foo = try(get(load(name)), silent = TRUE)
          if(class(foo) != 'try-error' & 'scenario' %in% names(foo)){
            bar = bind_cols(foo$scenario, Kappa = foo$model$results %>% slice_max(Kappa) %>% pull(Kappa))
            res = rbind(res, bar)
          }
        }
        
        res %<>% unique
        
        save(res, file = '~/SpatialNiche/Data/partial_data.rdata')
        
      }
      
      if(!seawulf){
        dat = get(load('~/SpatialNiche/Data/partial_data.rdata'))
        
        plot0 = 
          dat %>% 
          filter(scale == 0) %>%
          mutate(scale = factor(scale)) %>%
          ggplot(aes(x = scale, y = Kappa)) + 
          geom_boxplot(fill = 'grey') + 
          theme(legend.position = 'none') + 
          xlab('scale (meters)') +
          coord_cartesian(ylim = c(0, 1))
        
        plot1 = 
          dat %>% 
          filter(scale != 0) %>%
          ggplot(aes(x = scale, y = Kappa, group = scale, fill = factor(scale))) + 
          geom_boxplot() + 
          theme(legend.position = 'none') + 
          xlab('scale (meters)') +
          scale_x_log10(breaks = dat %>% pull(scale) %>% unique) +
          coord_cartesian(ylim = c(0, 1))
        
        plot_grid(plot0, plot1, nrow = 1, labels = 'AUTO', rel_widths = c(1, 5))
        
      }
      
    }
    
    ## plot null data by scales for visual interpretation
    if(FALSE){
      scales = parms %>% filter(scale > 0) %>% pull(scale) %>% unique
      data =
        sapply(scales, function(scale){
          unlist(fgperm(xy = base::cbind(dtf_wide$x, dtf_wide$y), z = dtf_wide$N, scale = scale, iter = 1))
        }) %>%
        as_tibble %>%
        bind_cols(dtf_wide %>% select(x, y)) %>%
        pivot_longer(-c(x, y), names_to = 'scale', values_to = 'N') %>%
        mutate(scale = scales[match(scale, paste0('V', 1:6))]) %>%
        bind_rows(
          dtf_wide %>%
            select(x, y, N) %>%
            mutate(scale = 0)
        )
      
      plot = 
        data %>%
        ggplot(aes(x, y, fill = N)) +
        geom_raster() +
        theme(aspect.ratio = .5) +
        scale_fill_gradientn(colors = terrain.colors(7)) +
        facet_wrap(~scale)
      
      plot %>% show
    }
    
  }
  
  
}


## Trait-based analysis
if(trait_based_analysis){
  data = readRDS('~/SpatialNiche/Data/bci_trait_pca_group_louvain.rds')
  
  dat = 
    data %>% 
    select(sp, type, pc1, group) %>% 
    unique %>% 
    pivot_wider(
      names_from = type, 
      values_from = pc1
    ) %>% 
    select(-sp)
  
  if(seawulf){
    
    ## Implement 100 iterations of the C5.0 classifier, 
    ## extract Cohen's Kappa from each
    cohens_kappa = function(index){
      set.seed(index)
      
      model = 
        train(
          group ~ ., 
          data = dat, 
          method = 'C5.0', 
          na.action = na.pass, 
          trControl = 
            trainControl(
              method = 'repeatedcv', 
              repeats = 10,
              selectionFunction = 'oneSE'
            ),
          metric = 'Kappa'
        )
      
      kappa = 
        model$results %>%
        slice_max(Kappa) %>%
        pull(Kappa)
      
      if(seawulf) writeLines(paste('run', index, ' Kappa = ', kappa))
      
      tibble(
        data = 'observed', 
        index = index, 
        kappa = kappa
      ) %>%
        return
    }
    
    ## Calculate a null set of Kappa's from randomized versions of the data
    ## where the group class is shuffled.
    null_kappa = function(nullrun){
      
      set.seed(nullrun)
      
      null_dtf = 
        dat %>%
        mutate(group = sample(group))
      
      model = 
        train(
          group ~ ., 
          data = null_dtf, 
          method = 'C5.0', 
          na.action = na.pass, 
          trControl = 
            trainControl(
              method = 'repeatedcv', 
              repeats = 10,
              selectionFunction = 'oneSE'
            ),
          metric = 'Kappa'
        )
      
      nullk = 
        model$results %>%
        slice_max(Kappa) %>%
        pull(Kappa)
      
      if(!seawulf) writeLines(paste('run', nullrun, ' Kappa = ', nullk))
      
      tibble(
        data = 'null', 
        index = nullrun, 
        kappa = nullk
      ) %>%
        return
    }
    
    Kappa = 
      1:100 %>%
      future_map_dfr(
        .f = cohens_kappa,
        .options = furrr_options(seed = NULL)
      )
    
    nullKappa = 
      1:1000 %>%
      future_map_dfr(
        .f = null_kappa,
        .options = furrr_options(seed = NULL)
      )
    
    kappa_tbl =
      Kappa %>% 
      bind_rows(nullKappa)
    
    saveRDS(kappa_tbl, '~/SpatialNiche/Data/bci_c50_cohens_kappa_louvain.rds')
    
  }
  
  if(!seawulf){
    
    kappa_tbl = readRDS('~/SpatialNiche/Data/bci_c50_cohens_kappa_louvain.rds')
    
    plot_kappa =
      kappa_tbl %>%
      ggplot(aes(data, kappa, fill = data)) +
      geom_violin() +
      theme(
        axis.title.x = element_blank(),
        legend.position = 'none'
      ) +
      ggtitle('Cohen\'s kappa')
    
    
    plot_kappa %>% show
  }
  
  if(do.pca){
    # pca_tibble = get(load('~/SpatialNiche/Data/trait_pca_tibble.rdata'))
    pca_tibble = readRDS('~/SpatialNiche/Data/bci_trait_pca_tibble_louvain.rds')
    
    dtf = 
      pca_tibble %>%
      select(-pc2) %>%
      pivot_wider(names_from = trait, values_from = pc1) %>%
      select(-sp)
    
    if(do.replicates){
      kappa = sapply(1:100, function(index){
        
        set.seed(index)
        
        model = 
          train(
            group ~ ., 
            data = dtf, 
            method = 'C5.0', 
            na.action = na.pass, 
            trControl = 
              trainControl(
                method = 'repeatedcv', 
                repeats = 10,
                selectionFunction = 'oneSE'
              ),
            metric = 'Kappa'
          )
        
        k = 
          model$results %>%
          slice_max(Kappa) %>%
          pull(Kappa)
        
        writeLines(paste('run = ', index, ' Kappa = ', k))
        
        return(k)
        
      })
      
      nullkappa = sapply(1:100, function(index){
        
        set.seed(index)
        
        null_dtf = 
          dtf %>%
          mutate(group = sample(group))
        
        model = 
          train(
            group ~ ., 
            data = null_dtf, 
            method = 'C5.0', 
            na.action = na.pass, 
            trControl = 
              trainControl(
                method = 'repeatedcv', 
                repeats = 10,
                selectionFunction = 'oneSE'
              ),
            metric = 'Kappa'
          )
        
        k = 
          model$results %>%
          slice_max(Kappa) %>%
          pull(Kappa) %>%
          unique
        
        writeLines(paste('run = ', index, ' Kappa = ', k))
        
        return(k)
        
      })
      
      kappa_tbl_pca = 
        tibble(
          observed = kappa,
          null = nullkappa
        ) %>%
        rowid_to_column() %>%
        pivot_longer(-rowid, names_to = 'type', values_to = 'kappa')
      
      save(kappa_tbl_pca, file = '~/SpatialNiche/Data/c50_cohens_kappa_pca.rdata')
      
      plot_kappa_pca = 
        kappa_tbl_pca %>%
        ggplot(aes(x = type, y = kappa, fill = type)) +
        geom_boxplot() +
        theme(
          axis.title.x = element_blank(),
          legend.position = 'none'
        ) +
        ggtitle('Cohen\'s kappa')
      
      
      plot_kappa_pca %>% show
      
    }
    
  }
  
}


## Matrix analysis
if(do.matrix.analysis){
  
  dtf = NULL
  for(i in 1:3){
    dtf = 
      dtf %>%
      bind_rows(
        read.csv(paste0('~/SpatialNiche/Data/soiltype', i, '_dist_10_abd_50.csv')) %>%
          as_tibble %>%
          mutate(soiltype = i)
      )
  }
  
  soiltype = 
    dtf %>%
    pivot_longer(-c(X, soiltype), names_to = 'Y') %>%
    mutate(Y = as.numeric(str_remove(Y, 'X'))) %>%
    mutate(x = (1 + X) * 20 - 10, y = (1 + Y) * 20 - 10) %>%
    select(X, Y, x, y, soiltype, value) %>%
    group_by(soiltype) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup
  
  
  bci = 
    get(load('~/BCI/bci.full7.rdata')) %>% 
    as_tibble %>% 
    filter(dbh >= 100) %>% 
    select(sp, x = gx, y = gy) %>% 
    mutate(
      gx = cut(x, seq(0, 1000, by = 20), labels = FALSE, include.lowest = TRUE), 
      gy = cut(y, seq(0, 1000, by = 20), labels = FALSE, include.lowest = TRUE),
      x = seq(10, 990, by = 20)[gx],
      y = seq(10, 990, by = 20)[gy]
    ) %>%
    inner_join(groups, by = 'sp') %>%
    select(sp, x, y, sniche = group)
  
  # stypes = 
  #   soiltype %>% 
  #   select(x, y, soiltype, standardized) %>% 
  #   pivot_wider(names_from = soiltype, values_from = standardized) %>% 
  #   group_by(x,y) %>% 
  #   summarize(stype = factor(which.max(c(`1`, `2`, `3`)), levels = 1:3), .groups = 'drop')
  
  stypes = 
    readRDS('~/SpatialNiche/Data/bci_20by20grid_soiltype_louvain.rds') %>%
    group_by(x, y) %>%
    slice_max(prob) %>%
    rename(stype = soiltype) %>%
    select(x, y, stype)
  
  abuns = 
    bci %>% 
    count(sp)
  
  matdat = 
    bci %>%
    left_join(stypes, by = c('x', 'y')) %>% 
    group_by(sp, sniche, stype) %>%
    summarize(nij = n(), .groups = 'drop') %>%
    left_join(abuns, by = 'sp') %>%
    mutate(Cij = nij / n) %>%
    select(sp, sniche, stype, Cij) %>%
    pivot_wider(names_from = stype, values_from = Cij) %>%
    replace_na(list(`1` = 0, `2` = 0, `3` = 0))
  
  matdat %<>% 
    filter(sniche == 1) %>% 
    arrange(`1`) %>% 
    bind_rows(
      matdat %>% 
        filter(sniche == 2) %>% 
        arrange(`2`)
    ) %>% 
    bind_rows(
      matdat %>% 
        filter(sniche == 3) %>% 
        arrange(`3`)
    )
  
  Cij_long = 
    matdat %>% 
    pivot_longer(-c(sp, sniche), names_to = 'stype') %>% 
    mutate(sp = factor(sp, levels = matdat$sp)) %>% 
    mutate(match = (sniche == stype)) 
  
  ## This works out to be Cd = 0.5 and Co = 0.25, so that Cd / Co = 2 which is much less than 
  ## the required amount of niche differentiation to see community-level non-neutral behavior 
  ## in D'Andrea, Gibbs and O'Dwyer 2020, Fig. 2E (who used 50 species and mean abundance 500, 
  ## compared to 77 species and mean abundance = 232 for BCI species >= 10cm dbh).
  Cd_Co =
    Cij_long %>% 
    group_by(match) %>% 
    summarize(C = mean(value), .groups = 'drop')
  
  printcolors = 
    Cij_long %>%
    select(sp, sniche) %>%
    unique %>%
    pull(sniche)
  
  Cij_plot = 
    Cij_long %>% 
    ggplot(aes(stype, sp, fill = value)) + 
    geom_raster() +
    scale_fill_gradient2(
      low = 'black', 
      mid = 'white', 
      high = 'red', 
      midpoint = mean(Cij)
    ) +
    theme(
      axis.text.y = 
        element_text(
          color = c('red', 'darkgreen', 'blue')[printcolors]
        ),
      aspect.ratio = 6
    ) +
    labs(
      x = 'inferred soil type', 
      y = 'species', 
      fill = expression(C[ij])
    ) +
    ggtitle('Resource use matrix')
  
  Cij_plot %>% show
          
  
  Cij = as.matrix(matdat[, 3:5])
  
  Aij = Cij %*% t(Cij) / rowSums(Cij ^ 2)
  
  Aij_long = 
    Aij %>%
    as_tibble %>%
    rowid_to_column(var = 'species.x') %>%
    mutate(species.x = paste0('V', species.x)) %>%
    pivot_longer(-species.x, names_to = 'species.y', values_to = 'Aij') %>%
    mutate(
      species.x = str_remove(species.x, 'V'),
      species.y = str_remove(species.y, 'V')
    ) %>%
    mutate(
      species.x = factor(species.x, levels = 1:nrow(matdat)),
      species.y = factor(species.y, levels = 1:nrow(matdat))
    )
  
  Aij_plot = 
    Aij_long %>%
    ggplot(aes(species.x, species.y, fill = Aij)) + 
    geom_raster() +
    scale_fill_gradient2(
      low = 'black', 
      mid = 'white', 
      high = 'red', 
      midpoint = mean(Aij)
    )
}
