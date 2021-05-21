## Trait and Soil by Cluster Machine Learning

## Determine whether working on SeaWulf (SBU hpc) or personal computer
seawulf = as.logical(Sys.info()['user'] == 'rdrocha')

## Read array argument from command line
index = as.numeric(commandArgs(TRUE))[1]
if(is.na(index)) index = 1

## operation specifiers
trait_based_analysis = 0
soil_based_analysis = 1
do.naive.bayes = 0
do.replicates = 0 ## repeatedly perform the C5.0 classifier on the raw trait data
do.pca = 1 ## perform the C5.0 classifier on pca data


## libraries
packages = 
  c(
    'e1071', 
    'tidyverse', 
    'readxl', 
    'magrittr', 
    'gmodels', 
    'caret', 
    'irr', 
    'C50', 
    'fgpt', 
    'pcaMethods',
    'future.apply'
  )
install.packages(setdiff(packages, rownames(installed.packages())))  
sapply(packages, library, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)

plan(multisession)

## plotting theme
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

## parameters
parms = expand_grid(scale = c(0, 30, 50, 100, 200, 500, 10000), run = 1:100)
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

groups = get(load('~/SpatialNiche/Data/species_by_soiltype.rdata'))

## Soil-based analysis

if(soil_based_analysis){
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
  
  nutrients = read_excel('~/SpatialNiche/Data/bci.block20.data-original.xls', sheet = 2)
  
  nutrients %<>%
    pivot_longer(-c(x, y), names_to = 'nutrient') %>% 
    group_by(nutrient) %>%
    mutate(standardized = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup
  
  unified = 
    soiltype %>%
    mutate(mark = as.factor(soiltype)) %>%
    select(x, y, mark, value = standardized) %>%
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
  
  ## compare C5 on nutrients vs coordinates
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
      
      parms = expand_grid(scale = c(0, 30, 50, 100, 200, 500, 10000), run = 1:100)
      
      filenames = 
        intersect(
          parms %>%
            mutate(filename = paste0('nutrient_null_scale_', scale, '_run_', run,'.rdata')) %>%
            pull(filename),
          list.files()
        )
      
      res = NULL
      for(name in filenames){
        foo = try(get(load(name)), silent = TRUE)
        if(class(foo) != 'try-error' & 'scenario' %in% names(foo)){
          bar = bind_cols(foo$scenario, Kappa = foo$model$results %>% slice_max(Kappa) %>% pull(Kappa))
          res = rbind(res, bar)
        }
      }
      
      save(res, file = '~/SpatialNiche/Data/partial_data.rdata')
      
    }
    
    if(!seawulf){
      dat = get(load('~/SpatialNiche/Data/partial_data.rdata'))
      
      dat %>% 
        # mutate(scale = factor(scale)) %>% 
        ggplot(aes(x = scale, y = Kappa, group = scale, fill = scale)) + 
        geom_boxplot() + 
        theme(legend.position = 'none') + xlab('scale (meters)') + 
        ggtitle('C5.0 classifier - Cohen\'s Kappa')
    }
    
  }
  
  ## plot null data by scales for visual interpretation
  if(TRUE){
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
      scale_fill_continuous(low = 'grey95', high = 'black') +
      facet_wrap(~scale)
    
    plot %>% show
  }
  
  
}



## Trait-based analysis

if(trait_based_analysis){
  trait_data = read_excel('~/BCI/BCITRAITS_20101220.xlsx')
  
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
  
  
  dat = 
    trait_data %>%
    mutate(sp = tolower(`SP$`)) %>%
    pivot_longer(-c(1:6, sp), names_to = 'trait') %>%
    right_join(groups) %>%
    filter(!is.na(value)) %>%
    filter(!str_detect(trait, '_N')) %>%
    filter(!str_detect(trait, 'N_')) %>%
    filter(!str_detect(trait, '_SE')) %>%
    filter(!str_detect(trait, 'SEM_')) %>%
    filter(value > 0) %>%
    mutate(logged_value = log(value))
  
  normality = 
    dat %>%
    group_by(trait) %>%
    summarize(
      normal = shapiro.test(value)$p.value > .05,
      lognormal = shapiro.test(logged_value)$p.value > .05
    ) %>%
    ungroup
  
  trait_dtf =
    dat %>%
    left_join(normality) %>%
    mutate(standardized = ifelse(normal, value, logged_value)) %>%
    group_by(trait) %>%
    mutate(standardized = scale(standardized)[, 1]) %>%
    ungroup %>%
    select(-c(1, 2, 3, 6, value, logged_value, normal, lognormal)) %>%
    rename(grwfrm1 = `GRWFRM1$`, grwfrm2 = `GRWFRM2$`) %>%
    mutate(grwfrm1 = as.factor(grwfrm1), grwfrm2 = as.factor(grwfrm2)) %>%
    select(sp, group, everything()) %>%
    filter(trait %in% unlist(traitlist))%>%
    pivot_wider(names_from = trait, values_from = standardized) %>%
    select(-c(sp, grwfrm2))
  
  
  if(do.replicates){
    ## Implement 100 iterations of the C5.0 classifier, 
    ## extract Cohen's Kappa from each
    Kappa = 
      sapply(1:100, function(index){
        
        set.seed(index)
        
        model = 
          train(
            group ~ ., 
            data = trait_dtf, 
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
        
        writeLines(paste('run', nullrun, ' Kappa = ', kappa))
        
        return(kappa)
      })
    
    
    ## Calculate a null set of Kappa's from randomized versions of the data
    ## where the group class is shuffled.
    nullKappa = sapply(1:100, function(nullrun){
      
      set.seed(nullrun)
      
      null_dtf = 
        trait_dtf %>%
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
      
      writeLines(paste('run', nullrun, ' Kappa = ', nullk))
      
      return(nullk)
    })
    
    kappa_tbl = 
      tibble(
        observed = Kappa,
        null = nullKappa
      ) %>%
      rowid_to_column() %>%
      pivot_longer(-rowid, names_to = 'type', values_to = 'kappa')
    
    save(kappa_tbl, file = '~/SpatialNiche/Data/c50_cohens_kappa.rdata')
    
    plot_kappa = 
      kappa_tbl %>%
      ggplot(aes(x = type, y = kappa, fill = type)) +
      geom_boxplot() +
      theme(
        axis.title.x = element_blank(),
        legend.position = 'none'
      ) +
      ggtitle('Cohen\'s kappa')
    
    
    plot_kappa %>% show
    
  }
  
  
  if(do.pca){
    pca_tibble = get(load('~/SpatialNiche/Data/trait_pca_tibble.rdata'))
    
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
