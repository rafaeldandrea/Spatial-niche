## Trait and Soil by Cluster Machine Learning

trait_based_analysis = 0
soil_based_analysis = 1
do.naive.bayes = 1


packages = c('e1071', 'tidyverse', 'readxl', 'magrittr', 'gmodels', 'caret', 'irr', 'C50')
install.packages(setdiff(packages, rownames(installed.packages())))  
sapply(packages, library, character.only = TRUE, quietly = TRUE)

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)

if(!exists('communities')) source('~/SpatialNiche/R Scripts/DistanceBasedAnalysisByCensus.R')

groups = 
  tibble(
    sp = communities %>% membership %>% names,
    group = communities %>% membership %>% unclass %>% as.factor
  )


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
  
  C5_model = train(soiltype ~ ., data = dtf_wide[, 6:19], method = 'C5.0')
  C5_predict = predict(C5_model, dtf_wide$soiltype)
  
}



## Trait-based analysis

if(trait_based_analysis){
  trait_data = read_excel('~/BCI/BCITRAITS_20101220.xlsx')
  
  groups = 
    tibble(
      sp = communities %>% membership %>% names,
      group = communities %>% membership %>% unclass %>% as.factor
    )
  
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
  
  dat %<>%
    left_join(normality) %>%
    mutate(standardized = ifelse(normal, value, logged_value)) %>%
    group_by(trait) %>%
    mutate(standardized = scale(standardized)[, 1]) %>%
    ungroup %>%
    select(-c(1, 2, 3, 6, value, logged_value, normal, lognormal)) %>%
    rename(grwfrm1 = `GRWFRM1$`, grwfrm2 = `GRWFRM2$`) %>%
    mutate(grwfrm1 = as.factor(grwfrm1), grwfrm2 = as.factor(grwfrm2)) %>%
    select(sp, group, everything()) %>%
    pivot_wider(names_from = trait, values_from = standardized)
  
  C5_model = train(group ~ ., data = dat %>% select(-sp, grwfrm2, WSG_CHAVE), method = 'C5.0')
  C5_predict = predict(C5_model, dtf_wide$soiltype)
}
