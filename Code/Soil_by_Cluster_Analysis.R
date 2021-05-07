library(tidyverse)
library(readxl)
library(magrittr)
library(fgpt)

theme_set(theme_bw())
theme_update(
  aspect.ratio = .5,
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)

do.fgpt = 1

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

plot_rasters = 
  unified %>%
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_wrap(~mark, scales = 'free') +
  scale_fill_gradient2(low = "#d7191c", mid = "#ffffbf", high = "#1a9641", midpoint = .5)

plot_histograms = 
  unified %>%
  ggplot(aes(value)) +
  geom_histogram(fill = rgb(153 / 255, 0, 0)) +
  facet_wrap(~mark)

plot_rasters

nutrient_soiltype_correlations = 
  nutrients %>% 
  left_join(soiltype, by = c('x', 'y')) %>% 
  group_by(nutrient, soiltype) %>% 
  summarize(
    cor = cor(standardized.x, standardized.y), 
    p.value = cor.test(standardized.x, standardized.y)$p.value, 
    p.adj = p.adjust(p.value, method = 'bonferroni'),
    significant = (p.adj < .05),
    cor_sig = ifelse(significant == TRUE, cor, NA),
    .groups = 'drop'
  ) 


plot_nutrient_soiltype_correlations =
  nutrient_soiltype_correlations %>%
  ggplot(aes(soiltype, nutrient, fill = cor_sig)) +
  geom_tile() +
  scale_fill_gradient2(low = "#d7191c", mid = "#ffffbf", high = "#1a9641", midpoint = 0) +
  theme(aspect.ratio = 2)

nutrient_nutrient_correlations = 
  nutrients %>% 
  left_join(nutrients, by = c('x', 'y')) %>% 
  group_by(nutrient.x, nutrient.y) %>% 
  summarize(
    cor = cor(standardized.x, standardized.y), 
    p.value = cor.test(standardized.x, standardized.y)$p.value, 
    p.adj = p.adjust(p.value, method = 'bonferroni'),
    significant = (p.adj < .05),
    cor_sig = ifelse(significant == TRUE, cor, NA),
    .groups = 'drop'
  ) 
  
plot_nutrient_nutrient_correlations = 
  nutrient_nutrient_correlations %>%
  ggplot(aes(nutrient.x, nutrient.y, fill = cor_sig)) +
  geom_tile() +
  scale_fill_gradient2(low = "#d7191c", mid = "#ffffbf", high = "#1a9641", midpoint = 0) +
  theme(aspect.ratio = 1)

gridExtra::grid.arrange(
  plot_nutrient_nutrient_correlations,
  plot_nutrient_soiltype_correlations,
  nrow = 1
)


## Floating Grid Permutation analysis (see vignette at https://cran.r-project.org/web/packages/fgpt/vignettes/intro-fgpt.pdf)
if(do.fgpt){
  
  unified_wider = 
    unified %>%
    pivot_wider(names_from = mark, values_from = value)
  
  fgcor_3Nmin = 
    with(
      unified_wider, 
      fgeasy(
        xy = cbind(x, y), 
        marks = cbind(`3`, `N(min)`), 
        pairwise = TRUE, 
        correlate = 'pearson', 
        iter = 99
      )
    )
  
  plot(fgcor_3Nmin, main = 'Soil type 3 versus Nitrogen (mineralized)')
  
  fgcor_23 = 
    with(
      unified_wider, 
      fgeasy(
        xy = cbind(x, y), 
        marks = cbind(`2`, `3`), 
        pairwise = TRUE, 
        correlate = 'pearson', 
        iter = 99
      )
    )
  
  plot(fgcor_23, main = 'Soil type 2 versus 3')
  
}


