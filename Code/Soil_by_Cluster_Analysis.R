library(tidyverse)
library(readxl)
library(magrittr)
library(fgpt)
library(pcaMethods)

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

plot_pca_colored_by_soiltype = 
  dtf_wide %>% 
  ggplot(aes(pc1, pc2, color = soiltype)) + 
  geom_point() +
  theme(legend.position = c(.9, .1))

plot_pca_colored_by_nutrient = 
  dtf_wide %>% 
  pivot_longer(6:18, names_to = 'nutrient') %>%
  ggplot(aes(pc1, pc2, color = value)) +
  geom_point() +
  facet_wrap(~nutrient) +
  scale_color_continuous(low = 'grey95', high = 'black')

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


