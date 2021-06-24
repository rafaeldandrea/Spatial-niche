## Compare results from envt niche analysis to species traits

library(tidyverse)
library(readxl)
library(pcaMethods)
library(magrittr)
library(nnet)

do.pca = 1

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange'),
  aspect.ratio = 1
)


groups = 
  readRDS('~/SpatialNiche/Data/bci_groups_cluster_louvain.rds') %>%
  rename(group = community)
  
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
  right_join(groups, by = 'sp') %>%
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
  ungroup

species_per_group = 
  dat %>% 
  group_by(trait, group) %>% 
  count %>% 
  ungroup %>%
  pivot_wider(
    names_from = group, 
    values_from = n
  )

trait_means = 
  dat %>%
  group_by(trait, group) %>%
  summarize_at('value', list(mean = mean, sd = sd)) %>%
  ungroup

plot_bars =
  trait_means %>%
  ggplot(aes(group, mean)) +
  geom_col() +
  facet_wrap(~trait, scales = 'free')

plot_violin =
  dat %>%
  ggplot(aes(group, value, group = group, fill = group)) +
  geom_violin() +
  facet_wrap(~trait, scales = 'free')


plot_violin_logged =
  dat %>%
  ggplot(aes(group, logged_value, group = group, fill = group)) +
  geom_violin() +
  facet_wrap(~trait, scales = 'free')


plot_violin_standardized =
  dat %>%
  filter(trait %in% unlist(traitlist)) %>%
  ggplot(aes(group, standardized, group = group, fill = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_violin() +
  facet_wrap(~trait, ncol = 8) +
  theme(legend.position = 'none') +
  ggtitle('Distribution of traits - normalized and standardized values')


## PCA analysis

pca_dat = 
  seq_along(traitlist) %>%
  map_dfr(
    .f = function(i){
      subdat = 
        dat %>%
        filter(trait %in% traitlist[[i]]) %>%
        select(sp, trait, standardized) %>%
        pivot_wider(names_from = trait, values_from = standardized) 
      
      subpca = 
        subdat %>%
        pca(method = 'ppca', scale = 'none', center = FALSE)
      
      tibble(
        sp = subdat$sp,
        trait = names(traitlist[i]),
        pc1 = subpca@scores[, 1],
        pc2 = subpca@scores[, 2]
      ) %>%
        return
    }
  ) %>%
  left_join(
    dat %>%
      select(sp, group) %>%
      unique,
    by = 'sp'
  )

traits_and_types = sapply(traitlist, length)

full_dat = 
  pca_dat %>% 
  rename(type = trait) %>% 
  left_join(
    traitlist %>% 
      unlist %>% 
      as_tibble %>% 
      mutate(
        type = rep(names(traits_and_types), traits_and_types)
      ) %>% 
      rename(trait = value), 
    by = 'type'
  ) %>% 
  inner_join(
    dat %>% 
      select(sp, trait, standardized), 
    by = c('sp', 'trait')
  )

correlations = 
  full_dat %>% 
  group_by(type, trait) %>% 
  summarize(
    cor_pc1 = cor(pc1, standardized), 
    cor_pc2 = cor(pc2, standardized),
    .groups = 'drop'
  ) %>%
  group_by(type) %>% 
  summarize_at(c('cor_pc1', 'cor_pc2'), mean)

full_dat %<>% 
  left_join(correlations, by = 'type') %>%
  mutate(
    pc1 = pc1 * sign(cor_pc1),
    pc2 = pc2 * sign(cor_pc2)
  ) %>%
  select(-c(cor_pc1, cor_pc2))

pca_wide = 
  pca_dat %>%
  select(-pc2) %>%
  pivot_wider(names_from = trait, values_from = pc1)

plot_pca =
  pca_dat %>%
  ggplot(aes(pc1, pc2, color = group)) +
  geom_point() +
  facet_wrap(~trait, scales = 'free')

## Multinomial logistic regression
baseline = relevel(pca_wide$group, ref = '2')
test = multinom(baseline ~ vital + leaf + seed + wood + size, data = pca_wide)

## 2-tailed z test
z = with(summary(test), coefficients / standard.errors)
p = (1 - pnorm(abs(z), 0, 1)) * 2

## MANOVA analysis
summary(manova(as.matrix(pca_wide[, 3:7]) ~ pca_wide$group))

comparisons = 
  apply(
    which(
      upper.tri(diag(groups %>% pull(group) %>% unique %>% length)), 
      arr.ind = TRUE
    ), 
    1, 
    function(x) paste(x, collapse = '')
  )

significant_pvalues = 
  pca_dat %>%
  group_by(trait) %>%
  summarize(
    p.value = as.numeric(pairwise.wilcox.test(pc1, group)$p.value)) %>%
  drop_na(p.value) %>%
  mutate(comparison = comparisons) %>%
  ungroup %>% 
  mutate(significant = p.value < .05)

## use significant_pvalues to create significance labels to be added to violin plot
labels_tbl =
  expand_grid(
    trait = c('leaf', 'seed', 'size', 'vital', 'wood'),
    group = as.factor(1:3)
  ) %>%
  # mutate(
  #   label = 
  #     c(
  #       'ab', 'ab', 'c', 
  #       'abc', 'ab', 'ac', 
  #       'abc', 'abc', 'abc',
  #       'ac', 'b', 'ac',
  #       'ac', 'bc', 'bc'
  #     )
  mutate(
    label = 
      c(
        'ac' ,  'b' ,  'ac', 
        'abc', 'abc', 'abc', 
        'abc', 'abc', 'abc',
        'a'  ,  'bc',  'bc',
        'ab' , 'abc',  'bc'
      )
  )

plot_violin_pc1 = 
  pca_dat %>% 
  ggplot(aes(group, pc1, fill = group)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') + 
  geom_violin() + 
  geom_text(aes(group, 0, label = label), data = labels_tbl) +
  facet_wrap(~trait, scales = 'free') +
  theme(legend.position = 'none')

plot_violin_pc2 = 
  pca_dat %>% 
  ggplot(aes(group, pc2, fill = group)) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_violin() + 
  facet_wrap(~trait, scales = 'free') +
  theme(legend.position = 'none')

plot_pc1_vs_traits = 
  full_dat %>%
  mutate(
    type = 
      factor(
        type, 
        levels = c('vital', 'size', 'leaf', 'seed', 'wood')
      )
  ) %>%
  ggplot(aes(standardized, pc1)) + 
  geom_point() + 
  facet_wrap(type ~ trait, ncol = 8)

plot_violin_standardized =
  full_dat %>%
  mutate(
    type = 
      factor(
        type, 
        levels = c('vital', 'size', 'leaf', 'seed', 'wood')
      )
  ) %>%
  ggplot(aes(group, standardized, group = group, fill = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_violin() +
  facet_wrap(type ~ trait, ncol = 8) +
  theme(legend.position = 'none') +
  ggtitle('Species trait distribution by group - normalized and standardized values')

plot_violin_pc1 %>% show
plot_pc1_vs_traits %>% show
plot_violin_standardized %>% show




