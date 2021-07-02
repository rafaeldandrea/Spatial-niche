lap_census = 
  read.table(
    '~/SpatialNiche/Data/relaplanadacensusdata/Censo 2_La Planada.txt', 
    header = TRUE
  ) %>%
  as_tibble %>%
  mutate(sp = tolower(sp)) %>%
  select(sp, gx, gy, dbh) %>%
  mutate(fdp = 'lap')

lap_groups = 
  readRDS('~/SpatialNiche/Data/laplanada_species_by_group_louvain.rds') |>
  mutate(fdp = 'lap')

bci_groups = 
  readRDS('~/SpatialNiche/Data/bci_species_by_group_louvain.rds')

bci_census = 
  get(
    load(
      url(
        'https://github.com/rafaeldandrea/BCI/blob/master/bci.full7.rdata?raw=true'
      )
    )
  ) %>%
  select(sp, gx, gy, dbh) %>%
  unique() %>%
  as_tibble %>%
  mutate(fdp = 'bci')

dat = 
  bci_census |>
  full_join(bci_groups, by = c('sp', 'fdp')) |>
  bind_rows(
    lap_census |>
      full_join(lap_groups, by = c('sp', 'fdp'))
  ) |>
  drop_na(gx, gy, dbh) |>
  filter(dbh > 0) |>
  mutate(
    x = cut(gx, breaks = seq(0, 1000, by = 20), labels = FALSE),
    y = cut(gy, breaks = seq(0, 1000, by = 20), labels = FALSE)
  ) |> 
  drop_na(x, y) |>
  mutate(
    x = 10 * (2 * x - 1),
    y = 10 * (2 * y - 1)
  )

common_sp = 
  dat |> 
  filter(fdp == 'bci') |> 
  pull(sp) |> 
  intersect(
    dat |> 
      filter(fdp == 'lap') |> 
      pull(sp) 
    )

lap_soiltype = 
  readRDS('~/SpatialNiche/Data/laplanada_20by20grid_soiltype_louvain.rds') |>
  mutate(fdp = 'lap') |>
  group_by(x, y) |>
  slice_max(prob) |>
  ungroup() |>
  select(-c(prob, value))

bci_soiltype = readRDS('~/SpatialNiche/Data/bci_20by20grid_soiltype_louvain.rds') |>
  mutate(fdp = 'bci') |>
  group_by(x, y) |>
  slice_max(prob) |>
  ungroup() |>
  select(-c(prob, value))

soiltypes = 
  lap_soiltype |>
  bind_rows(bci_soiltype)

dat %<>%
  left_join(
    soiltypes, 
    by = c('x', 'y', 'fdp')
  )

tally = 
  dat |> 
  drop_na(group) |> 
  filter(dbh > 100) |>
  count(fdp, soiltype) |>
  rename(tally = n)

dtf = 
  dat |> 
  drop_na(group) |> 
  filter(dbh > 100) |>
  count(fdp, group, soiltype) |> 
  left_join(tally, by = c('fdp', 'soiltype')) |>
  mutate(prob = n / tally) |>
  mutate(diagonal = 1 * (group == soiltype))

dtf |>
  group_by(fdp, group) |>
  summarize(
    ratio = sum(prob * diagonal) / sum(prob * !diagonal) * 2
  ) |>
  group_by(fdp) |>
  summarize_at('ratio', list(mean, sd))

plot_bars = 
  dtf |>
  mutate(fdp = ifelse(fdp == 'bci', 'BCI', 'La Planada')) |>
  ggplot(aes(group, prob, group = soiltype, fill = soiltype)) + 
  geom_col(position = 'dodge', color = 'black') +
  facet_wrap(~ fdp) +
  theme(strip.background = element_rect(fill = 'orange')) +
  ylab('proportion of soil type occupancy')


