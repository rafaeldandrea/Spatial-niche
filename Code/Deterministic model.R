## Deterministic Spatial Niche model

## ========== Load functions =========
loadfuns <- {
  lapply(
    list.files(
      path = "~/R_Functions",
      pattern = "[.][R]$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ),
    source
  )
}


RCmod <-
  function(t, y, parms) {
    S = parms$S
    NST = parms$NST
    N = matrix(y, S, NST)
    C = parms$C 
    P = parms$P
    m = parms$m
    b = parms$b
    R = pmax(0, P - colSums(N))
    dNdt = as.numeric(-m * N + b * C * outer(rowSums(N), R))
    return(list(dNdt = dNdt))
  }

## Analytic solution to N*
Nstar <-
  function(C, m, b, P){
    Rstar = m / b * rowSums(solve(C))
    X = t(t(C) * Rstar)
    Nstar = m / b * colSums((P - Rstar) * solve(X))
    return(pmax(0, Nstar))
  }

Rstar <-
  function(C, m, b, P){
    Rstar = m / b * rowSums(solve(C))
    X = t(t(C) * Rstar)
    Nstar = m / b * colSums((P - Rstar) * solve(X))
    return(pmax(0, Rstar))
  }

# Rstar = m / b * as.numeric( solve(C) %*% rep(1, nrow(C)))
# Nstar = m / b * as.numeric()

landscape.patches <- FALSE
landscape.randomfield <- TRUE

the_S = 10
the_NST = 10

scenarios <- {
  crossing(
    A = 120 ^ 2,
    clustering.method = c('hierarchical'),
    NST = the_NST,
    S  = the_S,
    landscape_autocor_parm  = 50,
    niche_index = seq(0.05, 1, by = 0.05),
    run = 1,
    noise_type = c('Beta', 'NormalizedNoise', 'FixedNoise')
  ) %>%
    arrange(run) %>%
    group_by(niche_index, run, noise_type) %>%
    mutate(noise = find_noise_value(
      niche_index,
      seed = run,
      S = the_S,
      NST = the_NST,
      option = noise_type
    )) %>%
    ungroup %>%
    rowid_to_column(var = 'scenario') %>%
    union(
      crossing(
        A = 120 ^ 2,
        clustering.method = c('hierarchical'),
        NST = the_NST,
        S = the_S,
        landscape_autocor_parm = 50,
        niche_index = seq(0.05, 1, by = 0.05),
        run = 1,
        noise_type = c('TruncatedNormal')
      ) %>%
        arrange(run) %>%
        group_by(niche_index, run, noise_type) %>%
        mutate(
          noise  =  find_noise_value(
            niche_index,
            seed = run,
            S = the_S,
            NST = the_NST,
            option = noise_type
          )
        ) %>%
        ungroup %>%
        mutate(scenario = 61:80)
    )
}

maxtime = 1e6

ind = make_index(3)

scen = scenarios[ind,]

list2env(scen, envir = environment())

Landscape <- Generate_Landscape(seed = run, landscape_autocor_parm = landscape_autocor_parm)

P = tabulate(Landscape)

m = 1e-2

b = 1

if(noise_type %in% c('Beta', 'NormalizedNoise')) true_C <- Generate_Cmatrix(seed = run, nrow = S, ncol = NST, noise = noise)
if(noise_type == 'FixedNoise') true_C <- Generate_FixedCmatrix(S = S, NST = NST, noise = noise)
if(noise_type == 'NormalizedNoise') true_C <- true_C / rowSums(true_C)
if(noise_type == 'TruncatedNormal') true_C <- Generate_NormnoiseCmatrix(seed = run, nrow = S, ncol = NST, noise = noise)

# true_C = true_C / rowSums(true_C)

initial_community = NULL
seed_initial_community = run
while(length(unique(initial_community)) != S){
  seed_initial_community <- seed_initial_community + 1
  initial_community <- Generate_Community(Landscape, true_C, seed = seed_initial_community)
}

N0 <-
  tibble(
    soil = Landscape,
    species = initial_community
  ) %>%
  group_by(soil) %>%
  count(species) %>%
  right_join(
    crossing(
      soil = 1:NST,
      species = 1:S),
    by = c("soil", "species")
    ) %>%
  replace_na(list(n = 0)) %>%
  pull(n)

N0 = rep(A / S / NST, S * NST)

out = ode(
  y = N0, 
  times = seq(maxtime), 
  RCmod, 
  parms = list(C = true_C, P = P, b = b, m = m, S = S, NST = NST)
)

Nij = matrix(out[maxtime, -1], S, NST)

years = round(seq(1, maxtime, length = 100))

N_time <- 
  t(apply(out[years, -1], 1, function(n){
    rowSums(matrix(n, S, NST))
  })) %>%
  as_tibble %>%
  mutate(year = years) %>%
  pivot_longer(cols = -year, names_to = 'species', values_to = 'abundance')

R_time <- 
  t(apply(out[years, -1], 1, function(n){
    P - colSums(matrix(n, S, NST))
  })) %>%
  as_tibble %>%
  mutate(year = years) %>%
  pivot_longer(cols = -year, names_to = 'resource', values_to = 'abundance')

plot_N_time <-
  N_time %>%
  ggplot(aes(year, abundance, color = species, group = species)) +
  geom_line() +
  theme(legend.position = 'none') +
  ggtitle('Species abundance by time')

plot_R_time <-
  R_time %>%
  ggplot(aes(year, abundance, color = resource, group = resource)) +
  geom_line() +
  theme(legend.position = 'none') +
  ggtitle('Resource abundance by time')


data <-
  tibble(
    N0 = rowSums(matrix(N0, S, NST)),
    N = round(rowSums(Nij)),
    R = P - colSums(Nij),
    K = rowSums(true_C),
    Kc = colSums(true_C),
    Nstar = Nstar(C = true_C, m = m, b = b, P = P),
    Rstar = Rstar(C = true_C, m = m, b = b, P = P)
  )

data2 <-
  tibble(
    pij = as.numeric(Nij/rowSums(Nij)),
    C = as.numeric(true_C)
  )

plot_N_vs_C <-
  data2 %>%
  ggplot(aes(C * mean(pij) / mean(C), pij)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color = red) +
  labs(
    x = 'Scaled true C matrix',
    y = 'Nij / Ni'
  ) +
  ggtitle('Nij / Ni vs true C matrix (double log scale)')

plot_K <-
  data %>%
  ggplot() +
  geom_hline(yintercept = 720, color = red, linetype = 'dashed') +
  geom_segment(aes(x = K, y = 0, xend = K, yend = N)) +
  #geom_point() +
  geom_point(aes(K, Nstar), color = blue) +
  #scale_y_log10() +
  labs(x = 'Row sum of C', y = 'Species abundance') +
  #ylim(c(0, max(data$N, data$Nstar))) +
  ggtitle('Abundance by row sum')

plot_Kc <-
  data %>%
  ggplot(aes(Kc, N)) +
  geom_hline(yintercept = 720, color = red, linetype = 'dashed') +
  geom_segment(aes(x = Kc, y = 0, xend = Kc, yend = N)) +
  #geom_point() +
  geom_point(aes(Kc, Nstar), color = blue) +
  #scale_y_log10() +
  labs(x = 'Column sum of C', y = 'Species abundance') +
  # ylim(c(0, max(data$N, data$Nstar))) +
  ggtitle('Abundance by column sum')

plot_KK <-
  data %>% 
  ggplot(aes(scale(K), scale(-Kc), color = N)) + 
  geom_point(size = 5) + 
  scale_color_gradient(low = 'lightgrey', high = 'black') +
  labs(
    x = 'Scaled C row sum -- Species intrinsic fitness',
    y = 'Scaled C column sum (x -1) -- Species niche size') +
  ggtitle('Abundance vs Row and Column Sums')

plot_sort <-
  data %>%
  arrange(N) %>%
  rowid_to_column(var = 'species') %>%
  ggplot(aes(species, N)) +
  geom_hline(yintercept = 720, color = red, linetype = 'dashed') +
  geom_segment(aes(x = species, y = 0, xend = species, yend = N)) +
  geom_point(aes(species, Nstar), color = blue) +
  #scale_y_log10() +
  # ylim(c(0, max(data$N, data$Nstar))) +
  labs(x = 'Species rank', y = 'Species abundance') +
  ggtitle('Abundance by rank')

plot_N_pred <-
  data %>%
  ggplot(aes(Nstar, N)) +
  geom_abline(aes(intercept = 0, slope = 1), color = red) +
  geom_point() +
  labs(x = 'Predicted abundance', y = 'Observed abundance') +
  ggtitle('Predicted vs Observed species abundances')

plot_R_pred <-
  data %>%
  ggplot(aes(Rstar, R)) +
  geom_abline(aes(intercept = 0, slope = 1), color = red) +
  geom_point() +
  labs(x = 'Predicted abundance', y = 'Observed abundance') +
  ggtitle('Predicted vs Observed resource abundances')

roundNij = round(Nij)
living_species = (1:S)[rowSums(roundNij) > 1]
true_C_pruned = true_C[living_species, ]
colnames(true_C_pruned) = LETTERS[1:NST]   ## names soil types alphabetically
rownames(true_C_pruned) = living_species


inferred_C = Estimate_C2(Nij = Nij[living_species, ], round = FALSE)
inferred_C = Nij[living_species, ] / rowSums(Nij[living_species, ])
colnames(inferred_C) = LETTERS[1:NST]       ## names soil types alphabetically
rownames(inferred_C) = living_species


plot_eC <-
  as_tibble(inferred_C) %>%
  mutate(species = living_species) %>%
  pivot_longer(cols = -species, names_to = 'resource', values_to = 'Cij') %>%
  ggplot(aes(resource, reorder(species, desc(species)))) +
  geom_raster(aes(fill = Cij)) +
  labs(x = 'Soil type', y = 'Species') +
  theme(legend.position = 'none') +
  ggtitle('Estimated C matrix')


plot_C <-
  as_tibble(true_C_pruned) %>%
  mutate(species = living_species) %>%
  pivot_longer(cols = -species, names_to = 'resource', values_to = 'Cij') %>%
  ggplot(aes(resource, reorder(species, desc(species)))) +
  geom_raster(aes(fill = Cij)) +
  labs(x = 'Soil type', y = 'Species') +
  theme(legend.position = 'none') +
  ggtitle('True C matrix')

colnames(Nij) = LETTERS[1:NST]
plot_Nij <-
  as_tibble(Nij[living_species, ] / rowSums(Nij[living_species, ])) %>%
  mutate(species = living_species) %>%
  pivot_longer(cols = -species, names_to = 'resource', values_to = 'Cij') %>%
  ggplot(aes(resource, reorder(species, desc(species)))) +
  geom_raster(aes(fill = Cij)) +
  labs(x = 'Soil type', y = 'Species') +
  theme(legend.position = 'none') +
  ggtitle('Observed Nij')

niche_index_true = round(1 - mean_cosine_rows(true_C), 2)
niche_index_true_pruned = round(1 - mean_cosine_rows(true_C_pruned), 2)
niche_index_inferred = round(1 - mean_cosine_rows(inferred_C), 2)

gridExtra::grid.arrange(
  plot_N_time,
  plot_R_time,
  # plot_KK,
  # plot_sort, 
  plot_N_pred,
  plot_R_pred, 
  plot_N_vs_C,
  plot_Nij,
  plot_C,
  plot_eC,
  ncol = 4,
  top = paste(
          'Niche index =', 
          round(niche_index, 2), 
          '         Niche index (true C) =',
          niche_index_true,
          '         Niche index (pruned true C) =',
          niche_index_true_pruned,
          '         Niche index (true C) =',
          niche_index_true,
          
          '         C matrix disribution:',
          noise_type,
          '         RV coefficient =',
          round(RV_coefficient(true_C_pruned, inferred_C)$RV, 2)
        )
  )

ls_Nij = Nij[living_species, ]
RHS = true_C_pruned * outer(rowSums(ls_Nij), colSums(ls_Nij) / colSums(rowSums(ls_Nij) * true_C_pruned)) + 1e-10
eRHS = inferred_C * outer(rowSums(ls_Nij), colSums(ls_Nij) / colSums(rowSums(ls_Nij) * inferred_C)) + 1e-10

max(abs(1 - RHS / ls_Nij))
max(abs(1 - eRHS / ls_Nij))
