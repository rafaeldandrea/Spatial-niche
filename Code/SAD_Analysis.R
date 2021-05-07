user = 'rafael'


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

library(dplyr)

datadir <-
  switch(
    7,
    '~/spatialniche/data/20200423', 
    '~/spatialniche/data/20200423-61scen-NST20S20-H', 
    '~/spatialniche/data/20200424-61scen-NST20S20-H-fixednoise', 
    '~/spatialniche/data/20200425-100scen-NST20S20-H-fixednoise-biotime1e5',
    '~/spatialniche/data/20200425-61scen-NST20S20-H-normalizednoise-biotime1e5',
    '~/spatialniche/data/20200429-20scen-NST20S20-H-allnoise-biotime1e6',
    '~/spatialniche/data/20200508-80scen-NST20S20-H-allnoise-biotime1e5'
  )

scenarios <- {
  crossing(
    A = 120 ^ 2,
    clustering.method = c('hierarchical'),
    NST = 20,
    S  = 20,
    landscape_autocor_parm  = 50,
    niche_index  = seq(0.05,1,by=0.05),
    run = 1,
    noise_type = c('Beta', 'NormalizedNoise', 'FixedNoise')
  ) %>%
  arrange(run) %>%
  group_by(niche_index, run, noise_type) %>%
  mutate(noise  =  find_noise_value(niche_index, seed = run, S = 20, NST = 20, option = noise_type)) %>%
  ungroup %>%
  rowid_to_column(var = 'scenario') %>%
  union(
    crossing(
      A = 120 ^ 2,
      clustering.method = c('hierarchical'),
      NST = 20,
      S  = 20,
      landscape_autocor_parm  = 50,
      niche_index  = seq(0.05,1,by=0.05),
      run = 1,
      noise_type = c('TruncatedNormal')
    ) %>%
      arrange(run) %>%
      group_by(niche_index, run, noise_type) %>%
      mutate(noise  =  find_noise_value(niche_index, seed = run, S = 20, NST = 20, option = noise_type)) %>%
      ungroup %>%
      mutate(scenario = 61:80)
  )
}

setwd(datadir)
lf=paste0('scenario_',1:76,'_run_1.RData')
df=NULL
for(char in lf){
  if(file.exists(char)){
    foo=get(load(char))
    oC = foo$C
    eC = foo$est_C
    rowsum_oC = rowSums(oC)[seq(nrow(oC)) %in% as.numeric(foo$observed_community$species)]
    rowsum_eC = rowSums(eC)
    off_diag_C = oC 
    diag(off_diag_C) = NA
    com <-
      # foo$com$abundances %>% 
      foo$abundances %>% 
      rowid_to_column(var = 'time') %>% 
      pivot_longer(cols = -time, names_to = 'species_id', values_to = 'abundance') %>% 
      mutate(
        species_id = factor(species_id, levels = paste0('species_', 1:20)), 
        species = match(species_id, levels(species_id))
      ) %>%
      select(time, species, abundance)
    parms <-
      tibble(
        scenario = foo$scenario, 
        noise = foo$noise, 
        mean_cosine_rows = foo$mean_cosine_rows,
        niche_index = 1 - mean_cosine_rows,
        RV = RV_coefficient(oC[rownames(eC), ], eC)$RV,
        rowsd_oC = apply(off_diag_C, 1, sd, na.rm = TRUE),
        rowsum_oC = rowsum_oC, 
        rowsum_eC = rowsum_eC,
        species = 1:20
      )
    bar <-
      parms %>%
      left_join(com, by = 'species')
    df = rbind(df, bar)
  }
}

burn_in = 50

census_interval = 5


thedata <-
  df %>%
  left_join(
    scenarios %>%
      select(scenario, noise_type),
    by = c('scenario')
  ) %>%
  filter(
    time > burn_in,
    mod(time, census_interval) == 0
  )

df_var=NULL
for(char in lf){
  if(file.exists(char)){
    foo=get(load(char))
    df_var = rbind(df_var, cbind(foo$scenario, foo$noise, foo$mean_cosine_rows, apply(foo$abundances, 1, var)))
  }
}
colnames(df_var) = c('scenario', 'noise', 'mean_cosine_rows', 'var')
df_var <- 
  df_var %>%
  as_tibble %>%
  mutate(niche_index = 1 - mean_cosine_rows)
  

df_ind=NULL
for(char in lf){
  if(file.exists(char)){
    foo = get(load(char))
    n = foo$abundances[-seq(burn_in),]
    n = n[mod(seq(nrow(n)), census_interval) == 0, ]
    df_ind <-
      rbind(
       df_ind, 
       cbind(
         foo$scenario, 
         nrow(foo$abundances), 
         nrow(n), 
         sapply(n, function(v) Box.test(v, type ='Ljung')$p.value)
        )
      )
  }
}
colnames(df_ind) = c('scenario', 'total_snapshots', 'used_snapshots', 'boxtest_pvalue')
df_ind = as_tibble(df_ind)

dat_ind <-
  df_ind %>%
  group_by(scenario) %>%
  summarize(independence = mean(boxtest_pvalue > .05)) %>%
  left_join(
    scenarios %>%
      select(scenario, noise_type, niche_index),
    by = c('scenario')
  )

plot_ind <-
  dat_ind %>%
  ggplot(aes(niche_index, independence, group = noise_type, color = noise_type)) +
  geom_line() +
  geom_point()

gridExtra::grid.arrange(plot_ind)


dat_var <-
  df_var %>%
  left_join(
    scenarios %>%
      select(scenario, noise_type),
    by = c('scenario')
  ) %>%
  group_by(scenario) %>%
  mutate(year = seq(n())) %>%
  ungroup


if(FALSE){
  res <-
    thedata %>%
    filter(niche_index < 1) %>%
    group_by(niche_index) %>%
    summarize(
      # poisson_pval = summary(vcd::goodfit(abundance))[3],
      poisson_pval_fitSAD = fitSAD(abundance, dbn = 'Pois', method = 'CVM')$p.value,
      logseries_pval = fitSAD(abundance, dbn = 'LS', method = 'CVM')$p.value
    )
}

dat_sad <-
  thedata %>%
  group_by(noise_type, niche_index) %>%
  count(abundance) %>%
  mutate(cum = (n + sum(n) - cumsum(n))/sum(n)) %>%
  ungroup

dat_sad_year <-
  thedata %>%
  group_by(noise_type, niche_index, time) %>%
  count(abundance) %>%
  mutate(cum = (n + sum(n) - cumsum(n))/sum(n)) %>%
  ungroup


chosen_ni <- dat_sad$niche_index[sapply(1:19/20, function(value) which.min(abs(dat_sad$niche_index - value)))]

fit_df <-
  dat_sad %>%
  select(noise_type, niche_index, abundance) %>%
  ddply(.(noise_type, niche_index), function(v){
    obs_abuns = min(v$abundance):max(v$abundance)
    ls_pmf = vcdExtra::dlogseries(obs_abuns, prob = LS_p(720))
    ls_cum = ls_pmf + vcdExtra::plogseries(obs_abuns, prob = LS_p(720), lower.tail = FALSE)
    pois_pmf = dpois(obs_abuns, lambda = 720)
    pois_cum = pois_pmf + ppois(obs_abuns, lambda = 720, lower.tail = FALSE)
    tibble(
      obs_abuns = obs_abuns,
      ls_pmf = ls_pmf,
      ls_cum = ls_cum,
      pois_pmf = pois_pmf,
      pois_cum = pois_cum
    )
  }) %>% 
  as_tibble
  
plot_sad_beta <-
  dat_sad %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  ggplot() +
  geom_point(aes(abundance, cum)) +
  geom_line(aes(obs_abuns, ls_cum), color = blue, data = fit_df %>% filter(noise_type == 'Beta', niche_index < 1)) +
  geom_line(aes(obs_abuns, pois_cum), color = red, data = fit_df %>% filter(noise_type == 'Beta', niche_index < 1)) +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  ggtitle('Beta-distributed off-diagonal matrix elements')

plot_sad_beta_year <-
  dat_sad_year %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  ggplot() +
  geom_line(aes(abundance, cum, color = time, group = time)) +
  geom_line(aes(obs_abuns, ls_cum), color = yellow, size = 2, data = fit_df %>% filter(noise_type == 'Beta', niche_index < 1)) +
  geom_line(aes(obs_abuns, pois_cum), color = red, data = fit_df %>% filter(noise_type == 'Beta', niche_index < 1)) +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  ggtitle('Beta-distributed off-diagonal matrix elements')


plot_sad_truncatednormal <-
  dat_sad %>%
  filter(noise_type == 'TruncatedNormal', niche_index < 1) %>%
  ggplot() +
  geom_point(aes(abundance, cum)) +
  geom_line(aes(obs_abuns, ls_cum), color = blue, data = fit_df %>% filter(noise_type == 'TruncatedNormal', niche_index < 1)) +
  geom_line(aes(obs_abuns, pois_cum), color = red, data = fit_df %>% filter(noise_type == 'TruncatedNormal', niche_index < 1)) +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  ggtitle('Truncated Normal-distributed off-diagonal matrix elements')

plot_sad_normalized <-
  dat_sad %>%
  filter(noise_type == 'NormalizedNoise', niche_index < 1) %>%
  ggplot() +
  geom_point(aes(abundance, cum)) +
  geom_line(aes(obs_abuns, ls_cum), color = blue, data = fit_df %>% filter(noise_type == 'NormalizedNoise', niche_index < 1)) +
  geom_line(aes(obs_abuns, pois_cum), color = red, data = fit_df %>% filter(noise_type == 'NormalizedNoise', niche_index < 1)) +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  ggtitle('Beta-distributed off-diagonal matrix elements -- Normalized row sums')

plot_sad_fixed <-
  dat_sad %>%
  filter(noise_type == 'FixedNoise', niche_index < 1) %>%
  ggplot() +
  geom_point(aes(abundance, cum)) +
  geom_line(aes(obs_abuns, ls_cum), color = blue, data = fit_df %>% filter(noise_type == 'FixedNoise', niche_index < 1)) +
  geom_line(aes(obs_abuns, pois_cum), color = red, data = fit_df %>% filter(noise_type == 'FixedNoise', niche_index < 1)) +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  ggtitle('Uniform off-diagonal matrix elements')


# gridExtra::grid.arrange(plot_sad, top = datadir)

plot_hist_beta <-
  thedata %>%
  filter(time > burn_in) %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  ggplot() +
  geom_histogram(aes(abundance), bins = 30) +
  geom_vline(aes(xintercept = 720), color = green, linetype = 'dashed', size = 1) +
  geom_line(aes(obs_abuns, 1e5 * ls_pmf / sum(ls_pmf)), color = blue, size = 1, data = fit_df %>% filter(noise_type == 'Beta', niche_index < 1)) +
  geom_line(aes(obs_abuns, 1e6 * pois_pmf / sum(pois_pmf)), color = red, size = 1, data = fit_df %>% filter(noise_type == 'Beta', niche_index < 1)) +
  scale_x_log10() +
  facet_wrap(~ niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  ggtitle('Beta-distributed off-diagonal matrix elements')

plot_hist_truncatednormal <-
  thedata %>%
  filter(time > burn_in) %>%
  filter(noise_type == 'TruncatedNormal', niche_index < 1) %>%
  ggplot() +
  geom_histogram(aes(log10(abundance))) +
  facet_wrap(~ niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  ggtitle('Truncated Normal-distributed off-diagonal matrix elements')

plot_hist_normalized <-
  thedata %>%
  filter(time > burn_in) %>%
  filter(noise_type == 'NormalizedNoise', niche_index < 1) %>%
  ggplot() +
  geom_histogram(aes(log10(abundance))) +
  facet_wrap(~ niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  ggtitle('Beta-distributed off-diagonal matrix elements -- Normalized row sums')

plot_hist_fixed <-
  thedata %>%
  filter(time > burn_in) %>%
  filter(noise_type == 'FixedNoise', niche_index < 1) %>%
  ggplot() +
  geom_histogram(aes(log10(abundance))) +
  facet_wrap(~ niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  ggtitle('Uniform off-diagonal matrix elements')

dat_RV <-
  thedata %>%
  select(noise_type, niche_index, RV) %>%
  unique

plot_RV_beta <-
  dat_RV %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  ggplot(aes(niche_index, RV)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylim(c(0, 1)) +
  ggtitle('Beta-distributed off-diagonal matrix elements')

plot_RV_truncatednormal <-
  dat_RV %>%
  filter(noise_type == 'TruncatedNormal', niche_index < 1) %>%
  ggplot(aes(niche_index, RV)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylim(c(0, 1)) +
  ggtitle('Truncated Normal-distributed off-diagonal matrix elements')

plot_RV_normalized <-
  dat_RV %>%
  filter(noise_type == 'NormalizedNoise', niche_index < 1) %>%
  ggplot(aes(niche_index, RV)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylim(c(0, 1)) +
  ggtitle('Beta-distributed off-diagonal matrix elements -- Normalized by row sums')

plot_RV_fixed <-
  dat_RV %>%
  filter(noise_type == 'FixedNoise', niche_index < 1) %>%
  ggplot(aes(niche_index, RV)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylim(c(0, 1)) +
  ggtitle('Uniform off-diagonal matrix elements')

# gridExtra::grid.arrange(plot_RV, top = datadir)

plot_n_rsO <-
  thedata %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  group_by(niche_index, abundance) %>%
  summarize(
    errbar = sd(rowsum_oC),
    rowsum_oC = mean(rowsum_oC)
  ) %>%
  ggplot(aes(abundance, rowsum_oC)) +
  geom_errorbar(aes(ymin = rowsum_oC - errbar, ymax = rowsum_oC + errbar)) +
  geom_point() +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = 'Species abundance', y = 'Row sums of true C matrix') +
  ggtitle('Beta-distributed off-diagonal matrix elements')

plot_n_rsE <-
  thedata %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  group_by(niche_index, abundance) %>%
  summarize(
    errbar = sd(rowsum_eC),
    rowsum_eC = mean(rowsum_eC)
  ) %>%
  ggplot(aes(abundance, rowsum_eC)) +
  geom_errorbar(aes(ymin = rowsum_eC - errbar, ymax = rowsum_eC + errbar)) +
  geom_point() +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = 'Species abundance', y = 'Row sums of estimated C matrix') +
  ggtitle('Beta-distributed off-diagonal matrix elements')
  
plot_rsO_rsE <-
  thedata %>%
  filter(noise_type == 'Beta', niche_index < 1) %>%
  ggplot(aes(rowsum_eC, rowsum_oC)) +
  geom_point() +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2), scales = 'free') +
  labs(x = 'Row sums of estimated C matrix', y = 'Row sums of true C matrix') +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle('Beta-distributed off-diagonal matrix elements')
  

# gridExtra::grid.arrange(plot_cor + scale_x_log10(), top = datadir)


plot_var_fixed <-
  dat_var %>%
  filter(noise_type == 'FixedNoise', year > burn_in) %>%
  ggplot(aes(year, sqrt(var)/720)) +
  geom_line() +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2)) +
  labs(x = 'Time (1k years)', y = 'Abundances coefficient of variation') +
  ggtitle('Uniform off-diagonal matrix elements')

plot_var_beta <-
  dat_var %>%
  filter(noise_type == 'Beta', year > burn_in) %>%
  ggplot(aes(year, sqrt(var)/720)) +
  geom_line() +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2)) +
  labs(x = 'Time (1k years)', y = 'Abundances coefficient of variation') +
  ggtitle('Beta-distributed off-diagonal matrix elements')

plot_var_normalized <-
  dat_var %>%
  filter(noise_type == 'NormalizedNoise', year > burn_in) %>%
  ggplot(aes(year, sqrt(var)/720)) +
  geom_line() +
  facet_wrap(~niche_index, labeller = function(x) round(x, 2)) +
  labs(x = 'Time (1k years)', y = 'Abundances coefficient of variation') +
  ggtitle('Beta-distributed off-diagonal matrix elements -- Normalized row sums')


plot_Cvar <-
  thedata %>%
  group_by(noise_type, niche_index) %>%
  summarize(
    var_rsO = var(rowsum_oC),
    var_rsE = var(rowsum_eC)
  ) %>%
  filter(noise_type == 'Beta') %>%
  pivot_longer(cols = c(var_rsO, var_rsE), names_to = 'matrix', values_to = 'var_rs') %>%
  ggplot(aes(niche_index, var_rs, group = matrix, color = matrix)) +
  geom_point() +
  geom_smooth(se = FALSE)


analysis <-
  thedata %>% 
  filter(time > burn_in, niche_index < 1) %>%
  ddply(.(noise_type, niche_index, species), function(v){
   print(paste(unique(v$noise_type), round(unique(v$niche_index), 2), unique(v$species)))
   ls_mod = fitSAD(v$abundance, dbn = 'LS', method = 'CVM')
   pois_mod = fitSAD(v$abundance, dbn = 'Pois', method = 'CVM')
   return(
     tibble(
       mean_abun = mean(v$abundance),
       ls_stat = ls_mod$statistic,
       ls_pval = ls_mod$p.value,
       pois_stat = pois_mod$statistic,
       pois_pval = pois_mod$p.value
     )
   )
  })


analysis_summary <-
  analysis %>%
  left_join(thedata, by = c('noise_type', 'niche_index', 'species')) %>%
  group_by(noise_type, niche_index) %>%
  summarize(
    ls_stat = mean(ls_stat),
    pois_stat = mean(pois_stat),
    ls_hits = mean(ls_pval > .05),
    pois_hits = mean(pois_pval > .05)
  )
  
plot_temporal_statistic <-
  analysis_summary %>%
  select(noise_type, niche_index, ls_stat, pois_stat) %>%
  pivot_longer(cols = c(ls_stat, pois_stat), values_to = 'statistic', names_to = 'distribution') %>%
  ggplot(aes(niche_index, -statistic, group = distribution, color = distribution)) +
  geom_smooth(se = FALSE) +
  geom_point() +
  facet_wrap(~ noise_type) +
  labs(x = 'Niche index', y = 'MLE fit deficit') +
  ggtitle('Average statistic of Logseries / Poisson fit to species dynamics')

plot_temporal_pval <-
  analysis_summary %>%
  select(noise_type, niche_index, ls_hits, pois_hits) %>%
  pivot_longer(cols = c(ls_hits, pois_hits), values_to = 'statistic', names_to = 'distribution') %>%
  ggplot(aes(niche_index, statistic, group = distribution, color = distribution)) +
  geom_smooth(se = FALSE) +
  geom_point() +
  facet_wrap(~ noise_type) +
  labs(x = 'Niche index', y = 'Proportion of non-rejected Logseries / Poisson fits') +
  ggtitle('Proportion of species following Logseries / Poisson dynamics')

gridExtra::grid.arrange(plot_temporal_statistic, plot_temporal_pval, nrow = 2)

plot_temporal_abundance <-
  analysis %>%
  filter(noise_type == 'Beta', mean_abun < 10) %>%
  select(-c(ls_pval, pois_pval)) %>%
  pivot_longer(
    cols = -c(noise_type, niche_index, species, mean_abun), 
    names_to = 'distribution', 
    values_to = 'statistic'
  ) %>%
  ggplot(aes(mean_abun, -statistic, group = distribution, color = distribution)) +
  geom_smooth(se = FALSE) +
  geom_point() +
  scale_x_log10() +
  labs(x = 'Mean abundance', y = 'Fit deficit') +
  facet_wrap(~ niche_index, labeller = function(x) round(x ,2), scales = 'free') +
  ggtitle('Beta-distributed off-diagonal matrix elements
           Statistic of Logseries / Poisson fit to species dynamics by abundance')

