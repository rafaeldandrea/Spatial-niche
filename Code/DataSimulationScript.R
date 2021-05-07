## =============== libraries ==========
library(tidyverse)
library(tibble)
library(vegan)
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps
library(fpc) ## for clustering analysis
library(numbers) ## for function divisors()
library(MASS) ## for generalized inverse (left inverse) of Cmatrix when S > NST
library(dbscan) ## for HDBSCAN clustering method
library(pracma) ## for function Mode()
library(plyr) ##
library(gdata) ## for function upperTriangle() called in function mean_cosine_rows()


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

## Read scenario from either manual input or command line from Bash file
ind <- make_index(1)

## ====== data folder ===========
dirname <- '~/SpatialNiche/Data/20200615-20scen-NST300S300-H-betanoise-biotime1e5'
setwd(dirname)

## ========= operation directives ==========
landscape.patches <- FALSE
landscape.randomfield <- TRUE


scenarios <- 
  tibble(
    A = 145 ^ 2,
    clustering.method = c('hierarchical'),
    NST = 300,
    S  = 300,
    landscape_autocor_parm  = 10,
    niche_index = seq(0.05,1,by = 0.05),
    run = 1,
    noise_type = 'Beta',
    annual_turnover_rate = .05,
    maxtime = 1e5,
    seed_initial_community = run - 1,
    census_interval = maxtime / 100,
    delta_t = 1,
    initial_biotime = 0
  ) %>%
    group_by(niche_index) %>%
    mutate(
      noise = find_noise_value(niche_index, seed = run, S = S, NST = NST, option = noise_type)
    ) %>%
    ungroup %>%
    rowid_to_column(var = 'scenario') 

## Choose row of 'parameters'data frame to specify scenario
scen = scenarios[ind, ]

## Assign each column of scen to a variable matching its name
list2env(as.list(scen), envir = globalenv())

do.references <- FALSE

Landscape <- Generate_Landscape(seed = run, landscape_autocor_parm = landscape_autocor_parm)

true_C <- Generate_Cmatrix(seed = run,nrow = S, ncol = NST,  noise = noise)

if(noise_type == 'NormalizedNoise') true_C = true_C / rowSums(true_C)

rownames(true_C)=paste0('species_',seq(S))
colnames(true_C)=paste0('resource_',seq(NST))


initial_community = NULL
while(length(unique(initial_community)) != S){
  seed_initial_community <- seed_initial_community + 1
  initial_community <- Generate_Community(Landscape, true_C, seed = seed_initial_community)
}

initial_species <- sort(unique(initial_community))
S <- length(unique(initial_species))

dat <- 
  Simulate_Dynamics(
    Landscape = Landscape, 
    Cmatrix = true_C, 
    Community = initial_community, 
    annual_turnover_rate = annual_turnover_rate, 
    delta_t = delta_t,                          
    biotime = initial_biotime,                          
    maxtime = maxtime,                      
    interval = census_interval,
    seed = run
  )

dat$scenario <- scen

if(hpcc) save(dat, file = paste0('scenario_',ind,'.RData'))

# R <- as.numeric(table(Landscape))
# 
# Community <- dat$Final_Community
# 
# abundance <- sapply(initial_species, function(species) sum(Community == species))
#   
