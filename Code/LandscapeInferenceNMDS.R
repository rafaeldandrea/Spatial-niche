
## libraries ==========
library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps
library(fpc) ## for clustering analysis
library(numbers) ## for function divisors()
library(MASS) ## for generalized inverse (left inverse) of Cmatrix when S > NST
library(dbscan) ## for HDBSCAN clustering method
library(pracma) ## for function Mode()


## Reader parameters ==============
## Indicator variable; 1 if running on High Power Computing cluster, 0 if running on my PC.
hpcc <- !(user <- Sys.info()["user"]) %in% c("rafael", "wangdz")


## Read scenario from either input on command line (if running on hpcc)
## or chosen scenario listed in the line below (if running on PC)
simind <- ifelse(hpcc, as.numeric(commandArgs(TRUE)[1]), 1601)

## Load functions =========
loadfuns <- {
  lapply(
    list.files(
      path = switch(
        match(user, c("rafael", "wangdz", "rdandrea", "dw9")),
        "~/R_Functions",
        "/Users/wangdz/Documents/UIUC/Odywer/R_Functions",
        "/data/rdandrea-dw9/R_Functions",
        "/data/rdandrea-dw9/R_Functions"
      ),
      pattern = "[.][R]$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = TRUE
    ),
    source
  )
}

## operation directives ==========
landscape.patches <- FALSE
landscape.randomfield <- TRUE

clustering.method <- "pam"

do.kmeansgap <- (clustering.method == "kmeans")
do.pamk <- (clustering.method == "pam")
do.hierarch <- (clustering.method == "hierarch")
do.NMDS <- (clustering.method == "nmds")
do.hdbscan <- (clustering.method == "hdbscan")

plot.landscape_community_C_matrices <- FALSE
plot.landscape_community_variogram <- FALSE
plot.clustering.statistics.against.cell.size <- FALSE

analyze.data <- FALSE
analyze.data2 <- FALSE

calculate.likelihood <- FALSE
simulate.likelihood <- FALSE
simulate.likelihood2 <- TRUE



## scenarios ===========
scenarios <- {
  crossing(
    tibble(
      S = c(10, 100), ## number of species
      NST = c(5, 10) ## number of soil types
    ),
    tibble(
      A = 120^2, ## number of sites
      NC = A / numbers::divisors(sqrt(A))^2, ## number of cells
      non_preferred_site_affinity = 0.1, ## affinity for non-preferred sites (between 0 and 1)
    ),
    tibble(
      C_sd = c(0, .01, .05, .1)
    ),
    tibble(
      landscape_autocor_parm = c(15, 30, 100)
    )
  ) %>%
    mutate(
      mode = "spec",
      C_mean = NA,
      C_cv = NA
    ) %>%
    rbind(
      crossing(
        tibble(
          S = c(10, 100), ## number of species
          NST = c(10, 100) ## number of soil types
        ),
        tibble(
          A = 120^2, ## number of sites
          NC = A / numbers::divisors(sqrt(A))^2, ## number of cells
          non_preferred_site_affinity = NA, ## affinity for non-preferred sites (between 0 and 1)
        ),
        tibble(
          mode = "gen",
          C_mean = 1,
          C_cv = c(.2, .4, .6),
          C_sd = NA
        ),
        tibble(
          landscape_autocor_parm = c(15, 30, 100)
        )
      )
    ) %>%
    filter(NC > NST & NC < A) %>% ## remove values of NC that do not make sense given A and NSTrd
    crossing(tibble(run = 1:100)) %>% ## seed for the random number generator
    arrange(desc(mode), desc(landscape_autocor_parm), C_sd, S, NC)
}

simind_vec <- 1 + (simind - 1) * 748 + 0:747
simind_vec <- simind_vec[simind_vec <= nrow(scenarios)]

## functions
{
  
  ReorderMatrix=function(data){
    order_matrix=NST+1-t(apply(data,1,rank))
    preferred=t(apply(order_matrix,1,order))
    mylist=NULL
    for(i in 1:nrow(data)){ 
      for(j in 1:ncol(data)){ 
        if(!preferred[i,j] %in% mylist){ 
          mylist=c(mylist,preferred[i,j])
          break
        }
      }
    }
    data=data[,mylist]
    return(data)
  }
  
  ## Landscape: vector whose A elements are the soil types of each site ========
  Generate_Landscape <- function(seed,landscape_autocor_parm) {
    if (landscape.patches) {
      patches <- Patches(side = sqrt(A / NC), number = NC, noise = 0, L = sqrt(A))
      set.seed(seed) ## change the parameter 'seed' to obtain different landscapes
      Landscape <- sample(NST, replace = TRUE, size = NC)[match(patches, unique(patches))]
    }
    if (landscape.randomfield) {
      ## create a Gaussian random field
      field <- RandomField(L = sqrt(A), rangepar = landscape_autocor_parm, sillpar = 1, nuggetpar = 0, seed = seed)
  
      ## simplify the field by boxing sites into one of NST soil types
      Landscape <- 1 + findInterval(field, quantile(field, seq(NST) / NST), rightmost.closed = TRUE)
    }
    return(Landscape)
  }
  
  ## define cell grid based on NC. (Uses pos2coord from R_Functions) ============
  Generate_Grid <- function(A, NC) {
    cellgrid <-
      as_tibble(pos2coord(pos = seq(A), dim.mat = c(sqrt(A), sqrt(A)))) %>%
      mutate(
        cx = 1 + (V1 - 1) %/% (max(V1) / sqrt(NC)),
        cy = 1 + (V2 - 1) %/% (max(V2) / sqrt(NC)),
        cellname = paste0(
          formatC(cx, width = nchar(sqrt(NC)), flag = "0"),
          formatC(cy, width = nchar(sqrt(NC)), flag = "0")
        ),
        cell = match(cellname, unique(cellname))
      ) %>%
      pull(cell)
    return(cellgrid)
  }
  
  ## affinity matrix Cmatrix ========
  Generate_Cmatrix <- function(seed,mode='spec') {
    set.seed(seed)
    if (mode == "spec") {
      foobar <- c(1, rep(non_preferred_site_affinity, NST - 1))
      Cmatrix <- matrix(non_preferred_site_affinity, S, NST)
      for (i in 1:S) Cmatrix[i, ] <- guyrot(v = foobar, n = i - 1)
      Cmatrix <- ReorderMatrix(Cmatrix)
      Cmatrix[Cmatrix != 1] <- Cmatrix[Cmatrix != 1] + rnorm(length(Cmatrix[Cmatrix != 1]), mean = 0, sd = C_sd)
      Cmatrix[Cmatrix < 0] <- 0
    }
    if (mode == "gen") {
      Cmatrix <- matrix(runif(NST * S, min = C_mean * (1 - sqrt(3) / 2 * C_cv), max = C_mean * (1 + sqrt(3) / 2 * C_cv)), S, NST)
      Cmatrix[Cmatrix < 0] <- 0
      Cmatrix <- ReorderMatrix(Cmatrix)
    }
    return(Cmatrix)
  }
  
  ## community matrix [sqrt(A) x sqrt(A)] =======
  Generate_Community <- function(Landscape, Cmatrix, seed) {
    set.seed(seed)
    R <- sapply(seq(NST), function(soiltype) sum(Landscape == soiltype))
    Community <-
      sapply(
        seq(A),
        function(site) {
          # sample(S,size=1,prob=pmax(0,rowSums(ginv(Cmatrix)))*Cmatrix[seq(S),Landscape[site]])
          sample(S, size = 1, prob = R[Landscape[site]] * Cmatrix[seq(S), Landscape[site]])
        }
      )
    return(Community)
  }
  
  ## Tabulation matrix [S x NC] ======
  Generate_Tabulation <- function(cellgrid, Community) {
    Tabulation <-
      tibble(
        cell = cellgrid,
        species = Community
      ) %>%
      group_by(cell) %>%
      dplyr::count(species) %>%
      spread(key = species, value = n, sep = "_") %>%
      ungroup()
    Tabulation[is.na(Tabulation)] <- 0
    return(Tabulation)
  }
  
  ## estimate Cmatrix matrix =========
  Generate_estC <- function(est_Nij, R) {
    est_C <- t(matrix(as.numeric(unlist(est_Nij[, -1])), nrow(est_Nij), ncol(est_Nij) - 1))
    est_C <- est_C / rowSums(est_C) / R
    # est_C=round(est_C/max(est_C),1)
    #est_C <- ReorderMatrix(est_C)
    return(est_C)
  }
  
  ## find the point in a curve with the maximum |2nd derivative| ====
  Max2Deriv <- function(x, y) {
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    indices <- which.max(abs(d2))
    return(indices)
  }
  
}

# for(simind in simind_vec){
#   try({

if (FALSE) {
  ## scenario for this sim
  scen=scenarios[simind,]
  list2env(scen,envir=.GlobalEnv)
  
  ## Generate landscape, cell grid, C matrix, Community map, tabulation (Nij) =======
  Landscape <- Generate_Landscape(seed = run,landscape_autocor_parm = 15)
  R <- as.numeric(table(Landscape))
  cellgrid <- Generate_Grid(A, NC)
  Cmatrix <- Generate_Cmatrix(seed = run)
  Community <- Generate_Community(Landscape, Cmatrix, seed = run)
  Tabulation <- Generate_Tabulation(cellgrid, Community)

  ## Clustering algorithm ==========
  ## NMDS
  if (do.NMDS) {
    as_matrix <- function(x) {
      if (!tibble::is_tibble(x)) stop("x must be a tibble")
      y <- as.matrix.data.frame(x[, -1])
      rownames(y) <- x[[1]]
      y
    }
    numdimensions <- min(S, NC - 1)
    matrix_Tabulation <- as_matrix(Tabulation)
    rownames(matrix_Tabulation) <- paste("c", 1:NC, sep = "")
    NMDS <- metaMDS(matrix_Tabulation, k = numdimensions, trymax = 100, trace = 0, autotransform = FALSE)
    df_cells <- as_tibble(NMDS$points) %>% mutate(cell = seq(nrow(.)))
    df_species <- as_tibble(NMDS$species) %>% mutate(species = seq(nrow(.)))
  }

  ## KmeansGap
  if (do.kmeansgap) {
    ## function for multi-d kmeans
    KmeansGapMD <- function(dat, traitcols, kvec = 1:(nrow(dat) / 2), itermax = 100, nstart = 1e4, seed) {
      set.seed(seed)
      Wk <- sapply(
        kvec,
        function(k) {
          kmeans(dat[, traitcols], centers = k, iter.max = itermax, nstart = nstart)$tot.within
        }
      )
      return(tibble(k = kvec, Wk = Wk))
    }

    ## null communities
    numnulls <- 1e2
    n.start <- 1e3
    new.nulls <- TRUE
    if (new.nulls) {
      null_species <- null_cells <- NULL
      set.seed(0)
      for (dummy in seq(numnulls)) {
        bar <-
          tibble(
            cell = cell,
            soiltype = Landscape,
            species = sample(Community)
          ) %>%
          group_by(cell) %>%
          dplyr::count(species) %>%
          spread(key = species, value = n, sep = "_") %>%
          ungroup()

        foo <- metaMDS(bar %>% select(-cell), k = numdimensions, trymax = 100, trace = 0)
        null_species <- rbind(null_species, as_tibble(foo$species) %>% mutate(index = dummy))
        null_cells <- rbind(null_cells, as_tibble(foo$points) %>% mutate(index = dummy))
      }
    }

    dat_species <-
      df_species %>%
      select(-species) %>%
      mutate(index = 0) %>%
      rbind(null_species)

    dat_cells <-
      df_cells %>%
      select(-cell) %>%
      mutate(index = 0) %>%
      rbind(null_cells)

    res_species <- ddply(dat_species, .(index), function(df) {
      KmeansGapMD(dat = df, traitcols = seq(numdimensions), nstart = n.start, seed = simind)
    }) %>%
      as_tibble()

    res_cells <- ddply(dat_cells, .(index), function(df) {
      KmeansGapMD(dat = df, traitcols = seq(numdimensions), nstart = n.start, seed = simind)
    }) %>%
      as_tibble()


    res_cells_summary <-
      res_cells %>%
      filter(index > 0) %>%
      mutate(logW = log(Wk)) %>%
      group_by(k) %>%
      summarize_at("logW", list(nullmean = mean, sd = sd)) %>%
      left_join(res_cells %>% filter(index == 0) %>% select(-index), by = "k") %>%
      mutate(gap = nullmean - log(Wk))

    res_species_summary <-
      res_species %>%
      filter(index > 0) %>%
      mutate(logW = log(Wk)) %>%
      group_by(k) %>%
      summarize_at("logW", list(nullmean = mean, sd = sd)) %>%
      left_join(res_species %>% filter(index == 0) %>% select(-index), by = "k") %>%
      mutate(gap = nullmean - log(Wk))

    gapplot_species <-
      res_species_summary %>%
      ggplot(aes(k, gap)) +
      geom_ribbon(aes(ymin = gap - sd, ymax = gap + sd, fill = red)) +
      geom_point() +
      geom_line() +
      theme(aspect.ratio = 1) +
      ggtitle("Species") +
      theme(legend.position = "none")

    gapplot_cells <-
      res_cells_summary %>%
      ggplot(aes(k, gap)) +
      geom_ribbon(aes(ymin = gap - sd, ymax = gap + sd, fill = red)) +
      geom_point() +
      geom_line() +
      theme(aspect.ratio = 1) +
      ggtitle("Cells") +
      theme(legend.position = "none")

    numclusters_cells <- res_cells_summary$k[which.max(res_cells_summary$gap)]
    numclusters_species <- res_species_summary$k[which.max(res_species_summary$gap)]

    df_cells <-
      df_cells %>%
      mutate(soiltype = kmeans(df_cells %>% select(-cell), centers = numclusters_cells, nstart = 1e3)$cluster)

    df_species <-
      df_species %>%
      mutate(niche = kmeans(df_species %>% select(-species), centers = numclusters_species, nstart = 1e3)$cluster)

    ## estimate matrix Nij based on the identification of cells to soil types
    est_Nij <-
      Tabulation %>%
      left_join(df_cells %>% select(cell, soiltype), by = "cell") %>%
      select(-cell) %>%
      group_by(soiltype) %>%
      summarize_all(sum)
  }

  ## PAM
  if (do.pamk) {
    if (FALSE) {
      pamk_cells <- fpc::pamk(df_cells %>% select(1:numdimensions), krange = 2:min(S, NC / 2), diss = FALSE)
      numclusters_cells <- pamk_cells$nc
      df_cells <-
        df_cells %>%
        mutate(soiltype = pamk_cells$pamobject$clustering)

      pamk_species <- fpc::pamk(df_species %>% select(1:numdimensions), krange = 2:(S - 1), diss = FALSE)
      numclusters_species <- pamk_species$nc
      df_species <-
        df_species %>%
        mutate(niche = pamk_species$pamobject$clustering)

      ## estimate matrix Nij based on the identification of cells to soil types
      est_Nij <-
        Tabulation %>%
        left_join(df_cells %>% select(cell, soiltype), by = "cell") %>%
        select(-cell) %>%
        group_by(soiltype) %>%
        summarize_all(sum)

      clustering.statistic <- as.numeric(clustermodel$pamobject$objective[1])
    }

    ## estimate matrix Nij based on the identification of cells to soil types
    # dis_sites <- vegdist(Tabulation %>% select(-1),method="bray")
    # clustermodel2=pamk(dis_sites,krange=2:min(S,NC-1),diss=TRUE)
    clustermodel <- pamk(Tabulation %>% select(-1), krange = 2:min(S, NC - 1), diss = FALSE)
    numclusters_cells <- clustermodel$nc
    cluster_index <- clustermodel$pamobject$clustering
    clustering.statistic <- as.numeric(clustermodel$pamobject$objective[2])
  }

  ## Hierarchical clustering
  if (do.hierarch) {
    # clustering sites
    dis_sites <- vegdist(Tabulation %>% select(-1), method = "jaccard")
    clus_sites <- hclust(dis_sites)
    plot(clus_sites)
    rect.hclust(clus_sites, NST)
    grp_sites <- cutree(clus_sites, NST)
    # clustering species
    dis_species <- vegdist(t(Tabulation %>% select(-1)), method = "jaccard")
    clus_species <- hclust(dis_species)
    plot(clus_species)
    rect.hclust(clus_species, NST)
    grp_species <- cutree(clus_species, NST)
  }

  ## HDBSCAN
  if (do.hdbscan) {
    dis_sites <- vegdist(Tabulation %>% select(-1), method = "jaccard")
    minPts_vec <- 5:max(5, min(NC / NST, 30))
    minPts_vec <- 2:30
    score <- t(sapply(minPts_vec, function(n) {
      cl <- hdbscan(Tabulation[, -1], minPts = n, xdist = dis_sites)
      return(c(
        badness = sum(cl$membership_prob < .05, na.rm = TRUE) / sum(!is.na(cl$membership_prob)),
        numclusters = max(cl$cluster),
        minPts = n
      ))
    }))
    score <- as_tibble(score) %>% filter(numclusters > 0)
    if (nrow(score) > 0) {
      numclusters_cells <- Mode(score$numclusters)
      minPts <- score$minPts[which(score$numclusters == numclusters_cells)[1]]
      clusteringmod <- hdbscan(Tabulation[, -1], minPts = minPts, xdist = dis_sites)
      cluster_index <- clusteringmod$cluster
      badness <- sum(clusteringmod$membership_prob < .05, na.rm = TRUE) / sum(!is.na(clusteringmod$membership_prob))
    } else {
      numclusters_cells <- 0
      cluster_index <- rep(1, NC)
      badness <- NA
    }
    clustering.statistic <- badness
  }

  est_Nij <-
    Tabulation %>%
    mutate(soiltype = cluster_index) %>%
    select(-cell) %>%
    group_by(soiltype) %>%
    summarize_all(sum) %>%
    filter(soiltype > 0)

  ## Estimate C ========
  est_C <- Generate_estC(est_Nij, R)
  if (!is.matrix(est_C)) est_C <- matrix(est_C, length(est_C), 1)

  ## RV index -- cosine between Cmatrix and est_C =========
  tr <- function(m) sum(diag(m))
  squarify <- function(m) m %*% t(m)
  center <- function(m) t(t(m) - colMeans(m))
  CC <- squarify(center(Cmatrix))
  eCeC <- squarify(center(est_C))
  RV <- tr(CC %*% eCeC) / sqrt(tr(CC %*% CC) * tr(eCeC %*% eCeC))

  if (!hpcc) {
    writeLines(paste("numclusters_cells =", numclusters_cells))
    writeLines(paste("RV coefficient =", round(RV, 3)))
  }


  ## plots ================
  ## plot landscape and community variograms
  if (!hpcc & plot.landscape_community_variogram) {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
    lanvar <- Variogram(stress = Landscape, plot = TRUE, cutoff = sqrt(A) * sqrt(2) / 3 * .9)
    comvar <- Variogram(stress = Community, plot = TRUE, cutoff = sqrt(A) * sqrt(2) / 3 * .9)
  }

  ## plot landscape and community plus Cmatrix matrices
  if (!hpcc & plot.landscape_community_C_matrices) {
    theme_set(theme_bw())

    colmin <- "darkgreen"
    colmax <- "gold"

    plot_Landscape <-
      as_tibble(matrix(Landscape, sqrt(A), sqrt(A))) %>%
      `colnames<-`(formatC(1:sqrt(A), width = 3, flag = "0")) %>%
      mutate(x = formatC(1:sqrt(A), width = 3, flag = "0")) %>%
      gather("y", "value", -x) %>%
      ggplot(aes(x, y, fill = value)) +
      geom_raster() +
      geom_hline(yintercept = c(0, seq(sqrt(NC)) * sqrt(A / NC)), size = 1) +
      geom_vline(xintercept = c(0, seq(sqrt(NC)) * sqrt(A / NC)), size = 1) +
      theme(aspect.ratio = 1) +
      theme(legend.position = 1) +
      scale_fill_gradient(low = colmin, high = colmax) +
      scale_y_discrete(
        breaks = as.character(formatC(seq(0, sqrt(A), by = 20), width = 3, flag = "0")),
        labels = as.character(seq(0, sqrt(A), by = 20))
      ) +
      scale_x_discrete(
        breaks = as.character(formatC(seq(0, sqrt(A), by = 20), width = 3, flag = "0")),
        labels = as.character(seq(0, sqrt(A), by = 20))
      ) +
      ggtitle(paste0("Landscape  (", NST, " soil types)"))

    plot_Community <-
      as_tibble(matrix(Community, sqrt(A), sqrt(A))) %>%
      `colnames<-`(formatC(1:sqrt(A), width = 3, flag = "0")) %>%
      mutate(x = formatC(1:sqrt(A), width = 3, flag = "0")) %>%
      gather("y", "value", -x) %>%
      ggplot(aes(x, y, fill = value)) +
      geom_raster() +
      #geom_hline(yintercept = c(0, seq(sqrt(NC)) * sqrt(A / NC)), size = 1) +
      #geom_vline(xintercept = c(0, seq(sqrt(NC)) * sqrt(A / NC)), size = 1) +
      theme(aspect.ratio = 1) +
      theme(legend.position = 1) +
      scale_fill_gradient(low = colmin, high = colmax) +
      scale_y_discrete(
        breaks = as.character(formatC(seq(0, sqrt(A), by = 20), width = 3, flag = "0")),
        labels = as.character(seq(0, sqrt(A), by = 20))
      ) +
      scale_x_discrete(
        breaks = as.character(formatC(seq(0, sqrt(A), by = 20), width = 3, flag = "0")),
        labels = as.character(seq(0, sqrt(A), by = 20))
      ) +
      ggtitle('Community')
      #ggtitle(paste0("Community  (", S, " species)"))

    if (do.NMDS) {
      plot_MDS_cells <-
        df_cells %>%
        mutate(soiltype = factor(soiltype)) %>%
        ggplot(aes(MDS1, MDS2, group = soiltype, color = soiltype)) +
        geom_point() +
        theme(aspect.ratio = 1) +
        theme(legend.position = "none") +
        ggtitle(paste("Cell MDS", ",  Estimated #soil types =", numclusters_cells))

      plot_MDS_species <-
        df_species %>%
        mutate(niche = factor(niche)) %>%
        ggplot(aes(MDS1, MDS2, group = niche, color = niche)) +
        geom_point() +
        theme(aspect.ratio = 1) +
        theme(legend.position = "none") +
        ggtitle(paste("Species MDS", ",  Estimated #niches =", numclusters_species))
    }

    plot_C <-
      as_tibble(Cmatrix) %>%
      `colnames<-`(formatC(1:NST, width = 3, flag = "0")) %>%
      mutate(species = formatC(1:S, width = 3, flag = "0")) %>%
      gather("soiltype", "value", -species) %>%
      ggplot(aes(species, soiltype, fill = value)) +
      geom_raster() +
      scale_fill_gradient(low = colmin, high = colmax) +
      theme(aspect.ratio = 1) +
      theme(legend.position = 1) +
      scale_y_discrete(
        # breaks=formatC(round(seq(1,NST,length=5)),width=3,flag="0"),
        labels = as.character(1:ncol(Cmatrix))
      ) +
      scale_x_discrete(
        # breaks=formatC(round(seq(1,S,length=5)),width=3,flag="0"),
        labels = as.character(1:nrow(Cmatrix))
      ) +
      ggtitle("Original Cmatrix matrix")

    plot_est_C <-
      as_tibble(est_C) %>%
      `colnames<-`(formatC(1:ncol(est_C), width = 3, flag = "0")) %>%
      mutate(species = formatC(1:nrow(est_C), width = 3, flag = "0")) %>%
      gather("soiltype", "value", -species) %>%
      ggplot(aes(species, soiltype, fill = value)) +
      geom_raster() +
      scale_fill_gradient(low = colmin, high = colmax) +
      theme(aspect.ratio = 1) +
      theme(legend.position = 1) +
      scale_y_discrete(
        # breaks=formatC(round(seq(1,numclusters_cells,length=5)),width=3,flag="0"),
        labels = as.character(1:ncol(est_C))
      ) +
      scale_x_discrete(
        # breaks=formatC(round(seq(1,S,length=5)),width=3,flag="0"),
        labels = as.character(1:nrow(est_C))
      ) +
      ggtitle("Estimated Cmatrix matrix")

    ## show cluster results, and MDS result
    if (do.kmeansgap) {
      gridExtra::grid.arrange(gapplot_species, gapplot_cells, plot_MDS_species, plot_MDS_cells, nrow = 2)
    }
    if (do.pamk & do.NMDS) {
      gridExtra::grid.arrange(
        plot_Landscape,
        plot_MDS_cells,
        plot_C,
        plot_Community,
        plot_MDS_species,
        plot_est_C,
        nrow = 2
      )
    }

    if (do.pamk & !do.NMDS) {
      gridExtra::grid.arrange(
        plot_Landscape,
        plot_C,
        plot_Community,
        plot_est_C,
        nrow = 2
      )
    }
  }

  ## plot a) estimated number of soil types, b) RV coefficient, and c) proportion of runs with correct
  ## estimate of soil types against cell size
  if (!hpcc & plot.clustering.statistics.against.cell.size) {
    results_df <-
      results %>%
      as_tibble() %>%
      mutate(
        index = rep(c("numclusters", "RV"), nrow(results) / 2),
        run = rep(seq(nrow(results) / 2), each = 2)
      ) %>%
      gather(key = "numcells", value = "value", -c(index, run)) %>%
      mutate(numcells = cellsizes[match(numcells, unique(numcells))]) %>%
      spread(index, value)

    results_summary <-
      results_df %>%
      group_by(numcells) %>%
      summarize_at(c("numclusters", "RV"), .funs = list(mean = mean, sd = sd, median = median)) %>%
      mutate(cellsize = sqrt(A / numcells))


    plot <-
      results_summary %>%
      ggplot() +
      geom_errorbar(
        aes(
          x = cellsize,
          ymin = numclusters_mean - numclusters_sd,
          ymax = numclusters_mean + numclusters_sd
        ),
        color = "gray"
      ) +
      geom_line(aes(cellsize, numclusters_mean)) +
      geom_point(aes(cellsize, numclusters_mean)) +
      scale_x_log10() +
      xlab("Cell size") +
      ylab("Estimated number of soil types")

    plot2 <-
      results_summary %>%
      ggplot() +
      geom_errorbar(
        aes(
          x = cellsize,
          ymin = RV_mean - RV_sd,
          ymax = RV_mean + RV_sd
        ),
        color = "gray"
      ) +
      geom_line(aes(cellsize, RV_mean)) +
      geom_point(aes(cellsize, RV_mean)) +
      scale_x_log10() +
      xlab("Cell size") +
      ylab("RV coefficient")

    plot3 <-
      results_df %>%
      group_by(numcells) %>%
      summarize(hits = sum(numclusters == 10) / n()) %>%
      mutate(cellsize = sqrt(A / numcells)) %>%
      ggplot(aes(cellsize, hits)) +
      geom_line() +
      geom_point() +
      scale_x_log10() +
      xlab("Cell size") +
      ylab("Proportion of runs with correct estimate of soil types")


    gridExtra::grid.arrange(plot, plot2, plot3, nrow = 1)
  }


  ## save data ===================
  if (hpcc) {
    data_to_save <-
      list(
        parms = cbind(
          scen,
          numclusters = numclusters_cells,
          clustering.method = clustering.method,
          RV = RV,
          clustering.statistic = clustering.statistic
        ),
        Landscape = Landscape,
        Community = Community,
        Cmatrix = Cmatrix,
        est_C = est_C,
        Nij = Tabulation,
        est_Nij = est_Nij
      )
    save(
      data_to_save,
      file = paste0(
        "/data/rdandrea-dw9/SpatialNiche/Data/20191108/",
        formatC(simind, format = "d", width = nchar(nrow(scenarios)), flag = "0"),
        ".RData"
      )
    )
  }
}
#  },
#  silent=TRUE)
# }

if (calculate.likelihood) {
  tN <- as.matrix(Tabulation %>% select(-cell))
  kvec <- 2:20
  loglik <- sapply(kvec, function(n) {
    cluster_index <- kmeans(tN, centers = n, nstart = 1e3)$cluster
    est_Nij <-
      Tabulation %>%
      mutate(soiltype = cluster_index) %>%
      select(-cell) %>%
      group_by(soiltype) %>%
      summarize_all(sum) %>%
      filter(soiltype > 0)

    est_C <- Generate_estC(est_Nij, R)
    C_prob <- t(t(est_C) / colSums(est_C))
    loglik <- tr(tN %*% log(1e-16 + C_prob[, cluster_index]))
  })
}

if (simulate.likelihood) {
  run <- 0
  A <- 120^2
  NST <- 5
  S <- 30
  landscape_autocor_parm <- 50
  C_sd <- 0.5
  C_mean <- 1
  C_cv <- 1
  Landscape <- Generate_Landscape(seed = run, landscape_autocor_parm = landscape_autocor_parm)
  R <- as.numeric(table(Landscape))
  non_preferred_site_affinity <- 0.025
  Cmatrix <- Generate_Cmatrix(seed = run, mode = 'spec')
  # Cmatrix=diag(S)
  Community <- Generate_Community(Landscape, Cmatrix, seed = run)

  ## knowing both C and the Landscape
  N <- as.numeric(table(Community))
  P0ij <- t(t(N*Cmatrix)/colSums(N*Cmatrix))
  foo <- tibble(i=Community,j=Landscape)
  bar <- apply(foo,1,function(v) P0ij[v[1],v[2]])
  loglik0 <- sum(log(bar))
  
  res <- 
    ddply(
      data.frame(NC = unique(scenarios$NC)) %>% 
      filter(NC > S & NC <= A/9), .(NC), function(v) {
        NC <- v$NC
        print(NC)
        cellgrid <- Generate_Grid(A, NC)
        Tabulation <- Generate_Tabulation(cellgrid, Community)
        abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
        tN <- t(abunds_by_cells)
        tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
    
        ## estimate matrix Nij based on the identification of cells to soil types
        # dis_sites <- vegdist(Tabulation %>% select(-1),method="jaccard")
        # clustermodel=pamk(dis_sites,krange=2:min(S,NC-1),diss=TRUE)
        # clustermodel=pamk(Tabulation %>% select(-1),krange=2:min(S,NC-1),diss=FALSE)
        # numclusters_cells=clustermodel$nc
        # cluster_index=clustermodel$pamobject$clustering
        # clustering.statistic=as.numeric(clustermodel$pamobject$objective[2])
    
    
        # clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
        # cluster_index=clustermodel$cluster
        # clustering.statistic=clustermodel$tot.withinss
        # numclusters_cells=S
    
        dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
        clus_sites <- hclust(dis_sites)
        cluster_index <- cutree(clus_sites, NST)
        clustering.statistic <- NA
        numclusters_cells <- S
    
        est_Nij <-
          Tabulation %>%
          mutate(soiltype = cluster_index) %>%
          select(-cell) %>%
          group_by(soiltype) %>%
          summarize_all(sum) %>%
          filter(soiltype > 0)
    
        est_C <- Generate_estC(est_Nij, R)
        Pij <- t(t(N*est_C)/colSums(N*est_C))
        
        
        ## estimating both C and the landscape
        logprob=log(Pij[, cluster_index])
        loglik <- sum(tN*logprob*is.finite(logprob),na.rm=TRUE)
        
        f=tibble(i=Community,j=Landscape)
        f$p=apply(f,1,function(v) Pij[v[1],v[2]])
        loglikprime=sum(log(f$p))
    
        ## knowing the landscape, estimating C
        foo <- tibble(i = Community, j = Landscape)
        foo$P <- apply(foo, 1, function(v) Pij[v[1], v[2]])
        loglik2 <- sum(log(foo$P)*is.finite(log(foo$P)),na.rm=TRUE)
    
        ## knowing C, estimating the landscape
        logprob0=log(P0ij[, cluster_index])
        loglik3 <- sum(tN*logprob0*is.finite(logprob0),na.rm=TRUE)
    
        
        tr <- function(m) sum(diag(m))
        squarify <- function(m) m %*% t(m)
        center <- function(m) t(t(m) - colMeans(m))
        CC <- squarify(center(Cmatrix))
        eCeC <- squarify(center(est_C))
        RV <- tr(CC %*% eCeC) / sqrt(tr(CC %*% CC) * tr(eCeC %*% eCeC))
    
        return(
          data.frame(
            NC = NC,
            numclusters = numclusters_cells,
            loglik = loglik,
            loglik2 = loglik2,
            loglik3 = loglik3,
            statistic = clustering.statistic,
            RV = RV
          )
        )
      }
    ) %>%
    mutate(cellsize = sqrt(A / NC))
  
  {
    NC=with(res,NC[which.max(loglik)])
    cellgrid <- Generate_Grid(A, NC)
    Tabulation <- Generate_Tabulation(cellgrid, Community)
    abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
    tN <- t(abunds_by_cells)
    tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
    
    dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
    clus_sites <- hclust(dis_sites)
    cluster_index <- cutree(clus_sites, NST)
    clustering.statistic <- NA
    numclusters_cells <- S
    
    est_Nij <-
      Tabulation %>%
      mutate(soiltype = cluster_index) %>%
      select(-cell) %>%
      group_by(soiltype) %>%
      summarize_all(sum) %>%
      filter(soiltype > 0)
    
    est_C <- ReorderMatrix(Generate_estC(est_Nij, R))
    
  }

  cellsize=res$cellsize
  plot0 <-
    res %>%
    ggplot() +
    #scale_x_log10() +
    ylab('Log Likelihood') +
    xlab('Cell Linear Size') +
    scale_x_continuous(breaks=cellsize, labels=cellsize)

  plotL <- tibble(landscape=Landscape,community=Community) %>% cbind(expand.grid(x=1:120,y=1:120)) %>%
    ggplot(aes(x,y,fill=landscape)) + geom_tile() + theme(legend.position = 'none') +
    ggtitle('Landscape')
  
  plotC <- tibble(landscape=Landscape,community=Community) %>% cbind(expand.grid(x=1:120,y=1:120)) %>%
    ggplot(aes(x,y,fill=community)) + geom_tile()  + theme(legend.position = 'none') +
    ggtitle('Community')
  
  plotCmatrix <- tibble(C=as.numeric(Cmatrix)) %>% cbind(expand.grid(species=seq(S),soiltype=seq(NST))) %>%
    ggplot(aes(x=species,y=soiltype,fill=C)) + geom_tile()  + theme(legend.position = 'none') +
    ggtitle('C')
  
  plotCest <- tibble(C=as.numeric(est_C)) %>% cbind(expand.grid(species=seq(S),soiltype=seq(NST))) %>%
    ggplot(aes(x=species,y=soiltype,fill=C)) + geom_tile()  + theme(legend.position = 'none') +
    ggtitle('estimated C')

  plot1 <- plot0 + geom_line(aes(cellsize, loglik/1e3)) + geom_point(aes(cellsize, loglik/1e3)) +
    ggtitle('Estimating both C and the Landscape') +
    geom_hline(yintercept=loglik0/1e3,color=red)

  plot2 <- plot0 + geom_line(aes(cellsize, loglik2)) + geom_point(aes(cellsize, loglik2)) +
    ggtitle('Knowing the Landscape, estimating C') +
    geom_hline(yintercept=loglik0,color=red)

  plot3 <- plot0 + geom_line(aes(cellsize, loglik3)) + geom_point(aes(cellsize, loglik3)) +
    ggtitle('Knowing C, estimating the Landscape') +
    geom_hline(yintercept=loglik0,color=red)

  plotstat <- plot0 + geom_line(aes(cellsize, statistic)) + geom_point(aes(cellsize, statistic))

  plotRV <- plot0 + geom_line(aes(cellsize, RV)) + geom_point(aes(cellsize, RV)) +
    ggtitle('Match between true C and estimated C') +
    ylab('RV') +
    ylim(c(0,1))

  # gridExtra::grid.arrange(plotL, plot1, plot2, plotC, plot3, plotRV, nrow = 2)
  gridExtra::grid.arrange(plotL, plotCmatrix, plot1, plotC, plotCest, plotRV, nrow = 2)
}

if (simulate.likelihood2) {
  clustering.mode <- 'hierarchical'
  
  run <- 10
  A <- 120^2
  NST <- 20
  S <- 20
  landscape_autocor_parm <- 50
  C_sd <- 0.5
  C_mean <- 1
  C_cv <- 1
  Landscape <- Generate_Landscape(seed = run, landscape_autocor_parm = landscape_autocor_parm)
  R <- as.numeric(table(Landscape))
  non_preferred_site_affinity <- 0.025
  Cmatrix <- Generate_Cmatrix(seed = run, mode = 'spec')
  # Cmatrix=diag(S)
  Community <- Generate_Community(Landscape, Cmatrix, seed = run)
  
  ## knowing both C and the Landscape
  N <- as.numeric(table(Community))
  P0ij <- t(t(N*Cmatrix)/colSums(N*Cmatrix))
  foo <- tibble(i=Community,j=Landscape)
  bar <- apply(foo,1,function(v) P0ij[v[1],v[2]])
  loglik0 <- sum(log(bar))
  
  res <- {
    ddply(
      crossing(tibble(NC = unique(scenarios$NC)),tibble(NST=3:30)) %>% 
        filter(NC > S & NC <= A/9), .(NST,NC), function(v) {
          NST <- v$NST
          NC <- v$NC
          print(paste(NST,NC))
          cellgrid <- Generate_Grid(A, NC)
          Tabulation <- Generate_Tabulation(cellgrid, Community)
          abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
          tN <- t(abunds_by_cells)
          tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
          
          ## estimate matrix Nij based on the identification of cells to soil types
          
          ## Using k-means
          if(clustering.mode=='kmeans'){
            clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
            cluster_index=clustermodel$cluster
            clustering.statistic=clustermodel$tot.withinss
            numclusters_cells <- NST
          }
          
          ## Using hierarchical clustering
          if(clustering.mode=='hierarchical'){
            dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
            clus_sites <- hclust(dis_sites)
            cluster_index <- cutree(clus_sites, NST)
            clustering.statistic <- NA
            numclusters_cells <- NST
          }
          
          ## Using pam clustering
          if(clustering.mode=='pam'){
            dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
            clustermodel=cluster::pam(x=dis_sites, k=NST)
            cluster_index <- clustermodel$clustering
            clustering.statistic <- as.numeric(clustermodel$objective['swap'])
            numclusters_cells <- NST
          }
            
          est_Nij <-
            Tabulation %>%
            mutate(soiltype = cluster_index) %>%
            select(-cell) %>%
            group_by(soiltype) %>%
            summarize_all(sum) %>%
            filter(soiltype > 0)
          
          est_C <- Generate_estC(est_Nij, R)
          Pij <- t(t(N*est_C)/colSums(N*est_C))
          
          
          logprob=log(Pij[, cluster_index])
          loglik <- sum(tN*logprob*is.finite(logprob),na.rm=TRUE)
          
          tr <- function(m) sum(diag(m))
          squarify <- function(m) m %*% t(m)
          center <- function(m) t(t(m) - colMeans(m))
          CC <- squarify(center(Cmatrix))
          eCeC <- squarify(center(est_C))
          RV <- tr(CC %*% eCeC) / sqrt(tr(CC %*% CC) * tr(eCeC %*% eCeC))
          
          return(
            data.frame(
              NC = NC,
              numclusters = numclusters_cells,
              loglik = loglik,
              statistic = clustering.statistic,
              RV = RV
            )
          )
        }
    ) %>%
    mutate(cellsize = sqrt(A / NC)) %>%
    as_tibble
  }
  
  {
    df_maxloglik <-
      res %>% 
      group_by(NST) %>%
      summarize(
        bestcellsize=cellsize[which.max(loglik)],
        maxloglik=loglik[which.max(loglik)],
        maxRV=RV[which.max(loglik)]
      )
    nst_continuum <- with(df_maxloglik,seq(min(NST),max(NST),l=1e3))
    smoothed_maxloglik <- with(df_maxloglik, predict(loess(maxloglik~NST,span=.3),nst_continuum))
    index <- Max2Deriv(nst_continuum, smoothed_maxloglik)
    bestNST <- unique(round(nst_continuum[index]))
    bestcellsize <- df_maxloglik %>% filter(NST==bestNST) %>% pull(bestcellsize)
    
    nst=df_maxloglik$NST
    
    NC=A/bestcellsize^2
    NST <- bestNST
    cellgrid <- Generate_Grid(A, NC)
    Tabulation <- Generate_Tabulation(cellgrid, Community)
    abunds_by_cells <- as.matrix(Tabulation %>% select(-cell))
    tN <- t(abunds_by_cells)
    tN_normalized <- t(abunds_by_cells) / colSums(abunds_by_cells) ## species x cells
    
    ## Using k-means
    if(clustering.mode=='kmeans'){
      clustermodel=kmeans(tN,centers=S,iter=100,nstart=1e3)
      cluster_index=clustermodel$cluster
      clustering.statistic=clustermodel$tot.withinss
      numclusters_cells <- NST
    }
    
    ## Using hierarchical clustering
    if(clustering.mode=='hierarchical'){
      dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
      clus_sites <- hclust(dis_sites)
      cluster_index <- cutree(clus_sites, NST)
      clustering.statistic <- NA
      numclusters_cells <- S
    }
    
    ## Using pam clustering
    if(clustering.mode=='pam'){
      dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
      clustermodel=cluster::pam(x=dis_sites, k=NST)
      cluster_index <- clustermodel$clustering
      clustering.statistic <- as.numeric(clustermodel$objective['swap'])
      numclusters_cells <- NST
    }
    
    est_Nij <-
      Tabulation %>%
      mutate(soiltype = cluster_index) %>%
      select(-cell) %>%
      group_by(soiltype) %>%
      summarize_all(sum) %>%
      filter(soiltype > 0)
    
    est_C <- ReorderMatrix(Generate_estC(est_Nij, R))
    
    }
  
  cellsize=unique(res$cellsize)
  plot0 <-
    res %>%
    ggplot() +
    theme_bw() +
    #scale_x_log10() +
    ylab('Log Likelihood') +
    xlab('Cell Linear Size') +
    scale_x_continuous(breaks=cellsize, labels=cellsize, minor_breaks = NULL)
  
  plot_Landscape <- 
    tibble(landscape=Landscape,community=Community) %>% 
    cbind(expand.grid(x=1:120,y=1:120)) %>%
    ggplot(aes(x,y,fill=landscape)) + 
    theme_bw() +
    geom_tile() + 
    theme(legend.position = 'none') +
    ggtitle('Landscape')  +
    scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
    scale_y_continuous(breaks=NULL,minor_breaks = NULL)
  
  plot_Community <- 
    tibble(landscape=Landscape,community=Community) %>% 
    cbind(expand.grid(x=1:120,y=1:120)) %>%
    ggplot(aes(x,y,fill=community)) + 
    theme_bw() +
    geom_tile()  + 
    theme(legend.position = 'none') +
    ggtitle('Community') +
    scale_x_continuous(breaks=NULL,minor_breaks = NULL) +
    scale_y_continuous(breaks=NULL,minor_breaks = NULL)
  
  plot_Cmatrix <- 
    tibble(C=as.numeric(Cmatrix)) %>% 
    cbind(expand.grid(species=seq(nrow(Cmatrix)),soiltype=seq(ncol(Cmatrix)))) %>%
    ggplot(aes(x=species,y=soiltype,fill=C)) + 
    theme_bw() +
    geom_tile()  + 
    theme(legend.position = 'none') +
    ggtitle('C matrix') 
  
  plot_Cest <- 
    tibble(C=as.numeric(est_C)) %>% 
    cbind(expand.grid(species=seq(S),soiltype=seq(NST))) %>%
    ggplot(aes(x=species,y=soiltype,fill=C)) + 
    theme_bw() +
    geom_tile()  + 
    theme(legend.position = 'none') +
    ggtitle('estimated C')
  
  plot_cellsize_by_loglik <- 
    plot0 + 
    geom_line(aes(cellsize, loglik/1e3,group=as.factor(NST),color=as.factor(NST))) + 
    geom_point(aes(cellsize, loglik/1e3,group=as.factor(NST),color=as.factor(NST))) +
    ggtitle('Estimating both C and the Landscape') +
    geom_hline(yintercept=loglik0/1e3,color=red) +
    theme(legend.position='none')
  
  plot_cellsize_by_RV <- 
    plot0 + 
    geom_line(aes(cellsize, RV,group=as.factor(NST),color=as.factor(NST))) + 
    geom_point(aes(cellsize, RV,group=as.factor(NST),color=as.factor(NST))) +
    ggtitle('Match between true C and estimated C') +
    ylab('RV') +
    ylim(c(0,1)) +
    theme(legend.position='none')
  
  plot_NST_by_maxloglik <-
    df_maxloglik %>%
    ggplot(aes(NST,maxloglik/1e3)) +
    geom_vline(aes(xintercept = bestNST),color=red) +
    geom_line(aes(x, y/1e3),data=tibble(x=nst_continuum,y=smoothed_maxloglik)) +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=nst, labels=nst, minor_breaks = NULL) +
    xlab('Number of soil types') +
    ylab('Max log likelihood (/ 1e3)')
  
  plot_NST_by_maxRV <-
    df_maxloglik %>%
    ggplot(aes(NST,maxRV)) +
    geom_vline(aes(xintercept = bestNST),color=red) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=nst, labels=nst, minor_breaks = NULL) +
    xlab('Number of soil types') +
    ylab('Max RV') +
    ylim(c(0,1))
  
  plot_NST_by_bestcellsize <-
    df_maxloglik %>%
    ggplot(aes(NST,bestcellsize)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks=nst, labels=nst, minor_breaks = NULL) +
    scale_y_continuous(breaks=cellsize, labels=cellsize, minor_breaks = NULL) +
    xlab('Number of soil types') +
    ylab('Best cell size') +
    ylim(range(cellsize))
  
  plotstat <- 
    plot0 + 
    geom_line(aes(cellsize, statistic)) + 
    geom_point(aes(cellsize, statistic))
  
  plotRV <- 
    plot0 + 
    geom_line(aes(cellsize, RV)) + 
    geom_point(aes(cellsize, RV)) +
    ggtitle('Match between true C and estimated C') +
    ylab('RV') +
    ylim(c(0,1))
  
  gridExtra::grid.arrange(
    plot_Landscape, 
    plot_Cmatrix, 
    plot_Community,
    plot_cellsize_by_loglik,
    plot_cellsize_by_RV,
    plot_NST_by_maxloglik,
    plot_NST_by_maxRV,
    plot_NST_by_bestcellsize,
    plot_Cest,
    nrow=3
  )
  
  # gridExtra::grid.arrange(plotL, plot1, plot2, plotC, plot3, plotRV, nrow = 2)
  # gridExtra::grid.arrange(plotL, plotCmatrix, plot1, plotC, plotCest, plotRV, nrow = 2)
}



## analyze data ===================
if (!hpcc & analyze.data) {
  datadir <- switch(
    match(user, c("rafael", "wangdz")),
    "~/SpatialNiche/Data/20191029/",
    "/Users/wangdz/Documents/UIUC/Odywer/SpatialNiche/Data/20191029/"
  )
  files <- formatC(seq(nrow(scenarios)), width = nchar(nrow(scenarios)), flag = "0")
  dtf <- NULL
  for (char in intersect(paste0(files, ".RData"), list.files(datadir))) dtf <- rbind(dtf, get(load(paste0(datadir, char)))$parms)
  dtf <- as_tibble(dtf)

  dtf %>%
    mutate(S = factor(S)) %>%
    ggplot(aes(NC, numclusters, group = S, color = S)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_continuous(name = "Estimated number of clusters", breaks = seq(1, 12, by = 2), labels = waiver()) +
    xlab("Number of cells") +
    facet_grid(landscape_autocor_parm ~ C_sd)

  # files=formatC(seq(336),width=nchar(nrow(scenarios)),flag='0')
  # for(char in intersect(paste0(files,'.RData'),list.files(datadir))){
  #   dtf=get(load(paste0(datadir,char)))
  #   dtf$parms$mode='spec'
  #   dtf$parms$C_mean=NA
  #   dtf$parms$C_cv=NA
  #   save(dtf,file=paste0(datadir,char))
  # }
}

if (!hpcc & analyze.data2) {
  library(tidyr)
  library(dplyr)
  data <- NULL
  for (char in list.files()) data <- rbind(data, get(load(char))$parms)
  data <- as_tibble(data)
  save(data, file = "/data/rdandrea-dw9/SpatialNiche/Data/SumDataPAM.RData")


  setwd("~/spatialniche/data/20191108/")
  data <- get(load("SumDataPAM.RData"))

  plot_RV <- {
    data %>%
      filter(mode == "spec") %>%
      group_by(A, S, NST, NC, C_sd, landscape_autocor_parm) %>%
      summarize_at(
        c("RV", "numclusters", "clustering.statistic"),
        .funs = list(mean = mean, sd = sd), na.rm = TRUE
      ) %>%
      ungroup() %>%
      mutate(cellsize = sqrt(A / NC), S = factor(S), NST = factor(NST)) %>%
      ggplot(aes(cellsize, RV_mean, group = S, color = S)) +
      geom_ribbon(aes(x = cellsize, ymin = RV_mean - RV_sd, ymax = RV_mean + RV_sd), alpha = .2, color = NA) +
      geom_line() +
      geom_point() +
      scale_x_log10() +
      facet_wrap(landscape_autocor_parm ~ C_sd) +
      ggtitle("Specialists")
  }

  plot_nclus <- {
    data %>%
      filter(mode == "spec") %>%
      group_by(S, NST, NC, C_sd, landscape_autocor_parm) %>%
      summarize_at(
        c("RV", "numclusters", "clustering.statistic"),
        .funs = list(mean = mean, sd = sd), na.rm = TRUE
      ) %>%
      ungroup() %>%
      mutate(cellsize = sqrt(A / NC), S = factor(S), NST = factor(NST)) %>%
      ggplot(aes(cellsize, numclusters_mean, group = NST, color = NST)) +
      geom_ribbon(
        aes(
          x = cellsize,
          ymin = numclusters_mean - numclusters_sd,
          ymax = numclusters_mean + numclusters_sd
        ),
        alpha = .2, color = NA
      ) +
      geom_line() +
      geom_point() +
      scale_x_log10() +
      facet_wrap(landscape_autocor_parm ~ C_sd) +
      ggtitle("Specialists")
  }

  plot_nclus2 <- {
    data %>%
      filter(mode == "spec" & !is.na(numclusters)) %>%
      mutate(hit = (numclusters == NST)) %>%
      group_by(A, S, NST, NC, C_sd, landscape_autocor_parm) %>%
      summarize(hits = sum(hit, na.rm = TRUE) / n()) %>%
      ungroup() %>%
      mutate(cellsize = sqrt(A / NC), S = factor(S), NST = factor(NST)) %>%
      ggplot(aes(cellsize, hits, group = NST, color = NST)) +
      geom_line() +
      geom_point() +
      scale_x_log10() +
      facet_wrap(landscape_autocor_parm ~ C_sd) +
      ggtitle("Specialists")
  }
}

if (FALSE) {
  dtf %>%
    mutate(
      cell_size = sqrt(A / NC),
      envt_range = c(21, 30, 71)[factor(landscape_autocor_parm)]
    ) %>%
    mutate(S = factor(S)) %>%
    ggplot(aes(cell_size, numclusters, group = S, color = S)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_continuous(name = "Estimated number of clusters", breaks = seq(1, 12, by = 2), labels = waiver()) +
    xlab("Cell linear size") +
    geom_vline(aes(xintercept = envt_range / 4), color = "gray") +
    facet_grid(envt_range ~ C_sd)
}

