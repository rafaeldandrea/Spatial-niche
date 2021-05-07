
## Load functions
#library(dplyr)
if(TRUE){
library(tidyverse)
library(vegan)
library(ggplot2)
library(wavethresh) ## for function guyrot() -- rotates a vector v by n steps
library(fpc) ## for clustering analysis
library(numbers) ## for function divisors()
library(MASS) ## for generalized inverse (left inverse) of C when S > NST
library(Rtsne)
select = dplyr::select
}
##
ReorderMatrix=function(data){
  mylist=NULL
  order_matrix=matrix(NA,nrow(data),ncol(data))
  for(k in 1:ncol(data)) order_matrix[k,]=order(data[,k],decreasing=TRUE)
  for(i in 1:ncol(data)){
    for(j in 1:nrow(data)){
      if(!order_matrix[i,j] %in% mylist){
        mylist=c(mylist,order_matrix[i,j])
        break
      }
    }
  }
  data=data[mylist,]
}
as_matrix <- function(x){
  if(!tibble::is_tibble(x) ) stop("x must be a tibble")
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}

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
if(!hpcc){
wd ='/Users/wangdz/Documents/UIUC/Odywer/BCI/' 
}else{
  wd = '/home/dw9/SpatialNiche/BCI'
}
setwd(wd)
filenames = gtools::mixedsort(list.files(pattern = 'bci*'))
scenarios <- {
  crossing(ncx = seq(50,300,by=50),
           ncy = seq(50,300,by=50),
           norm_distance = c(TRUE,FALSE))
}
scenarios = add_row(scenarios, ncx = 10,ncy = 10, norm_distance = TRUE)

theind = 73
ind <- make_index(theind)
## Choose row of 'parameters'data frame to specify scenario
scen = scenarios[ind, ]

## Assign each column of scen to a variable matching its name
list2env(as.list(scen), envir = globalenv())

NC = ncx*ncy
threshold = 0
cell_by_species = matrix(0, NC, 300)
for (char in filenames){
  dat = get(load(char))
  bci= filter(bci.full1,bci.full1$dbh>100&!is.na(bci.full1$gx)&!is.na(bci.full1$gy))
  #seeting a threshhold
  bci_amounts= filter(bci %>% group_by(bci$sp) %>% summarize(count=n()),count>threshold)
  #griding 
  gxmin=range(bci$gx)[1]
  gxmax=range(bci$gx)[2]
  gxbreaks=seq(gxmin,gxmax,length=ncx+1)
  cellindex_x=findInterval(bci$gx,gxbreaks,rightmost.closed=TRUE)
  gymin=range(bci$gy)[1]
  gymax=range(bci$gy)[2]
  gybreaks=seq(gymin,gymax,length=ncy+1)
  cellindex_y=findInterval(bci$gy,gybreaks,rightmost.closed=TRUE)
  bci=bci %>%
    mutate(
      cx=cellindex_x,
      cy=cellindex_y,
      cellindex=cx+(cy-1)*ncx
    )
  Tabulation <-
    tibble(
      cell=bci$cellindex,
      # soiltype=Landscape,
      species=bci$sp
    ) %>%
    group_by(cell) %>%
    count(species) %>%
    spread(key=species,value=n,sep='_') %>%
    ungroup
  Tabulation[is.na(Tabulation)]=0
  colnames(Tabulation)= unlist(strsplit(colnames(Tabulation),split='species_'))[c(TRUE,FALSE)]
  Tabulation_selected = Tabulation[bci_amounts$`bci$sp`]
  tabula = Tabulation_selected
  names(tabula) = paste0('species_', 1:ncol(tabula))
  presences = match(colnames(tabula), paste0('species_', 1:ncol(tabula)))
  absences = setdiff(1:ncol(tabula), presences)
  if(length(absences) > 0){
    tabula = cbind(tabula, matrix(0, NC, length(absences)))
    colnames(tabula) = c('cell', paste0('species_', presences), paste0('species_', absences))
  }
  tabula = select(tabula, paste0('species_', 1:ncol(tabula)))
  cell_by_species = cell_by_species + tabula
} 

#clustermodel=pamk(dis_sites,krange=2:min(S,ncx*ncy-1),diss=TRUE)
do.nmds = FALSE
if(do.nmds){
  #rownames(matrix_Tabulation) <-paste("c",1:25,sep="")
  b = character(0)
  for (it in 2:2){
    NMDS = metaMDS(Tabulation_selected,k=it,trymax=200)
    jpeg(paste(it,".jpg",sep=""))
    plot(NMDS,type='t',main=paste(it,"D",sep=""))
    dev.off()
    print(NMDS$stress)
    b <- c(b,NMDS$stress)
  }
 # plot(2:20,b)
}
## pam clustering analysis
do.pamk = FALSE
if(do.pamk){
  ## estimate matrix Nij based on the identification of cells to soil types
  dis_sites <- vegdist(Tabulation_selected ,method="bray")
  clustermodel=pamk(dis_sites,krange=2:30,diss=TRUE)
  #clustermodel=pamk(Tabulation %>% select(-1),krange=2:min(S,NC-1),diss=FALSE)
  
  numclusters_cells=clustermodel$nc
  est_Nij <-
    Tabulation_selected %>% 
    mutate(soiltype=clustermodel$pamobject$clustering) %>%
    group_by(soiltype) %>% 
    summarize_all(sum)
}
do.hierarch=TRUE
if(do.hierarch){
  linkage = 'ward.D2'
  tN <- t(cell_by_species)
  tN_normalized <- tN / rowSums(tN) ## species x cells
  if(norm_distance){
    dis_sites <- vegdist(t(tN_normalized), method = "jaccard")
  }else{
    dis_sites <- vegdist(t(tN), method = "jaccard")
  }
  #clus_sites <- hclust(dis_sites, method = linkage)
  #dis_sites <- vegdist(t(tN), method = "jaccard")
  
  clustresult <- 
    NbClust::NbClust(
      data = cell_by_species, 
      diss = dis_sites, 
      distance = NULL, 
      min.nc = 20, 
      max.nc = min(300, ncx*ncy - 1),
      method = linkage, 
      index = 'sdindex'
    )
  result = c(clustresult,ncx=ncx,ncy=ncy,norm = norm_distance)
  if (hpcc){
    datadir <- paste0('~/SpatialNiche/Data/', Sys.Date(),'-BCIanalysis', '/')
    dir.create(datadir, showWarnings = FALSE)
    savename = paste0(datadir,'BCI','-ncx-',ncx,'-ncy-',ncy,'-norm-',norm_distance,'.RData')
    save(result, file = savename)
  } else {
    plot(20: 400, clustresult$All.index, t='l',main = paste0("threshold=",threshold,',NC=',NC))
}
}
if(FALSE){
## estimate C matrix
est_C=t(matrix(as.numeric(unlist(est_Nij[,-1])),nrow(est_Nij),ncol(est_Nij)-1))
est_C=est_C/rowSums(est_C)
#est_C=round(est_C/max(est_C),1)
est_C=ReorderMatrix(est_C)
#plot est_c
plot_est_C <- 
  as_tibble(est_C) %>% 
  `colnames<-`(formatC(1:ncol(est_C),width=3,flag="0")) %>% 
  mutate(species=formatC(1:nrow(est_C),width=3,flag="0")) %>% 
  gather('soiltype','value',-species) %>% 
  ggplot(aes(species,soiltype,fill=value)) + 
  geom_raster() +
 # scale_fill_gradient(low=colmin,high=colmax) +
  theme(aspect.ratio = 1) +
  theme(legend.position = 1) +
  scale_y_discrete(
    # breaks=formatC(round(seq(1,numclusters_cells,length=5)),width=3,flag="0"),
    labels=as.character(1:ncol(est_C))
  ) +
  scale_x_discrete(
    #breaks=formatC(round(seq(1,S,length=5)),width=3,flag="0"),
    labels=as.character(1:nrow(est_C))
  )  +
  ggtitle('Estimated C matrix')
plot_est_C
plot(NMDS,type='t',main=paste(it,"D",sep=""))
ordihull(NMDS, grp_sites, lty = 2, col = "red")
}
                                                                                           