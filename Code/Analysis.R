## Code for reading and analysing data

## set data folder
datadir='~/spatialniche/data/20191001/'

## read all data files
lf=setdiff(list.files(datadir),c('Parameters.RData','readme.txt'))

scenario=snapshot=pval=stat=NULL

for(char in lf){
  
  ## read the corresponding data file
  dat=get(load(paste0(datadir,char)))
  
  ## select the abundance matrix from the data file
  N=dat$N
  
  ## fit the log-series
  lsfit=fitSAD(data=rowSums(N),dbn='LS',method='CVM')
  
  ## select the data frame containing the fitted values
  fit=lsfit$dtf %>%
    as_tibble
  
  ## process the observed abundances into a data frame with cumulative counts
  data=plyr::count(rowSums(N)) %>%
    mutate(
      cdf=cumsum(freq),
      cum=(sum(freq)-cdf+freq)/sum(freq)
    ) %>%
    as_tibble
  
  ## merge the observed and fitted data frames
  dtf <-
    data %>%
    inner_join(fit,by='x')
  
  ## plot the merged results as points (data) and red curve (fitted LS)
  myplot <- 
    dtf %>% 
    ggplot() + 
    geom_line(aes(x,cum.y),color=red) +
    geom_point(aes(x,cum.x)) + 
    scale_x_log10() + 
    scale_y_log10() +
    theme_bw() +
    #theme(aspect.ratio = 1) +
    xlab('Abundance') +
    ylab('Cumulative probability')+
    ggtitle(paste0('data file = ',char,'  p-value of LS fit = ',round(lsfit$p.value,2)))
  
  ## assign a name to the plot variable based on the data file used
  assign(paste0('plot',match(char,lf)),myplot)
  
  ## get p-value and test statistic from the LS fit
  foo=as.numeric(unlist(strsplit(unlist(strsplit(char,split='.RData')),split='_')))
  scenario=c(scenario,foo[1])
  snapshot=c(snapshot,foo[2])
  pval=c(pval,lsfit$p.value)
  stat=c(stat,lsfit$statistic)
}

dtf <-
  tibble(
    scenario=scenario,
    snapshot=snapshot,
    stat=stat,
    pval=pval
  )

## show the plots on a grid
gridExtra::grid.arrange(plot1,plot2,plot3,plot4,plot5,nrow=3)
