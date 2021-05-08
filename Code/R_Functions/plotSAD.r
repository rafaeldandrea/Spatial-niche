plotSAD=function(data,dbn){
  n=if(is.numeric(data)) plyr::count(data) else data 
  mod=fitSAD(n,dbn=dbn)
  p0 <-
    n %>% 
    mutate(cum=(sum(freq)-cumsum(freq)+freq)/sum(freq)) %>% 
    ggplot(aes(x,cum)) + 
    geom_line(aes(x,cum),data=mod$dtf,color=red) + 
    geom_point() +
    xlim(range(n$x)) +
    theme_bw() +
    theme(aspect.ratio=1) +
    annotate(x=.9*max(n$x),y=.8,'text',label=paste('p-value = ',round(mod$p.value,2))) +
	ylab('Proportion of species') +
	xlab('Abundance')
  
  return(p0)
}