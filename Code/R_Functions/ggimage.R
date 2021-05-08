ggimage=function(C,col.min='black',col.max='red'){
  p1<-
    reshape2::melt(C) %>% 
    as_tibble %>% 
    ggplot(aes(x=Var2,y=Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient(low=col.min,high=col.max) + 
    theme_bw() +
    xlab('Column') +
    ylab('Row')
  return(p1)
}