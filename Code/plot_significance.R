library(dplyr)
library(ggplot2)
datadir='/Users/wangdzh/Downloads/20200219'
setwd(datadir)
lf=setdiff(list.files(datadir),c('data.RData','readme.txt'))
dat=NULL
for(char in lf){
  foo=get(load(char))
  dat=rbind(dat,cbind(foo$scen,testresult=foo$testresult))
}
dat = dat %>%
  mutate(significant=1*(testresult>.05)) %>%
  as_tibble
save(dat,file='AllResults.RData')             
results <-
  dat %>%
  group_by(clustering.mode,C_sd) %>%
  summarize(percent_significant=sum(significant)/length(significant))
red <- "#DD5144"
plot <-
  ggplot() +
  geom_point(aes(C_sd,percent_significant,color=red),data=results[27:52,]) +
  geom_point(aes(C_sd,percent_significant),data=results[1:26,]) +
  labs(
    x='Std deviation in C matrix',
    y='Percent of runs with significant difference between data and inferemce')+
  theme_bw()
