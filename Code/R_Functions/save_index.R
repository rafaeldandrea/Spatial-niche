save_index=function(data,index,scenarios,folder){
  save(data,file=paste0(folder,formatC(index,width=nchar(nrow(scenarios)),flag="0"),'.RData'))
}