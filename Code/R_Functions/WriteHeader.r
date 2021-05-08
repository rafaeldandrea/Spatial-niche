## Saves dtf as a file under the name filename, then attaches header to top and saves again.
WriteHeader=function(dtf,filename,header=NULL,append=FALSE){
	if(!append | !file.exists(filename)){
		write.table(dtf,filename,quote=FALSE,row.names=FALSE)
		header=paste0('## Created on ',date(),'\n',header,'\n')
		txt=c(header,readLines(filename))
		write(txt,filename)
	}else write.table(dtf,filename,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
}