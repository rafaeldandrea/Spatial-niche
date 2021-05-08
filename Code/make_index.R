## Read scenario from either manual input or command line
make_index <- function(manual = NULL) {
  
  read_index <- as.numeric(commandArgs(TRUE))[1]
  
  if (is.na(read_index) & is.null(manual)) {
  
    stop('R is not being called from command line and manual index is not provided.')
  
  }
  
  return(ifelse(!is.na(read_index), read_index, manual))

}
