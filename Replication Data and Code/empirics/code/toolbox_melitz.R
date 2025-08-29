#' @export

checkdim = function(x,x_name,N,S=NULL,destinationspecific=F){
  
  
  if(is.null(S)){
    
    if(is.vector(x) & length(as.vector(x))==1){
      x = rep(x,length=N)
    }
    
  }else if(destinationspecific){
    
    if(is.vector(x) & length(as.vector(x))==1){
      x = array(x,dim=c(N,N,S))
    }
    else if(dim(x)!=c(N,N,S)){
      stop(paste0("Wrong input dimension, ",x_name," should be a uni-dimensional vector or an array of dimension (",N,",",N,",",S,")"))
    }
    
  }else{
    
    if(is.vector(x) & length(as.vector(x))==1){
      x = array(x,dim=c(N,1,S))
    }
    else if(dim(x)!=c(N,1,S)){
      stop(paste0("Wrong input dimension, ",x_name," should be a uni-dimensional vector or an array of dimension (",N,",",1,",",S,")"))
    }
    
  }
  
  return(x)
}

#' @export

checkdiag = function(x){
  
  int = dim(x)
  x = lapply(int[3],FUN = function(i,data){
    diag(data[,,i]) = 1 
    return(data[,,i])
  },data=x)
  x=array(as.numeric(unlist(x)), dim=int)
  return(x)
  
}

#' export

oppdiag = function(x){
  return(x[(row(x) - col(x))!=0])
}

