oppdiag <- function(x){
  
  y = rev(x[(row(x) + col(x)) == (nrow(x) + 1)])
  return(y)
}


test_boot <- function(test = "nmad_test", x = NULL, n = NULL, dist, prior, coeff, B = 99, est = T, quant = NULL, lowertrunc=NULL, uppertrunc=NULL,  ...) {
  
  # Calculate original test statistic
  if (is.null(x)) {
    if (is.null(n)) {
      stop("n should be provided")
    }
    boottype <- 1
    stat <- do.call(test, c(..., list(dist = dist, prior = prior, coeff = coeff , lowertrunc=lowertrunc, uppertrunc=uppertrunc)))
  } else {
    boottype <- 2
    n <- length(x)
    stat <- do.call(test, c(..., list(x = x, dist = dist, prior = prior, coeff = coeff , lowertrunc=lowertrunc, uppertrunc=uppertrunc)))
  }

  stat <- as.numeric(stat)

  stat_star <- sapply(1:B, function(i, n, dist, prior, coeff, est, boottype, ...) {
    
    
    x_star <- rcombdist(n = n, dist = dist, prior = prior, coeff = coeff, lowertrunc=lowertrunc, uppertrunc=uppertrunc )
    x_star <- sort(x_star)

    if (est) {
      coeff.out = combdist.mle(dist=dist,x=x_star,start=pars$coefficients,steps=1)$pars[[1]]
      prior_star
      coeff_star
      extra_args[which(c("prior","coeff") %in% names(extra_args))] = list(coeff = coeff_star, prior = prior_star)[which(c("prior","coeff") %in% names(extra_args))]
    } else {
      prior_star <- prior
      coeff_star <- coeff
    }
    

    if (boottype == 1) {
      dist_star <- "empirical"
      prior_star <- 1
      coeff_star <- list(data = x_star)
      stat_star <- do.call(test, c(..., list(dist = dist_star, prior = prior_star, coeff = coeff_star, lowertrunc=lowertrunc, uppertrunc=uppertrunc)))
    } else {
      dist_star <- dist
      stat_star <- do.call(test, c(..., list(x = x_star, dist = dist_star, prior = prior_star, coeff = coeff_star, lowertrunc=lowertrunc, uppertrunc=uppertrunc)))
    }

    return(as.numeric(stat_star))
  }, n = n, dist = dist, prior = prior, coeff = coeff, est = est, boottype = boottype, ...)

  stat_star <- matrix(stat_star, nrow = length(stat), ncol = B)

  p.value <- (rowSums(stat_star >= stat) + 1) / (B + 1)

  if (is.null(quant)) {
    return_list <- cbind.data.frame(p.value = p.value, stat = stat, stat_star = stat_star)
  } else {
    conf.int <- t(apply(cbind(stat, stat_star), 1, function(x) {
      quantile(x = x, quant, na.rm = T)
    }))
    return_list <- cbind.data.frame(p.value = p.value, stat = stat, conf.int = conf.int)
  }

  return(return_list)
}


coeff.out = function(dist,prior = NULL,coeff,nested=F){
  

  dist <- unlist(strsplit(dist, split = "_"))
  
  if (length(dist) > 1) {
    
    # Composite distribution
    
    names_int = names(coeff)
    
    names_int=str_replace_all(string=names_int,pattern="meanlog", replacement="\\\\mu")
    names_int=str_replace_all(string=names_int,pattern="sdlog", replacement="SD")
    names_int=str_replace_all(string=names_int,pattern="shape", replacement="k")
    names_int=str_replace_all(string=names_int,pattern="scale", replacement="s")
    
    int = grep("coeff1.",names_int)
    names_int[int] = paste0(str_replace_all(string=names_int[int],pattern="coeff1.", replacement=""),"_1")
    int = grep("coeff2.",names_int)
    names_int[int] = paste0(str_replace_all(string=names_int[int],pattern="coeff2.", replacement=""),"_2")
    int = grep("coeff3.",names_int)
    if(length(int)>0){
      names_int[int] = paste0(str_replace_all(string=names_int[int],pattern="coeff3.", replacement=""),"_3")
    }
    
    prior.out = paste0("$\\pi_",1,"$=",formatC(prior,digits=2,format="f"))
    coeff.out = paste0("$",names_int,"$=",formatC(coeff,digits=2,format="f"),collapse=", ")
    
  }else if(!is.matrix(coeff)){
    
    #Single distributions
    
    if(dist=="doubleparetolognormal"){
      coeff = coeff[c("shape1","meanlog", "sdlog","shape2")]
    }else if(dist=="leftparetolognormal"){
      coeff = coeff[c("shape1","meanlog", "sdlog")]
    } else if(dist=="leftparetolognormal"){
      coeff = coeff[c("shape1","meanlog", "sdlog")]
    } else if(dist == "loglogis") {
      log_likelihood <- function(params, data) {
        shape <- params[1]
        scale <- params[2]
        -sum(dloglogis(data, shape, scale, log = TRUE))  # Negative log-likelihood
      }
    }
    
    
    names_int = names(coeff)
    
    names_int=str_replace_all(string=names_int,pattern="meanlog", replacement="\\\\mu")
    names_int=str_replace_all(string=names_int,pattern="sdlog", replacement="SD")
    names_int=str_replace_all(string=names_int,pattern="xmin", replacement="x_{min}")
    names_int=str_replace_all(string=names_int,pattern="shape1", replacement="k_1")
    names_int=str_replace_all(string=names_int,pattern="shape2", replacement="k_2")
    names_int=str_replace_all(string=names_int,pattern="shape", replacement="k")
    names_int=str_replace_all(string=names_int,pattern="scale", replacement="s")
    
    prior.out = paste0("$\\pi_",1,"$=",formatC(prior,digits=2,format="f"))
    coeff.out = paste0("$",names_int,"$=",formatC(coeff,digits=2,format="f"),collapse=", ")
    
  } else{
    
    coeff.out = sapply(1:length(prior), function(i,coeff,prior){
      coeff.int = coeff[,i]
      names_int = rownames(coeff)
      
      names_int=str_replace_all(string=names_int,pattern="meanlog", replacement="\\\\mu")
      names_int=str_replace_all(string=names_int,pattern="sdlog", replacement="SD")
      names_int=str_replace_all(string=names_int,pattern="shape", replacement="k")
      names_int=str_replace_all(string=names_int,pattern="shape1", replacement="k_1")
      names_int=str_replace_all(string=names_int,pattern="shape2", replacement="k_2")
      names_int=str_replace_all(string=names_int,pattern="scale", replacement="s")
      
      
      if(i==length(prior)){
        out_prior = paste0("$\\pi_",i,"$=",formatC(prior[i],digits=2,format="f"))
        out_coeff = paste0("$",names_int,"$=",formatC(coeff.int,digits=2,format="f"),collapse=", ")
      }else{
        out_prior = paste0(paste0("$\\pi_",i,"$=",formatC(prior[i],digits=2,format="f"))," \\newline")
        out_coeff = paste0(paste0("$",names_int,"$=",formatC(coeff.int,digits=2,format="f"),collapse=", ")," \\newline")
      }
      
      return(c(out_prior=out_prior,out_coeff=out_coeff))
      
    },coeff=coeff,prior=prior)
    
    prior.out = paste0(coeff.out[1,], collapse=" ")
    coeff.out = paste0(coeff.out[2,], collapse=" ")
    
  }
  
  
  #Return values
  if (is.null(prior)) {
    return_list <- list(coefficients = coeff.out)
  }  else {
    if (nested) {
      return_list <- list(prior = list(prior.out), coefficients = list(coeff.out))
    }
    else {
      return_list <- list(prior = list(prior.out), coefficients = (coeff.out))
    }
  }
  return(return_list)
  
  
}