#' The Melitz-model
#'
#' This function calibrates the Melitz model
#' @param N number of countries
#' @param S number of sectors
#' @param delta exogenous death shock
#' @param alpha sectorial expenses, defaults to 1/S
#' @param sigma elasticity of substitution, defaults to 4
#' @param L Labour endowment
#' @param W wage
#' @param tau,ffixed variable and fixed cost, defaults to 1.5
#' @param fixede fixed entry cost, defaults to 0.545
#' @param phi cutoff
#' @param dist character vector indicating the used distribution, defaults to "pareto"
#' @param coeff coefficient list corresponding to the provided distribution, defaults to list(k=4.25, xmin=1)
#' 
#' @name melitz


#Examples
# test = melitz()

#' @rdname melitz
#' @export

melitz <- function(N=2, S=1, delta = 1, alpha =1/S, sigma = 4, L=1, W = 1, tau = 1.54, fixed = 1, 
                   fixede = 0.545, phi = 1, dist =c("pareto"),prior=1, coeff=list(k=4.25,xmin=1),
                   general=FALSE){
  
  
  require(rootSolve)
  require(distributionsrd)
  
  # Initialize values
  alpha = checkdim(x=alpha,x_name="alpha",N=N,S=S)
  sigma = checkdim(x=sigma,x_name="sigma",N=N,S=S) 
  delta = checkdim(x=delta,x_name="delta",N=N,S=S) 
  fixede = checkdim(x=fixede,x_name="fixed entry costs, fixede,",N=N,S=S)  
  
  L = checkdim(x=L,x_name="labour, L,",N=N) 
  W = checkdim(x=W,x_name="wages, W,",N=N)  
  
  tau = checkdim(x=tau,x_name="variable costs, tau,",N=N,S=S,destinationspecific = T) 
  fixed = checkdim(x=fixed,x_name="fixed costs, fixed,",N=N,S=S,destinationspecific = T) 
  
  if(N>1){
    tau = checkdiag(tau)
    fixed = checkdiag(fixed)
  }
  
  phi = checkdim(x=phi,x_name="cutoffs, phi,",N=N,S=S,destinationspecific = T)  
  
  varlist_destinationspecific = c("X","m_0","m_sigma","phi_av")
  varlist_destinationspecific = lapply(X=varlist_destinationspecific, FUN=function(i){
    return(array(NA,dim=c(N,N,S)))
  })
  names(varlist_destinationspecific) = c("X","m_0","m_sigma","phi_av")
  
  varlist = c("info","Me","M","Mx","P","R")
  varlist = lapply(X=varlist, FUN=function(i){
    return(array(NA,dim=c(N,1,S)))
  })
  names(varlist) =  c("info","Me","M","Mx","P","R")
  
  vars = c(varlist_destinationspecific,varlist)
  
  #Check feasibility of provided distribution
  check = as.numeric(sigma)
  check = sapply(1:length(sigma),function(i){
      return(suppressWarnings(mcombdist(r=(check[i]-1),dist=dist,prior=prior,coeff=coeff,truncation=0,lower.tail=F)))
    })
  if(any(is.na(check))){
      warning(paste0("Model can not be run as the (sigma-1) moment can't be calculated for distribution ",dist),call.=F)
      
    vars$phi = array(NA,dim=dim(phi))
      
    vars$phi = apply(vars$phi,c(1,2),sum)
    vars$P = apply(vars$P^alpha,c(1,2),prod)
    vars$M = apply(vars$M,c(1,2),sum)
    vars$R = apply(vars$R,c(1,2),sum)
    vars$X = apply(vars$X,c(1,2),sum)
    vars$m_0 = apply(vars$m_0,c(1,2),sum)
    vars$m_sigma = apply(vars$m_sigma,c(1,2),sum)
      
      return_list = vars
      return(return_list)
  }
  

  for(s in 1:S){
    out = melitz_equi_sector(N=N,S=S,delta = delta[,,s], alpha = alpha[,,s], sigma = sigma[,,s], L=L, W=W,
                             tau=tau[,,s], fixed = fixed[,,s], fixede=fixede[,,s], 
                             phi=phi[,,s],dist=dist,prior=prior,coeff=coeff)
    #Fill out next spots for better starting values too
    vars$phi[,,(s:S)] = out$phi
    vars$phi_av[,,s] = out$phi_av
    vars$Me[,,s] = out$Me
    vars$M[,,s] = out$M
    vars$Mx[,,s] = out$Mx
    vars$P[,,s] = out$P
    vars$R[,,s] = out$R
    vars$X[,,s] = out$X
    vars$m_0[,,s] = out$m_0
    vars$m_sigma[,,s] = out$m_sigma
    vars$info[,,s] = out$info
  }

    if(S==1){
      vars$phi = apply(vars$phi,c(1,2),sum)
      vars$P = apply(vars$P^alpha,c(1,2),prod)
      vars$M = apply(vars$M,c(1,2),sum)
      vars$R = apply(vars$R,c(1,2),sum)
      vars$X = apply(vars$X,c(1,2),sum)
      vars$m_0 = apply(vars$m_0,c(1,2),sum)
      vars$m_sigma = apply(vars$m_sigma,c(1,2),sum)
    }

  return_list = vars
  return(return_list)
  
}


#' @rdname melitz
#' @export

melitz_equi_sector = function(delta, alpha, sigma , L, W, tau, fixed , fixede, phi,dist,prior,coeff,N,S){
  
   # relation between domestic and exporting cutoffs
  if(N==1){
    cutoff_rel = matrix(1,1,1)
  }else{
    cutoff_rel = (W%*%t(1/W))^(sigma/(sigma - 1)) *  (fixed/diag(fixed))^(1/(sigma - 1)) *  (tau/diag(tau))
  }
  
  if(dist=="empirical"){
    
    assign("m_0_f",mcombdist(r=0,dist=dist,coeff=coeff,truncation=NULL,lower.tail=F),env=.GlobalEnv)
    sapply(1:length(sigma),function(i){
      assign(paste0("m_sigma_f_",i),mcombdist(r=(sigma[i]-1),dist=dist,coeff=coeff,lower.tail=F),env=.GlobalEnv)
    })
    
  }
  
  #Retrieve cutoffs from Free entry equilibrium
  roots = multiroot(melitz_cutoffs, start = diag(phi), cutoff_rel = cutoff_rel, sigma=sigma,
                    fixede=fixede,fixed=fixed, delta=delta, dist=dist,prior=prior, coeff = coeff)
  phi = roots$root * cutoff_rel
  
  if(dist=="empirical"){
    
    m_0 =  m_0_f(phi)
    m_0 = matrix(m_0,dim(phi))
    
    m_sigma = sapply(1:length(sigma),function(i){
      return(do.call(paste0("m_sigma_f_",i),list(phi[i,])))
    })
    m_sigma = t(m_sigma)
    
  }else{
    
    m_0 =  suppressWarnings(mcombdist(r=0,dist=dist,prior=prior,coeff=coeff,truncation=phi,lower.tail=F))
    m_0 = matrix(m_0,dim(phi))
    
    m_sigma = sapply(1:length(sigma),function(i){
      return(suppressWarnings(mcombdist(r=(sigma[i]-1),dist=dist,prior=prior,coeff=coeff,
                                        truncation=phi[i,],lower.tail=F)))
    })
    m_sigma = t(m_sigma)
    
  }
  
  #Mass of firms
  Me = c((delta * alpha * L) / rowSums((phi^(1-sigma)*m_sigma * sigma *fixed )))
  M = Me * m_0 / delta
  Mx = rowSums(M) - diag(M)
  M = diag(M)
  
  # Price index
  P = rowSums((sigma/(sigma - 1))^(1 - sigma) * Me / delta * (W * tau)^(1 - sigma) * m_sigma)^(1/(1 - sigma))

  # Output
  X = Me/delta * sigma * W * fixed * phi^(1-sigma) * m_sigma
  R = rowSums(X)
  
  phi_av = (delta*fixede)/m_0
  
  return_list = list(phi = phi,Me = Me,M = M,Mx = Mx,P = P,R = R, X=X, m_0 = m_0, m_sigma = m_sigma,
                     info = roots$estim.precis, phi_av = phi_av)
  
}

#' @rdname melitz
#' @export

melitz_cutoffs = function(phi_dom, cutoff_rel, sigma, fixede, fixed, dist,prior, coeff,delta){
  
  phi = phi_dom * cutoff_rel
  
  if(dist=="empirical"){
    
    m_0 =  m_0_f(phi)
    m_0 = matrix(m_0,dim(phi))
    
    m_sigma = sapply(1:length(sigma),function(i){
      return(do.call(paste0("m_sigma_f_",i),list(phi[i,])))
    })
    m_sigma = t(m_sigma)
    
  }else{
    
    m_0 =  suppressWarnings(mcombdist(r=0,dist=dist,prior=prior,coeff=coeff,truncation=phi,lower.tail=F))
    m_0 = matrix(m_0,dim(phi))
    
    m_sigma = sapply(1:length(sigma),function(i){
      return(suppressWarnings(mcombdist(r=(sigma[i]-1),dist=dist,prior=prior,coeff=coeff,
                                        truncation=phi[i,],lower.tail=F)))
    })
    m_sigma = t(m_sigma)
    
    
  }
  
  J = phi^(1-sigma) * m_sigma - m_0
  
  return(fixede*delta - rowSums(fixed*J))
  
}




