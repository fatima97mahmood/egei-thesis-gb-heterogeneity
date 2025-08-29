##### Adding code needed for incorporation of new distributions

`%||%` <- function(a, b) if (!is.null(a)) a else b

#################################################################################################################
# 1 Parameter Distributions 
#################################################################################################################

#################################################################################################################
# Loglogistic
#################################################################################################################
# Log-Logistic PDF (Fisk) via GB parameters: a=α, c=1, p=1, q=1, b=β
dloglogis <- function(x, alpha, beta) {
  (alpha / beta) * (x / beta)^(alpha - 1) / (1 + (x / beta)^alpha)^2
}
dloglogistic <- function(x, alpha, beta = 1) {
  (alpha/beta) * (x/beta)^(alpha-1) / (1 + (x/beta)^alpha)^2
}

# Log-likelihood
loglik_loglogis <- function(params, x){
  alpha <- params[1]; beta <- params[2]
  if (alpha <= 0 || beta <= 0) return(1e6)
  if (any(x <= 0)) return(1e6)
  -sum(log(alpha) - log(beta) + (alpha - 1)*log(x/beta) - 2*log1p((x/beta)^alpha))
}


# MLE
loglogis.mle <- function(x) {
  optim(c(alpha = 1, beta = median(x)), loglik_loglogis, x = x,
        method = "L-BFGS-B", lower = c(1e-8, 1e-8))
}

#################################################################################################################
# Rayleigh
#################################################################################################################
# Custom Rayleigh functions
drayleigh <- function(x, sigma) {
  if (sigma <= 0) return(rep(NA, length(x)))
  (x / sigma^2) * exp(-x^2 / (2 * sigma^2)) * (x >= 0)
}

loglik_rayleigh <- function(sigma, x){
  if (sigma <= 0) return(1e6)
  if (any(x < 0)) return(1e6)
  -sum(log(x) - 2*log(sigma) - (x^2)/(2*sigma^2))
}

rayleigh.mle <- function(x) {
  sigma_hat <- sqrt(mean(x^2) / 2)  # MLE for Rayleigh
  loglik <- loglik_rayleigh(sigma_hat, x)
  list(sigma = sigma_hat, loglik = loglik)
}

#################################################################################################################
# Uniform
#################################################################################################################
# Uniform PDF (1-parameter, derived from Power with θ=1)
dunif1 <- function(y, beta) {
  if (beta <= 0) stop("beta must be positive")
  pdf <- rep(1/beta, length(y))
  pdf[y < 0 | y > beta] <- 0  # Support: [0, β]
  return(pdf)
}

# Log-likelihood function
loglik_unif1 <- function(beta, x) {
  if (beta <= 0) return(1e6)
  if (any(x < 0 | x > beta)) return(1e6)  # All data must be in [0, β]
  -sum(log(1/beta))  # Negative log-likelihood
}

# MLE function
unif1.mle <- function(x, ...) {
  # MLE for β is simply max(x) (optimization not needed)
  beta_hat <- max(x)
  list(
    par = c(beta = beta_hat),
    value = loglik_unif1(beta_hat, x),
    convergence = 0,
    message = "MLE is analytical (max(x))"
  )
}

#################################################################################################################
# t-dist
#################################################################################################################
# t-Distribution PDF
dt_1param <- function(x, nu) {
  # Standard t-distribution (μ=0, σ=1)
  (gamma((nu + 1)/2) / (sqrt(nu * pi) * gamma(nu/2))) * 
    (1 + (x^2)/nu)^(-(nu + 1)/2)
}

# Log-likelihood function
loglik_t1 <- function(nu, x) {
  if (nu <= 0) return(1e10)  # Penalize invalid ν
  -sum(log(dt_1param(x, nu)))
}

t1.mle <- function(x) {
  # Initialize ν using moment matching (ν ≈ 2*σ² / (σ² - 1))
  sigma2 <- var(x)
  nu_start <- ifelse(sigma2 > 1, 2*sigma2 / (sigma2 - 1), 3)  # Fallback if σ² ≤ 1
  
  optim(
    par = nu_start,
    fn = loglik_t1,
    x = x,
    method = "L-BFGS-B",
    lower = 1e-8,
    upper = 100  # Constrain ν to avoid extreme tails
  )
}

#################################################################################################################
# 2 Parameter Distributions 
#################################################################################################################

#################################################################################################################
# Power
#################################################################################################################
dpower <- function(y, beta, theta) {
  # Input validation
  if (beta <= 0 || theta <= 0) {
    stop("Both beta and theta must be positive")
  }
  
  # Initialize pdf vector
  pdf <- numeric(length(y))
  
  # Calculate pdf for values in [0, beta]
  in_support <- (y >= 0) & (y <= beta)
  pdf[in_support] <- theta * (y[in_support]^(theta - 1)) / (beta^theta)
  
  # Values outside support have pdf = 0
  pdf[!in_support] <- 0
  
  return(pdf)
}

# Power log-likelihood function
loglik_power <- function(params, x){
  beta <- params[1]; theta <- params[2]
  if (beta <= 0 || theta <= 0) return(1e6)
  if (any(x < 0 | x > beta)) return(1e6)
  # all x are in [0,beta] now
  if (any(x == 0 & theta < 1)) return(1e6)  # avoid -Inf when theta<1
  -sum(log(theta) + (theta - 1)*log(x) - theta*log(beta))
}


# Power MLE function
power.mle <- function(x, start = NULL, lower = NULL, upper = NULL, ...) {
  if (is.null(start)) {
    # Reasonable starting values
    start <- c(beta = max(x), theta = 1)
  }
  
  if (is.null(lower)) {
    lower <- c(beta = 1e-8, theta = 1e-8)
  }
  
  if (is.null(upper)) {
    upper <- c(beta = Inf, theta = Inf)
  }
  
  # Estimate parameters
  res <- optim(
    par = start,
    fn = loglik_power,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
  
  return(res)
}


#################################################################################################################
# Fisk
#################################################################################################################
# Fisk PDF (Log-Logistic)
dfisk <- function(y, beta, theta) {
  if (beta <= 0 || theta <= 0) stop("beta and theta must be positive")
  term <- (y / beta)^theta
  pdf <- (theta / beta) * term^(1 - 1/theta) / (1 + term)^2
  pdf[y <= 0] <- 0  # Support: y > 0
  return(pdf)
}

# Log-likelihood 
loglik_fisk <- function(params, x) {
  beta <- params[1]
  theta <- params[2]
  if (beta <= 0 || theta <= 0) return(1e6)
  if (any(x <= 0)) return(1e6)
  term <- (x / beta)^theta
  log_pdf <- log(theta / beta) + (theta - 1) * log(x / beta) - 2 * log(1 + term)
  -sum(log_pdf)  # Return negative log-likelihood
}

# MLE Estimator
fisk.mle <- function(x, start = c(beta = median(x), theta = 1), ...) {
  optim(
    par = start,
    fn = loglik_fisk,
    x = x,
    method = "L-BFGS-B",
    lower = c(beta = 1e-8, theta = 1e-8),
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# Lomax
#################################################################################################################
# Lomax PDF (GB parameterization: a=1, c=1, p=1, q=θ)
dlomax <- function(y, beta, theta) {
  if (beta <= 0 || theta <= 0) stop("beta and theta must be positive")
  pdf <- theta / (beta * (1 + y/beta)^(theta + 1))
  pdf[y < 0] <- 0  # Support: y >= 0
  return(pdf)
}

# Log-likelihood 
loglik_lomax <- function(params, x) {
  beta <- params[1]
  theta <- params[2]
  if (beta <= 0 || theta <= 0) return(1e6)
  -sum(log(theta) - log(beta) - (theta + 1) * log(1 + x/beta))
}

# MLE 
lomax.mle <- function(x, ...) {
  optim(c(beta = 1, theta = 1), loglik_lomax, x = x, 
        method = "L-BFGS-B", lower = c(1e-8, 1e-8))
}

#################################################################################################################
# Inverse Lomax
#################################################################################################################
# Inverse Lomax PDF (GB parameters: a=-1, c=1, p=1, q=θ)
dinvlomax <- function(y, beta, theta){
  if (beta <= 0 || theta <= 0) stop("beta and theta must be positive")
  pdf <- (theta * beta^theta) / ( y^(theta+1) * (1 + beta/y)^(theta+1) )
  pdf[y <= 0] <- 0
  pdf
}

#loglik invlomax
loglik_invlomax <- function(params, x){
  beta <- params[1]; theta <- params[2]
  if (beta <= 0 || theta <= 0) return(1e6)
  if (any(x <= 0)) return(1e6)
  -sum( log(theta) + theta*log(beta) - (theta+1)*log(x) - (theta+1)*log1p(beta/x) )
}


# MLE function 
invlomax.mle <- function(x, start = c(beta = 1, theta = 1), ...) {
  optim(
    par = start,
    fn = loglik_invlomax,
    x = x,
    method = "L-BFGS-B",
    lower = c(beta = 1e-8, theta = 1e-8),
    control = list(maxit = 5000)
  )
}


#################################################################################################################
# 3 Parameter Distributions
#################################################################################################################

#################################################################################################################
# Inverse Beta Type 1 (IB1)
#################################################################################################################
# (3-parameter: b, p, q)
# Fixed: a=-1, c=0
# Density function (from GB formula with a=-1, c=0)
dib1 <- function(y, b, p, q) {
  num <- y^(-p - 1) * (1 - (1 - 0)*(y/b)^-1)^(q - 1)
  denom <- b^(-p) * beta(p, q) * (1 + 0*(y/b)^-1)^(p + q)
  pdf <- num / denom
  pdf[(y <= 0) | (y^-1 >= b^-1/(1 - 0))] <- 0  # Support: y > b
  return(pdf)
}

# Simplified form (after algebraic reduction)
dib1 <- function(y, b, p, q) {
  pdf <- (y^(-p - 1) * (1 - b/y)^(q - 1)) / (b^(-p) * beta(p, q))
  pdf[y <= b] <- 0  # Support: y > b
  return(pdf)
}

# Log-likelihood
loglik_ib1 <- function(params, x) {
  b <- params[1]
  p <- params[2]
  q <- params[3]
  if (b <= 0 || p <= 0 || q <= 0) return(1e6)
  
  y <- x[x > b]  # Support: y > b
  if (length(y) == 0) return(1e6)
  
  pdf <- dib1(y, b, p, q)
  if (any(pdf <= 0 | !is.finite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

# MLE fitting
ib1.mle <- function(x, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) start <- c(b = min(x)/2, p = 1, q = 1)
  if (is.null(lower)) lower <- c(b = 1e-8, p = 1e-8, q = 1e-8)
  if (is.null(upper)) upper <- c(b = max(x)-1e-8, p = Inf, q = Inf)
  
  optim(
    par = start,
    fn = loglik_ib1,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# Unit Gamma (UG) Distribution
#################################################################################################################
# Parameters: b > 0, d > 0, q > 0
# Derived from GB with c=0, a→0, p=d/a
# UG density function
dug <- function(y, b, d, q){
  if (b <= 0 || d <= 0 || q <= 0) return(rep(0, length(y)))
  log_pdf <- log(q) - (d/q)*log(b) - lgamma(d/q) + (d-1)*log(log(y)) - (log(y)/b)^q - log(y)
  pdf <- exp(log_pdf)
  pdf[y <= 1 | !is.finite(pdf)] <- 0
  pdf
}


# Log-likelihood
loglik_ug <- function(params, x){
  b <- params[1]; d <- params[2]; q <- params[3]
  if (b <= 0 || d <= 0 || q <= 0) return(1e10)
  if (any(x <= 1)) return(1e10)
  log_pdf <- log(q) - (d/q)*log(b) - lgamma(d/q) + (d-1)*log(log(x)) - (log(x)/b)^q - log(x)
  if (any(!is.finite(log_pdf))) return(1e10)
  -sum(log_pdf)
}

# Fitting function (matches your 3-param pattern)
ug.mle <- function(x, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) start <- c(b = 1, d = 1, q = 1)
  if (is.null(lower)) lower <- c(b = 1e-8, d = 1e-8, q = 1e-8)
  if (is.null(upper)) upper <- c(b = Inf, d = Inf, q = Inf)
  
  optim(
    par = start,
    fn = loglik_ug,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# B1 Beta 1 distribution
#################################################################################################################
# Beta 1 PDF
db1 <- function(y, b, p, q) {
  beta_pq <- beta(p, q)
  pdf <- (y^(p - 1) * (1 - (y/b))^(q - 1)) / (b^p * beta_pq)
  pdf[(y <= 0) | (y >= b)] <- 0
  return(pdf)
}

loglik_b1 <- function(params, x) {
  b <- params[1]; p <- params[2]; q <- params[3]
  if (any(x <= 0 | x >= b)) return(1e6)
  if (b <= 0 || p <= 0 || q <= 0) return(1e6)
  y <- x[(x > 0) & (x < b)]
  if (length(y) == 0) return(1e6)
  pdf <- db1(y, b, p, q)
  if (any(pdf <= 0 | !is.finite(pdf))) return(1e6)  # Combined check
  return(-sum(log(pdf)))
}

b1.mle <- function(x, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) start <- c(b = max(x)*1.01, p = 1, q = 1)  # Slightly safer b
  if (is.null(lower)) lower <- c(1e-8, 1e-8, 1e-8)
  if (is.null(upper)) upper <- c(Inf, Inf, Inf)
  
  optim(par = start, fn = loglik_b1, x = x,
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = 5000))
}


#################################################################################################################
# GG Generalized Gamma Distribution 
#################################################################################################################
# 3-parameter version matching your format
# Parameters: b (scale), p (shape), q (shape)
#=============================================

# GG density (3-param form)
dgg <- function(y, b, p, q) {
  num <- abs(q) * y^(q*p - 1) * exp(-(y/b)^q)
  denom <- b^(q*p) * gamma(p)
  pdf <- num / denom
  pdf[y <= 0] <- 0
  return(pdf)
}

# Log-likelihood 
loglik_gg <- function(params, x) {
  b <- params[1]
  p <- params[2]
  q <- params[3]
  
  if (b <= 0 || p <= 0 || q <= 0) return(1e10)
  if (any(x <= 0)) return(1e10)
  y <- x[x > 0]
  if (length(y) == 0) return(1e10)
  
  log_pdf <- log(abs(q)) + (q*p - 1)*log(y) - (y/b)^q - q*p*log(b) - lgamma(p)
  if (any(!is.finite(log_pdf))) return(1e10)
  -sum(log_pdf)
}

# Fitting function (matches your style)
gg.mle <- function(x, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) start <- c(b = median(x), p = 1, q = 1)
  if (is.null(lower)) lower <- c(b = 1e-8, p = 1e-8, q = 1e-8)
  if (is.null(upper)) upper <- c(b = Inf, p = Inf, q = Inf)
  
  optim(
    par = start,
    fn = loglik_gg,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# B2 Beta 2 Distribution
#################################################################################################################
# PDF for B2
db2 <- function(y, b, p, q) {
  beta_pq <- beta(p, q)
  num <- y^(p - 1)
  denom <- b^p * beta_pq * (1 + y / b)^(p + q)
  pdf <- num / denom
  pdf[y <= 0] <- 0
  return(pdf)
}

loglik_b2 <- function(params, x) {
  b <- params[1]; p <- params[2]; q <- params[3]
  if (any(x <= 0)) return(1e6)
  if (b <= 0 || p <= 0 || q <= 0) return(1e6)
  y <- x[x > 0]
  if (length(y) == 0) return(1e6)
  pdf <- db2(y, b, p, q)
  if (any(pdf <= 0 | is.nan(pdf) | is.infinite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

b2.mle <- function(x, start = NULL, lower = NULL, upper = NULL, ...) {
  if (is.null(start)) {
    start <- c(b = mean(x), p = 2, q = 2)
  }
  if (is.null(lower)) {
    lower <- c(1e-8, 1e-8, 1e-8)
  }
  if (is.null(upper)) {
    upper <- c(Inf, Inf, Inf)
  }
  
  optim(
    par = start,
    fn = loglik_b2,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# BR12 Beta Residual 12 distribution 
#################################################################################################################
# (3 free parameters: a, b, q)
# Fixed: c=1, p=1
# =============================================

# Density function (GB with c=1, p=1)
dbr12 <- function(y, a, b, q) {
  beta_1q <- beta(1, q)  # beta(1,q) = 1/q
  num <- abs(a) * y^(a*1 - 1) * (1 - 0*(y/b)^a)^(q - 1)  # (1-c)=0 since c=1
  denom <- b^(a*1) * beta_1q * (1 + 1*(y/b)^a)^(1 + q)
  pdf <- num / denom
  pdf[(y <= 0) | (y^a >= b^a/(1 - 1))] <- 0  # Support: 0 < y^a < b^a/0 → y > 0
  return(pdf)
}

# Log-likelihood
loglik_br12 <- function(params, x) {
  a <- params[1]
  b <- params[2]
  q <- params[3]
  if (b <= 0 || q <= 0 || a == 0) return(1e6)  # a can be negative
  
  y <- x[x > 0]  # Support: y > 0 (since c=1)
  if (length(y) == 0) return(1e6)
  if (any(x <= 0)) return(1e6)
  pdf <- dbr12(y, a, b, q)
  if (any(pdf <= 0 | !is.finite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

# MLE fitting
br12.mle <- function(x, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) start <- c(a = 1, b = median(x), q = 1)
  if (is.null(lower)) lower <- c(a = -Inf, b = 1e-8, q = 1e-8)
  if (is.null(upper)) upper <- c(a = Inf, b = Inf, q = Inf)
  
  optim(
    par = start,
    fn = loglik_br12,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# BR3 Beta Residual 3 distribution
#################################################################################################################
# (3-parameter: a, b, p)
# Fixed: c=1, q=1
# =============================================

# Density function (from GB formula with c=1, q=1)
dbr3 <- function(y, a, b, p) {
  num <- abs(a) * y^(a*p - 1)
  denom <- b^(a*p) * beta(p, 1) * (1 + (y/b)^a)^(p + 1)
  pdf <- num / denom
  pdf[(y <= 0) | (y^a >= b^a/(1 - 1))] <- 0  # Support: y > 0 (since c=1)
  return(pdf)
}

# Log-likelihood
loglik_br3 <- function(params, x) {
  a <- params[1]
  b <- params[2]
  p <- params[3]
  if (b <= 0 || p <= 0 || a == 0) return(1e6)
  if (any(x <= 0)) return(1e6)
  y <- x[x > 0]  # Support: y > 0
  if (length(y) == 0) return(1e6)
  
  pdf <- dbr3(y, a, b, p)
  if (any(pdf <= 0 | !is.finite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

# MLE fitting
br3.mle <- function(x, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) start <- c(a = 1, b = median(x), p = 1)
  if (is.null(lower)) lower <- c(a = -Inf, b = 1e-8, p = 1e-8)
  if (is.null(upper)) upper <- c(a = Inf, b = Inf, p = Inf)
  
  optim(
    par = start,
    fn = loglik_br3,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# 4 parameter distributions
#################################################################################################################

#################################################################################################################
# Generalized Beta Type I (GB1) Distribution
#################################################################################################################

# GB1 PDF
dgb1 <- function(y, a, b, p, q) {
  beta_pq <- beta(p, q)
  num <- abs(a) * y^(a * p - 1) * (1 - (y / b)^a)^(q - 1)
  denom <- b^(a * p) * beta_pq
  pdf <- num / denom
  pdf[(y <= 0) | (y >= b)] <- 0
  return(pdf)
}

loglik_gb1 <- function(params, x) {
  a <- params[1]; b <- params[2]; p <- params[3]; q <- params[4]
  if (any(x <= 0 | x >= b)) return(1e6)
  if (b <= 0 || p <= 0 || q <= 0 || a == 0) return(1e6)
  y <- x[(x > 0) & (x < b)]
  if (length(y) == 0) return(1e6)
  pdf <- dgb1(y, a, b, p, q)
  if (any(pdf <= 0 | is.nan(pdf) | is.infinite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

# GB1 MLE function
gb1.mle <- function(x, start = NULL, lower = NULL, upper = NULL, ...) {
  if (is.null(start)) {
    start <- c(a = 1, b = mean(x), p = 2, q = 2)
  }
  if (is.null(lower)) {
    lower <- c(a = -Inf, b = 1e-8, p = 1e-8, q = 1e-8)
  }
  if (is.null(upper)) {
    upper <- c(a = Inf, b = Inf, p = Inf, q = Inf)
  }
  
  optim(
    par = start,
    fn = function(params) {
      a <- params[1]; b <- params[2]; p <- params[3]; q <- params[4]
      loglik_gb(c(a, b, 0, p, q), x)  # Fix c = 0
    },
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}



#################################################################################################################
# Generalized Beta of the Second Kind (GB2)
#################################################################################################################
# GB2 PDF
dgb2 <- function(y, a, b, p, q) {
  beta_pq <- beta(p, q)
  num <- abs(a) * y^(a * p - 1)
  denom <- b^(a * p) * beta_pq * (1 + (y / b)^a)^(p + q)
  pdf <- num / denom
  pdf[y <= 0] <- 0
  return(pdf)
}
loglik_gb2 <- function(params, x) {
  a <- params[1]; b <- params[2]; p <- params[3]; q <- params[4]
  if (b <= 0 || p <= 0 || q <= 0 || a == 0) return(1e6)
  y <- x[x > 0]
  if (length(y) == 0) return(1e6)
  pdf <- dgb2(y, a, b, p, q)
  if (any(pdf <= 0 | is.nan(pdf) | is.infinite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

# GB2 MLE
gb2.mle <- function(x, start = NULL, lower = NULL, upper = NULL, ...) {
  if (is.null(start)) {
    start <- c(a = 1, b = mean(x), p = 2, q = 2)
  }
  if (is.null(lower)) {
    lower <- c(a = -Inf, b = 1e-8, p = 1e-8, q = 1e-8)
  }
  if (is.null(upper)) {
    upper <- c(a = Inf, b = Inf, p = Inf, q = Inf)
  }
  
  optim(
    par = start,
    fn = function(params) {
      a <- params[1]; b <- params[2]; p <- params[3]; q <- params[4]
      loglik_gb(c(a, b, 1, p, q), x)  # Fix c = 1
    },
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}


#################################################################################################################
# Beta4
#################################################################################################################
# PDF
dbeta4 <- function(y, b, c, p, q) {
  beta_pq <- beta(p, q)
  num <- y^(p - 1) * (1 - (1 - c)*(y / b))^(q - 1)
  denom <- b^p * beta_pq * (1 + c * (y / b))^(p + q)
  pdf <- num / denom
  pdf[(y <= 0) | (y >= b / (1 - c))] <- 0
  return(pdf)
}


loglik_beta4 <- function(params, x) {
  b <- params[1]; c <- params[2]; p <- params[3]; q <- params[4]
  if (any(!(x > 0 & x < b/(1-c)))) return(1e6)
  if (b <= 0 || c < 0 || c >= 1 || p <= 0 || q <= 0) return(1e6)
  support_ok <- (x > 0) & (x < b / (1 - c))
  if (sum(support_ok) == 0) return(1e6)
  y <- x[support_ok]
  pdf <- dbeta4(y, b, c, p, q)
  if (any(pdf <= 0 | is.nan(pdf) | is.infinite(pdf))) return(1e6)
  return(-sum(log(pdf)))
}

# mle

beta4.mle <- function(x, start = NULL, lower = NULL, upper = NULL, ...) {
  if (is.null(start)) {
    start <- c(b = mean(x), c = 0.5, p = 2, q = 2)
  }
  if (is.null(lower)) {
    lower <- c(1e-8, 0, 1e-8, 1e-8)
  }
  if (is.null(upper)) {
    upper <- c(Inf, 1 - 1e-8, Inf, Inf)
  }
  
  optim(
    par = start,
    fn = function(params) {
      b <- params[1]; c <- params[2]; p <- params[3]; q <- params[4]
      loglik_gb(c(1, b, c, p, q), x)  # Fix a = 1
    },
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
}

#################################################################################################################
# 5 parameter distribution
#################################################################################################################
# Generalized Beta (5-parameter)
dgb <- function(y, a, b, c, p, q) {
  beta_pq <- beta(p, q)
  num <- abs(a) * y^(a*p - 1) * (1 - (1 - c) * (y/b)^a)^(q - 1)
  denom <- (b^(a*p)) * beta_pq * (1 + c * (y/b)^a)^(p + q)
  pdf <- num / denom
  
  # Set pdf = 0 for values outside support
  pdf[(y^a <= 0) | (y^a >= (b^a)/(1 - c))] <- 0
  
  return(pdf)
}

loglik_gb <- function(params, x) {
  a <- params[1]
  b <- params[2]
  c <- params[3]
  p <- params[4]
  q <- params[5]
  
  # Check parameter boundaries
  if (b <= 0 || c < 0 || c > 1 || p <= 0 || q <= 0 || a == 0) {
    return(1e6)
  }
  
  support_ok <- (x^a > 0) & (x^a < (b^a)/(1 - c))
  if (sum(support_ok) == 0) {
    return(1e6)
  }
  
  y <- x[support_ok]
  beta_pq <- beta(p, q)
  num <- abs(a) * y^(a * p - 1) * (1 - (1 - c) * (y / b)^a)^(q - 1)
  denom <- (b^(a * p)) * beta_pq * (1 + c * (y / b)^a)^(p + q)
  pdf <- num / denom
  
  if (any(pdf <= 0 | is.nan(pdf) | is.infinite(pdf))) {
    return(1e6)
  }
  
  return(-sum(log(pdf)))
}

# mle function
gb.mle <- function(x, start = NULL, lower = NULL, upper = NULL, ...) {
  if (is.null(start)) {
    start <- c(a = 1, b = mean(x), c = 0.5, p = 2, q = 2)
  }
  
  if (is.null(lower)) {
    lower <- c(a = 1e-5, b = 1e-8, c = 0, p = 1e-8, q = 1e-8)
  }
  
  if (is.null(upper)) {
    upper <- c(a = Inf, b = Inf, c = 1, p = Inf, q = Inf)
  }
  
  res <- optim(
    par = start,
    fn = loglik_gb,
    x = x,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 5000)
  )
  
  return(res)
}

#################################################################################################################
# Updating combdist functions
#################################################################################################################
create_failed_result_1param <- function(dist, components, n) {
  params <- list(
    unif1 = "beta",
    exponential = "rate",
    rayleigh = "sigma",
    chisq = "df",
    loglogistic = "alpha",
    normal = "sd",
    t = "df"
  )[[dist]]
  
  list(
    dist = dist,
    components = components,
    prior = 1,
    coefficients = setNames(NA_real_, params),
    np = length(params),
    n = n,
    convergence = FALSE,
    loglik = NA_real_,
    NLL = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    message = "Fit failed"
  )
}

create_failed_result_2param <- function(dist, components, n) {
  params <- list(
    power = c("beta", "theta"),
    fisk = c("beta", "theta"),
    lomax = c("beta", "theta"),
    invlomax = c("beta", "theta")
  )[[dist]]
  
  list(
    dist = dist,
    components = components,
    prior = 1,
    coefficients = setNames(rep(NA_real_, length(params)), params),
    np = length(params),
    n = n,
    convergence = FALSE,
    loglik = NA_real_,
    NLL = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    message = "Fit failed"
  )
}

create_failed_result_3param <- function(dist, components, n) {
  params <- list(
    b1 = c("b", "p", "q"),
    b2 = c("b", "p", "q"),
    br12 = c("a", "b", "q"),
    br3 = c("a", "b", "p"),
    ib1 = c("b", "p", "q"),
    gg = c("b", "p", "q"),
    ug = c("b", "d", "q")
  )[[dist]]
  
  list(
    dist = dist,
    components = components,
    prior = 1,
    coefficients = setNames(rep(NA_real_, length(params)), params),
    np = length(params),
    n = n,
    convergence = FALSE,
    loglik = NA_real_,
    NLL = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    message = "Fit failed"
  )
}

create_failed_result_4_5_param <- function(dist, components, n) {
  params <- list(
    gb = c("a", "b", "c", "p", "q"),
    gb1 = c("a", "b", "p", "q"),
    gb2 = c("a", "b", "p", "q"),
    beta4 = c("b", "c", "p", "q")
  )[[dist]]
  
  list(
    dist = dist,
    components = components,
    prior = 1,
    coefficients = setNames(rep(NA_real_, length(params)), params),
    np = length(params),
    n = n,
    convergence = FALSE,
    loglik = NA_real_,
    NLL = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    message = "Fit failed"
  )
}

# Update the combdist functions to include NLL, AIC, BIC calculations for consistency
update_combdist_output <- function(fit_result, dist, components, n) {
  if (is.null(fit_result) || is.null(fit_result$par)) {
    # choose factory by expected param count (simple map)
    np_map <- list(
      unif1=1, exponential=1, rayleigh=1, chisq=1, loglogistic=1, normal=1, t=1,
      power=2, fisk=2, lomax=2, invlomax=2,
      b1=3, b2=3, br12=3, br3=3, ib1=3, gg=3, ug=3,
      gb=5, gb1=4, gb2=4, beta4=4
    )
    np_guess <- np_map[[dist]] %||% 1
    return(if (np_guess==1) create_failed_result_1param(dist, components, n)
           else if (np_guess==2) create_failed_result_2param(dist, components, n)
           else if (np_guess==3) create_failed_result_3param(dist, components, n)
           else create_failed_result_4_5_param(dist, components, n))
  }
  coefficients <- fit_result$par
  loglik <- ifelse(is.null(fit_result$value), NA, -fit_result$value)
  
  list(
    dist = dist,
    components = components,
    prior = 1,
    coefficients = coefficients,
    np = length(coefficients),
    n = n,
    convergence = ifelse(is.null(fit_result$convergence), 
                         FALSE,  # Default to FALSE if missing
                         fit_result$convergence == 0),  # Convert to logical
    loglik = loglik,
    NLL = ifelse(is.na(loglik), NA, -loglik),
    AIC = ifelse(is.na(loglik), NA, 2 * length(coefficients) - 2 * loglik),
    BIC = ifelse(is.na(loglik), NA, log(n) * length(coefficients) - 2 * loglik),
    message = ifelse(is.null(fit_result$message),
                     ifelse(is.null(fit_result$convergence) || fit_result$convergence != 0,
                            "Convergence issue", "Success"),
                     fit_result$message)
  )
}


#################################################################################################################
# 1 parameter distributions
#################################################################################################################
combdist.mle_1param <- function(x, dist, start = NULL, lower = NULL, upper = NULL,
                                components = 1, nested = FALSE, steps = 1,
                                lowertrunc = 0, uppertrunc = Inf, ...) {
  # Supported 1-parameter distributions and their information
  dist_info <- list(
    unif1 = list(
      params = "beta",
      start_fun = function(x) max(x),
      fit_fun = function(x) unif1.mle(x)
    ),
    exponential = list(
      params = "rate",
      start_fun = function(x) 1/mean(x),
      fit_fun = function(x) {
        rate_hat <- 1/mean(x)
        list(par = c(rate = rate_hat),
             value = -sum(dexp(x, rate = rate_hat, log = TRUE)),
             convergence = 0,
             message = "MLE is analytical (1/mean(x))")
      }
    ),
    rayleigh = list(
      params = "sigma",
      start_fun = function(x) sqrt(mean(x^2)/2),
      fit_fun = function(x) {
        sigma_hat <- sqrt(mean(x^2)/2)
        list(par = c(sigma = sigma_hat),
             value = -sum(log(x) - 2*log(sigma_hat) - (x^2)/(2*sigma_hat^2)),
             convergence = 0,
             message = "MLE is analytical")
      }
    ),
    chisq = list(
      params = "df",
      start_fun = function(x) mean(x),
      fit_fun = function(x) {
        df_hat <- mean(x)
        list(par = c(df = df_hat),
             value = -sum(dchisq(x, df = df_hat, log = TRUE)),
             convergence = 0,
             message = "MLE is analytical (mean(x))")
      }
    ),
    loglogistic = list(
      params = "alpha",
      start_fun = function(x) 1,
      fit_fun = function(x) {
        # Simplified 1-parameter version (fix beta=1)
        loglik <- function(alpha, x) {
          if (alpha <= 0) return(1e6)
          -sum(log(alpha) + (alpha - 1)*log(x) - 2*log(1 + x^alpha))
        }
        opt <- optimize(loglik, c(1e-8, 100), x = x)
        list(par = c(alpha = opt$minimum),
             value = opt$objective,
             convergence = 0,
             message = "Optimization complete")
      }
    ),
    normal = list(
      params = "sd",
      start_fun = function(x) sd(x),
      fit_fun = function(x) {
        sd_hat <- sd(x)
        list(par = c(sd = sd_hat),
             value = -sum(dnorm(x, mean = 0, sd = sd_hat, log = TRUE)),
             convergence = 0,
             message = "MLE is analytical (sd(x))")
      }
    ),
    t = list(
      params = "df",
      start_fun = function(x) 10,
      fit_fun = function(x) {
        # Fit t-distribution with mean=0, sd=1 (standardized)
        x_scaled <- (x - mean(x))/sd(x)  # Standardize the data
        loglik <- function(df, x) {
          if (df <= 0) return(1e6)
          -sum(dt(x, df = df, log = TRUE))
        }
        opt <- optimize(loglik, c(1, 100), x = x_scaled)
        list(par = c(df = opt$minimum),
             value = opt$objective,
             convergence = 0,
             message = "Optimization complete")
      }
    )
  )
  
  if (!dist %in% names(dist_info)) {
    stop(paste("Unsupported 1-parameter distribution:", dist))
  }
  
  # Apply truncation
  x <- x[x >= lowertrunc & x <= uppertrunc]
  if (length(x) < 3) {
    warning("Insufficient data after truncation")
    return(create_failed_result_1param(dist, components, length(x)))
  }
  
  # Set default starting value if not provided
  if (is.null(start)) {
    start <- dist_info[[dist]]$start_fun(x)
    names(start) <- dist_info[[dist]]$params
  }
  
  # Set default bounds if not provided
  if (is.null(lower)) {
    lower <- 1e-8
    names(lower) <- dist_info[[dist]]$params
  }
  
  if (is.null(upper)) {
    upper <- Inf
    names(upper) <- dist_info[[dist]]$params
  }
  
  # Fit the distribution
  fit <- tryCatch({
    fit_fun <- dist_info[[dist]]$fit_fun
    result <- fit_fun(x)
    
    if (is.null(result$par) || any(is.na(result$par))) {
      stop("Fit returned NA parameters")
    }
    result
  }, error = function(e) {
    warning(paste("Fitting failed for", dist, ":", e$message))
    NULL
  })
  
  # Replace the output preparation with:
  if (!is.null(fit) && !is.null(fit$par)) {
    result <- update_combdist_output(fit, dist, components, length(x))
  } else {
    result <- create_failed_result_1param(dist, components, length(x))
  }
  
  if (nested) {
    result <- list(
      dist = result$dist,
      components = result$components,
      prior = list(result$prior),
      coefficients = list(result$coefficients),
      np = result$np,
      n = result$n,
      convergence = result$convergence,
      loglik = result$loglik,
      NLL = result$NLL,
      AIC = result$AIC,
      BIC = result$BIC,
      message = result$message
    )
  }
  
  return(result)
}

#################################################################################################################
# 2 parameter distributions
#################################################################################################################
combdist.mle_2param <- function(x, dist, start = NULL, lower = NULL, upper = NULL,  
                                components = 1, nested = FALSE, steps = 1,  
                                lowertrunc = 0, uppertrunc = Inf, ...) {
  
  # Supported distributions and their parameter names
  dist_info <- list(
    power = list(params = c("beta", "theta"), start = c(max(x)*1.1, 1)),
    fisk = list(params = c("beta", "theta"), start = c(median(x), 1)),
    lomax = list(params = c("beta", "theta"), start = c(median(x), 1)),
    invlomax = list(params = c("beta", "theta"), start = c(median(x), 1))
  )
  
  if (!dist %in% names(dist_info)) {
    stop(paste("Unsupported distribution:", dist))
  }
  
  # Apply truncation
  x <- x[x >= lowertrunc & x <= uppertrunc]
  if (length(x) < 2) {
    warning("Insufficient data after truncation")
    return(create_failed_result_2param(dist, components, length(x)))
  }
  
  # Set default starting values if not provided
  if (is.null(start)) {
    start <- dist_info[[dist]]$start
    names(start) <- dist_info[[dist]]$params
  }
  
  # Set default bounds if not provided
  if (is.null(lower)) {
    lower <- rep(1e-8, length(start))
    names(lower) <- dist_info[[dist]]$params
  }
  
  if (is.null(upper)) {
    upper <- rep(Inf, length(start))
    names(upper) <- dist_info[[dist]]$params
  }
  
  # Fit the distribution
  fit <- tryCatch({
    # Get the correct fitting function
    fit_fun <- get(paste0(dist, ".mle"))
    
    # Perform the fit
    result <- fit_fun(
      x = x,
      start = start,
      lower = lower,
      upper = upper,
      ...
    )
    
    # Validate the result
    if (is.null(result$par) || any(is.na(result$par))) {
      stop("Fit returned NA parameters")
    }
    result
  }, error = function(e) {
    warning(paste("Fitting failed for", dist, ":", e$message))
    NULL
  })
  
  # Prepare the output structure
  if (!is.null(fit) && !is.null(fit$par)) {
    result <- update_combdist_output(fit, dist, components, length(x))
  } else {
    result <- create_failed_result_2param(dist, components, length(x))
  }
  
  if (nested) {
    result <- list(
      dist = result$dist,
      components = result$components,
      prior = list(result$prior),
      coefficients = list(result$coefficients),
      np = result$np,
      n = result$n,
      convergence = result$convergence,
      loglik = result$loglik,
      NLL = result$NLL,
      AIC = result$AIC,
      BIC = result$BIC,
      message = result$message
    )
  }
  
  return(result)
}


#################################################################################################################
# 3-parameter distributions
#################################################################################################################

combdist.mle_3param <- function(x, dist, start = NULL, lower = NULL, upper = NULL,  
                                components = 1, nested = FALSE, steps = 1,  
                                lowertrunc = 0, uppertrunc = Inf, ...) {
  
  # Supported distributions and their parameter names
  dist_info <- list(
    b1 = list(params = c("b", "p", "q"), start = c(max(x)*1.1, 1, 1)),
    b2 = list(params = c("b", "p", "q"), start = c(median(x), 2, 2)),
    br12 = list(params = c("a", "b", "q"), start = c(1, median(x), 1)),
    br3 = list(params = c("a", "b", "p"), start = c(1, median(x), 1)),
    ib1 = list(params = c("b", "p", "q"), start = c(min(x)*0.9, 1, 1)),
    gg = list(params = c("b", "p", "q"), start = c(median(x), 1, 1)),
    ug = list(params = c("b", "d", "q"), start = c(1, 1, 1))
  )
  
  if (!dist %in% names(dist_info)) {
    stop(paste("Unsupported distribution:", dist))
  }
  
  # Apply truncation
  x <- x[x >= lowertrunc & x <= uppertrunc]
  if (length(x) < 3) {
    warning("Insufficient data after truncation")
    return(create_failed_result(dist, components, length(x)))
  }
  
  # Special case handling
  if (dist == "ug") {
    x <- x[x > 1]
    if (length(x) < 3) {
      warning("Insufficient data > 1 for UG distribution")
      return(create_failed_result(dist, components, length(x)))
    }
  }
  
  # Set default starting values if not provided
  if (is.null(start)) {
    start <- dist_info[[dist]]$start
    names(start) <- dist_info[[dist]]$params
  }
  
  # Set default bounds if not provided
  if (is.null(lower)) {
    lower <- rep(1e-8, 3)
    names(lower) <- dist_info[[dist]]$params
  }
  
  if (is.null(upper)) {
    upper <- rep(Inf, 3)
    names(upper) <- dist_info[[dist]]$params
  }
  
  # Fit the distribution
  fit <- tryCatch({
    # Get the correct fitting function
    fit_fun <- get(paste0(dist, ".mle"))
    
    # Perform the fit
    result <- fit_fun(
      x = x,
      start = start,
      lower = lower,
      upper = upper,
      ...
    )
    
    # Validate the result
    if (is.null(result$par) || any(is.na(result$par))) {
      stop("Fit returned NA parameters")
    }
    result
  }, error = function(e) {
    warning(paste("Fitting failed for", dist, ":", e$message))
    NULL
  })
  
  # Prepare the output structure
  if (!is.null(fit) && !is.null(fit$par)) {
    result <- update_combdist_output(fit, dist, components, length(x))
  } else {
    result <- create_failed_result_3param(dist, components, length(x))
  }
  
  if (nested) {
    result <- list(
      dist = result$dist,
      components = result$components,
      prior = list(result$prior),
      coefficients = list(result$coefficients),
      np = result$np,
      n = result$n,
      convergence = result$convergence,
      loglik = result$loglik,
      NLL = result$NLL,
      AIC = result$AIC,
      BIC = result$BIC,
      message = result$message
    )
  }
  
  return(result)
}


#################################################################################################################
# 4 and 5-parameter distributions
#################################################################################################################

combdist.mle_4_5_params <- function(x, dist, start = NULL, lower = NULL, upper = NULL,
                                    components = 1, nested = FALSE, steps = 1,
                                    lowertrunc = 0, uppertrunc = Inf, ...) {
  
  param_dists <- c("gb", "gb1", "gb2", "beta4")
  
  if (!(dist %in% param_dists)) {
    stop(paste("Unsupported 4 and 5-parameter distribution:", dist))
  }
  
  # Apply truncation if specified
  x <- x[x >= lowertrunc & x <= uppertrunc]
  if (length(x) < 3) {
    warning("Insufficient data after truncation")
    return(list(dist = dist, components = components, coefficients = NA, 
                np = 0, n = length(x), convergence = 1, loglik = NA))
  }
  
  output <- tryCatch({
    fit <- do.call(paste0(dist, ".mle"),
                   c(list(x = x, start = start, lower = lower, upper = upper), list(...)))
    
    # Ensure the output contains log-likelihood value
    if (!is.null(fit) && !is.null(fit$par)) {
      result <- update_combdist_output(fit, dist, components, length(x))
    } else {
      result <- create_failed_result_4_5_param(dist, components, length(x))
    }
    
    if (nested) {
      result <- list(
        dist = result$dist,
        components = result$components,
        prior = list(result$prior),
        coefficients = list(result$coefficients),
        np = result$np,
        n = result$n,
        convergence = result$convergence,
        loglik = result$loglik,
        NLL = result$NLL,
        AIC = result$AIC,
        BIC = result$BIC,
        message = result$message
      )
    }
    
    return(result)
  }, error = function(e) {
    warning("Fitting failed: ", conditionMessage(e))
    return(create_failed_result_1param(dist, components, length(x)))
  })
  
  return(output)
}


#################################################################################################################
# Updating existing one
#################################################################################################################

FLXMCdist1 <- function(dist = "lnorm", lowertrunc = NULL, uppertrunc = NULL, 
                       lower = c(-Inf, 1e-10), upper = c(Inf, Inf), ...) {
  
  # Handle special cases first
  if (dist == "norm") {
    return(FLXMRglm(family = "gaussian"))
  }
  
  if (dist == "lnormtrunc") {
    if (is.null(lowertrunc) || is.null(uppertrunc)) {
      stop("For lnormtrunc, both lowertrunc and uppertrunc must be specified")
    }
    return(FLXMClnormtrunc(lowertrunc = lowertrunc, uppertrunc = uppertrunc, 
                           lower = lower, upper = upper, ...))
  }
  
  if (dist == "frechet") {
    return(FLXMCfrechet(lower = lower, upper = upper, ...))
  }
  
  # Default case uses FLXMCnorm1 for other distributions
  FLXMCnorm1()
}

# FLXMCnorm1 implementation (needed as base case)
FLXMCnorm1 <- function(formula = .~., ...) {
  retval <- new("FLXMC", 
                weighted = TRUE,
                formula = formula,
                name = "FLXMCnorm1")
  
  retval@defineComponent <- function(para) {
    predict <- function(x, ...) {
      matrix(para$center, nrow = nrow(x), ncol = length(para$center), 
             byrow = TRUE)
    }
    
    logLik <- function(x, y) {
      dnorm(y, mean = para$center, sd = para$scale, log = TRUE)
    }
    
    new("FLXcomponent", 
        parameters = list(center = para$center, scale = para$scale),
        logLik = logLik, 
        predict = predict, 
        df = para$df)
  }
  
  retval@fit <- function(x, y, w, ...) {
    para <- cov.wt(cbind(y), wt = w)[c("center", "cov")]
    para$df <- ncol(x) + 2
    para$scale <- sqrt(para$cov)
    with(para, retval@defineComponent(para))
  }
  retval
}

combdist.mle_update <- function(
    x, dist, start = NULL, lower = NULL, upper = NULL,
    components = 1, nested = FALSE, steps = 1,
    lowertrunc = 0, uppertrunc = Inf, ...
) {
  create_failed_result <- function(dist, components, n, msg = "Fit failed") {
    list(dist = dist, components = components, prior = NA, coefficients = NA,
         np = 0, n = n, convergence = FALSE,
         loglik = NA, NLL = NA, AIC = NA, BIC = NA, message = msg)
  }
  
  raw_dist <- dist
  dist <- unlist(strsplit(dist, split = "_"))
  
  # -----------------------------
  # 1) Composite distributions
  # -----------------------------
  if (length(dist) > 1) {
    output <- tryCatch(
      composite.mle(x = x, dist = dist, start = start, ...),
      error = function(e) return(create_failed_result(paste(dist, collapse = "_"),
                                                      length(dist), length(x), e$message))
    )
    if (is.list(output) && is.null(output$coefficients)) return(output)
    
    dist <- paste0(dist, collapse = "_")
    components <- length(strsplit(dist, "_")[[1]])
    coeff <- output$coefficients
    np <- length(output$coefficients)
    prior <- 1
    n_used <- length(x)
    loglik <- tryCatch(
      sum(dcombdist(dist = dist, x = x, prior = prior, coeff = coeff, log = TRUE)),
      error = function(e) NA_real_
    )
    convergence <- output$convergence
    message <- "Success"
    
    # -----------------------------
    # 2) Special Pareto / hybrids
    # -----------------------------
  } else if (dist %in% c("pareto", "invpareto",
                         "leftparetolognormal", "doubleparetolognormal",
                         "rightparetolognormal")) {
    dist <- dist[1]
    output <- tryCatch(
      do.call(paste0(dist, ".mle"),
              c(list(x = x), lower = lower, upper = upper, ...)),
      error = function(e) return(create_failed_result(dist,
                                                      ifelse(dist %in% c("leftparetolognormal","rightparetolognormal"), 2,
                                                             ifelse(dist %in% c("pareto","invpareto"), 1, 3)),
                                                      length(x), e$message))
    )
    if (is.list(output) && is.null(output$coefficients)) return(output)
    
    components <- if (dist %in% c("leftparetolognormal","rightparetolognormal")) 2
    else if (dist %in% c("pareto","invpareto")) 1
    else 3
    coeff <- output$coefficients
    np <- length(output$coefficients)
    prior <- 1
    n_used <- length(x)
    loglik <- tryCatch(
      sum(dcombdist(dist = dist, x = x, prior = prior, coeff = coeff, log = TRUE)),
      error = function(e) NA_real_
    )
    convergence <- output$convergence
    message <- "Success"
    
    # -----------------------------
    # 3) Base-R handled distributions
    # -----------------------------
  } else {
    dist <- dist[1]
    
    # truncation
    if (!is.null(lowertrunc)) {
      x <- x[(x >= lowertrunc) & (x <= uppertrunc)]
    }
    if (length(x) < 3) {
      return(create_failed_result(dist, components, length(x),
                                  "Insufficient data after truncation"))
    }
    
    prior <- 1
    n_used <- length(x)
    coeff <- NULL
    loglik <- NA
    np <- NA
    
    fit_result <- list(convergence = 0, message = NULL)  # local fit result
    
    if (dist == "lnorm") {
      mu <- mean(log(x)); sigma <- sd(log(x))
      coeff <- rbind(meanlog = mu, sdlog = sigma)
      loglik <- sum(dlnorm(x, meanlog = mu, sdlog = sigma, log = TRUE))
      np <- 2
    } else if (dist == "gamma") {
      fit <- tryCatch(MASS::fitdistr(x, densfun = "gamma"),
                      error = function(e) NULL)
      if (is.null(fit)) return(create_failed_result(dist, components, length(x), "Gamma fit failed"))
      shape <- fit$estimate["shape"]; rate <- fit$estimate["rate"]
      coeff <- rbind(shape = shape, rate = rate)
      loglik <- sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
      np <- 2
    } else if (dist == "weibull") {
      fit <- tryCatch(MASS::fitdistr(x, densfun = "weibull"),
                      error = function(e) NULL)
      if (is.null(fit)) return(create_failed_result(dist, components, length(x), "Weibull fit failed"))
      shape <- fit$estimate["shape"]; scale <- fit$estimate["scale"]
      coeff <- rbind(shape = shape, scale = scale)
      loglik <- sum(dweibull(x, shape = shape, scale = scale, log = TRUE))
      np <- 2
    } else if (dist == "exp") {
      rate <- 1 / mean(x)
      coeff <- rbind(rate = rate)
      loglik <- sum(dexp(x, rate = rate, log = TRUE))
      np <- 1
    } else if (dist == "norm") {
      mu <- mean(x); sigma <- sd(x)
      coeff <- rbind(mean = mu, sd = sigma)
      loglik <- sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
      np <- 2
    } else if (dist == "frechet") {
      # keep Frechet via FLXMCfrechet
      output <- tryCatch(
        flexmix::stepFlexmix(
          x ~ 1,
          k = components,
          model = FLXMCfrechet(),
          control = list(minprior = 0, iter.max = 10000),
          nrep = steps
        ),
        error = function(e) NULL
      )
      if (is.null(output)) {
        return(create_failed_result(dist, components, length(x), "Frechet fit failed"))
      }
      prior <- output@prior
      coeff <- flexmix::parameters(output)
      np <- nrow(coeff)
      loglik <- tryCatch(
        sum(dcombdist(dist = dist, x = x, prior = prior, coeff = coeff, log = TRUE)),
        error = function(e) NA_real_
      )
      convergence <- output@converged
      message <- "Success"
    } else {
      return(create_failed_result(dist, components, length(x),
                                  "Unsupported distribution in base-R version"))
    }
    
    # Convergence + message logic (only for base-R dists)
    if (dist %in% c("lnorm","gamma","weibull","exp","norm")) {
      convergence <- ifelse(is.null(fit_result$convergence),
                            FALSE,
                            fit_result$convergence == 0)
      message <- ifelse(is.null(fit_result$message),
                        ifelse(is.null(fit_result$convergence) || fit_result$convergence != 0,
                               "Convergence issue", "Success"),
                        fit_result$message)
    }
  }
  
  # -----------------------------
  # Info criteria
  # -----------------------------
  NLL <- ifelse(is.na(loglik), NA, -loglik)
  AIC <- ifelse(is.na(loglik), NA, 2 * np - 2 * loglik)
  BIC <- ifelse(is.na(loglik), NA, log(n_used) * np - 2 * loglik)
  
  # -----------------------------
  # Return
  # -----------------------------
  return_list <- list(
    dist = if (length(dist) == 1) dist else paste(dist, collapse = "_"),
    components = components,
    prior = prior,
    coefficients = coeff,
    np = np,
    n = n_used,
    convergence = convergence,
    loglik = loglik,
    NLL = NLL,
    AIC = AIC,
    BIC = BIC,
    message = message
  )
  
  if (nested) {
    return_list <- list(
      dist = return_list$dist,
      components = components,
      prior = list(prior),
      coefficients = list(coeff),
      np = np,
      n = n_used,
      convergence = convergence,
      loglik = loglik,
      NLL = NLL,
      AIC = AIC,
      BIC = BIC,
      message = message
    )
  }
  
  return(return_list)
}

#################################################################################################################
# KS Test
#################################################################################################################

# ========== 1-parameter families ==========

# Uniform(0, beta)
punif1 <- function(q, beta) {
  res <- q / beta
  res[q < 0] <- 0
  res[q > beta] <- 1
  pmin(pmax(res, 0), 1)
}

# Rayleigh(sigma)
prayleigh <- function(q, sigma) 1 - exp(-q^2 / (2 * sigma^2))

# Log-logistic(alpha, beta)
ploglogistic <- function(q, alpha, beta = 1) {
  1 / (1 + (beta / q)^alpha)
}

# Normal(0, sd)
pnormal <- function(q, sd) stats::pnorm(q, mean = 0, sd = sd)

# Student-t(df) (standardized)
pt_1param <- function(q, df) stats::pt(q, df = df)


# ========== 2-parameter families ==========

# Power(beta, theta), support [0,beta]
ppower <- function(q, beta, theta) {
  res <- (q / beta)^theta
  res[q < 0] <- 0
  res[q > beta] <- 1
  pmin(pmax(res, 0), 1)
}

# Fisk(beta, theta) ≡ Log-logistic
pfisk <- function(q, beta, theta) {
  1 / (1 + (beta / q)^theta)
}

# Lomax(beta, theta)
plomax <- function(q, beta, theta) {
  1 - (1 + q / beta)^(-theta)
}

# Inverse Lomax(beta, theta), support (0,∞)
pinvlomax <- function(q, beta, theta) {
  res <- (q / (beta + q))^theta
  res[q <= 0] <- 0
  res
}


# --------- tiny helper used by all numeric CDFs (no external deps) ----------
.numeric_cdf_from_pdf <- function(x, pdf_fun, args, lower, upper = Inf, n = 2048L) {
  # Build an integration grid that covers the query x
  xmin <- max(lower, suppressWarnings(min(x[x > 0], na.rm = TRUE)))
  if (!is.finite(xmin)) xmin <- lower
  xmax_target <- suppressWarnings(max(x, na.rm = TRUE))
  xmax <- if (is.finite(upper)) upper else xmax_target
  if (!is.finite(xmax) || xmax <= xmin) xmax <- xmin * 10
  
  # Log grid if wide range, else linear
  gx <- if (xmin > 0 && is.finite(xmax) && (xmax / xmin) > 50) {
    exp(seq(log(xmin), log(xmax), length.out = n))
  } else {
    seq(max(xmin, .Machine$double.eps), xmax, length.out = n)
  }
  
  dens <- do.call(pdf_fun, c(list(gx), args))
  dens[!is.finite(dens)] <- 0
  dens[dens < 0] <- 0
  
  # Trapezoidal cum integral, normalize to 1
  dx   <- diff(gx)
  avg  <- (dens[-1] + dens[-length(dens)]) / 2
  Fg   <- c(0, cumsum(dx * avg))
  tot  <- Fg[length(Fg)]
  if (!is.finite(tot) || tot <= 0) return(rep(NA_real_, length(x)))
  Fg <- Fg / tot
  
  # Interpolate to the requested x; clamp to [0,1]
  approx_vals <- stats::approx(gx, Fg, xout = x, yleft = 0, yright = 1, ties = "ordered")$y
  pmin(pmax(approx_vals, 0), 1)
}

# GB1: dgb1(y, a, b, p, q), support: 0 < y < b
pgb1 <- function(x, a, b, p, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, dgb1, list(a = a, b = b, p = p, q = q), lower = 0, upper = b, n = n)
}

# GB2: dgb2(y, a, b, p, q), support: y > 0
pgb2 <- function(x, a, b, p, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, dgb2, list(a = a, b = b, p = p, q = q), lower = 0, upper = Inf, n = n)
}

# Generalized Beta (5-param): dgb(y, a, b, c, p, q)
# Effective upper support is b / (1 - c)^(1/a) when 0 <= c < 1
pgb <- function(x, a, b, c, p, q, n = 2048L) {
  up <- if (is.finite(c) && c >= 0 && c < 1) b / (1 - c)^(1 / a) else Inf
  .numeric_cdf_from_pdf(x, dgb, list(a = a, b = b, c = c, p = p, q = q), lower = 0, upper = up, n = n)
}

# Beta4: dbeta4(y, b, c, p, q), support: 0 < y < b/(1-c), 0 <= c < 1
pbeta4 <- function(x, b, c, p, q, n = 2048L) {
  up <- b / (1 - c)
  .numeric_cdf_from_pdf(x, dbeta4, list(b = b, c = c, p = p, q = q), lower = 0, upper = up, n = n)
}
# BR12: dbr12(y, a, b, q), support: y > 0
pbr12 <- function(x, a, b, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, dbr12, list(a = a, b = b, q = q), lower = 0, upper = Inf, n = n)
}

# BR3: dbr3(y, a, b, p), support: y > 0
pbr3 <- function(x, a, b, p, n = 2048L) {
  .numeric_cdf_from_pdf(x, dbr3, list(a = a, b = b, p = p), lower = 0, upper = Inf, n = n)
}

# B1: db1(y, b, p, q), support: 0 < y < b
pb1 <- function(x, b, p, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, db1, list(b = b, p = p, q = q), lower = 0, upper = b, n = n)
}

# B2: db2(y, b, p, q), support: y > 0
pb2 <- function(x, b, p, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, db2, list(b = b, p = p, q = q), lower = 0, upper = Inf, n = n)
}

# IB1: dib1(y, b, p, q), support: y > b
pib1 <- function(x, b, p, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, dib1, list(b = b, p = p, q = q), lower = b, upper = Inf, n = n)
}

# UG (Unit Gamma): dug(y, b, d, q), support: y > 1
pug <- function(x, b, d, q, n = 2048L) {
  .numeric_cdf_from_pdf(x, dug, list(b = b, d = d, q = q), lower = 1, upper = Inf, n = n)
}

# Generalized Gamma CDF (matches your dgg with params b, p, q; y > 0, q > 0)
pgg <- function(x, b, p, q) {
  u <- (x / b)^q
  res <- stats::pgamma(u, shape = p, rate = 1)
  res[x <= 0] <- 0
  res[!is.finite(res)] <- 0
  res
}

.map_cdf_name <- function(nm) {
  nm <- tolower(nm)
  switch(nm,
         "loglogis" = "loglogistic",   # alias used elsewhere
         "t"        = "t_1param",      # your CDF is pt_1param()
         # keep canonical names as-is:
         "fisk"     = "fisk",
         "lomax"    = "lomax",
         "invlomax" = "invlomax",
         "power"    = "power",
         "gb1"      = "gb1",
         "gb2"      = "gb2",
         "beta4"    = "beta4",
         "b1"       = "b1",
         "b2"       = "b2",
         "br12"     = "br12",
         "br3"      = "br3",
         "ib1"      = "ib1",
         "ug"       = "ug",
         nm)
}

calculate_ks_test <- function(dist_name, params, data) {
  actual <- .map_cdf_name(dist_name)
  cdf_fun <- tryCatch(get(paste0("p", actual), mode = "function"),
                      error = function(e) NULL)
  if (is.null(cdf_fun)) {
    warning(paste("No CDF function found for distribution:", dist_name))
    return(list(statistic = NA, p.value = NA))
  }
  args <- c(list(data), as.list(params))
  suppressWarnings({
    tryCatch({ do.call(ks.test, args) },
             error = function(e) {
               warning(paste("KS test failed for", dist_name, ":", e$message))
               list(statistic = NA, p.value = NA)
             })
  })
}


#################################################################################################################
# Printing output and running the distributions
#################################################################################################################

process_distribution_fit <- function(y, j, x, fit_function) {
  # Convert input to data frame if it's a tibble
  x <- as.data.frame(x)
  
  # Initialize output data frame
  output <- data.frame()
  
  # Get unique combinations of varind and year
  groups <- unique(x[, c("varind", "year")])
  
  for(i in 1:nrow(groups)) {
    current_group <- groups[i, ]
    group_data <- x[x$varind == current_group$varind & x$year == current_group$year, ]
    x_values <- group_data$value
    
    # Fit distribution
    fit_result <- fit_function(
      x = x_values,
      components = j,
      dist = y,
      steps = 1,
      nested = TRUE
    )
    
    # Calculate KS test
    ks_result <- if (!is.null(fit_result$coefficients) && length(fit_result$coefficients) > 0) {
      calculate_ks_test(y, fit_result$coefficients, x_values)
    } else {
      list(statistic = NA, p.value = NA)
    }
    
    # Create output row
    result_row <- data.frame(
      varind = current_group$varind,
      year = current_group$year,
      dist = fit_result$dist,
      components = fit_result$components,
      np = fit_result$np,
      n = fit_result$n,
      convergence = fit_result$convergence,
      loglik = fit_result$loglik,
      NLL = fit_result$NLL,
      AIC = fit_result$AIC,
      BIC = fit_result$BIC,
      ks_statistic = ks_result$statistic,
      ks_pvalue = ks_result$p.value,
      stringsAsFactors = FALSE
    )
    
    # Add coefficients as a list column
    result_row$coefficients <- list(fit_result$coefficients)
    result_row$prior <- list(fit_result$prior)
    result_row$message <- fit_result$message
    
    # Combine with output
    output <- rbind(output, result_row)
  }
  
  
  # Make sure output folder exists
  output_dir <- file.path("empirics", "output", "domsales")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  var_name <- paste0(y, j)
  assign(var_name, output, envir = .GlobalEnv)
  save(list = var_name, file = file.path(output_dir, paste0(var_name, ".Rdata")))
  
  return(output)
}