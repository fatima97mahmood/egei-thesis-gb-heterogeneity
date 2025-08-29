# Augments the FLXM components from the flexmix package with the truncated lognormal distribution

FLXMClnormtrunc <- function(formula = . ~ ., method = "Nelder-Mead", warn = -1, lowertrunc, uppertrunc, lower = c(-Inf, 1e-10), upper = c(Inf, Inf), ...) {
  z <- new("FLXMC",
    weighted = TRUE, formula = formula,
    name = "model-based truncated Lognormal clustering",
    dist = "lnormtrunc"
  )
  z@preproc.y <- prepoc.y.pos.1

  z@defineComponent <- function(para) {
    predict <- function(x, ...) {
      matrix(para$meanlog,
        nrow = nrow(x), ncol = length(para$meanlog),
        byrow = TRUE
      )
    }


    logLik <- function(x, y) {
      (dtruncdist(y, dist = "lnormtrunc", coeff = list(meanlog = para$meanlog, sdlog = para$sdlog), lowertrunc = lowertrunc, uppertrunc = uppertrunc, log = TRUE))
    }

    new("FLXcomponent",
      parameters = list(meanlog = para$meanlog, sdlog = para$sdlog), df = para$df, logLik = logLik,
      predict = predict
    )
  }


  z@fit <- function(x, y, w, component, ...) {
    if (!length(component)) {
      start <- list(mean(log(y)), sd(log(y)))
    } else {
      start <- unname(unlist(component))
    }


    f <- function(coeff, w, lowertrunc, uppertrunc) {
      out <- -sum(dtruncdist(y, dist = "lnormtrunc", coeff = list(meanlog = coeff[1], sdlog = coeff[2]), lowertrunc = lowertrunc, uppertrunc = uppertrunc, log = TRUE) * w)
      if (is.infinite(out)) {
        out <- 1e20
      }
      return(out)
    }
    oop <- options(warn = warn)
    on.exit(options(oop))
    # Adapt optimization method to allow for upper and lower bounds
    # Original code: parms <- optim(start, f, method = "Nelder-Mead")$par
    if (any(start < lower)) {
      start[start < lower] <- lower[start < lower]
    }
    if (any(start > upper)) {
      start[start > upper] <- upper[start > upper]
    }

    parms <- suppressWarnings(nlminb(start, f, lower = lower, upper = upper, w = w, lowertrunc = lowertrunc, uppertrunc = uppertrunc, control = list(maxit = 1e5)))$par
    z@defineComponent(list(meanlog = parms[1], sdlog = parms[2], df = 2))
  }
  z
}

# Augments the FLXM components from the flexmix package with the frechet distribution

FLXMCfrechet <- function(formula = . ~ ., method = "Nelder-Mead", warn = -1, lower = c(1e-10, 1e-10), upper = c(Inf, Inf), ...) {
  z <- new("FLXMC",
    weighted = TRUE, formula = formula,
    name = "model-based frechet clustering",
    dist = "frechet"
  )
  z@preproc.y <- prepoc.y.pos.1

  z@defineComponent <- function(para) {
    predict <- function(x, ...) {
      matrix(para$shape,
        nrow = nrow(x), ncol = length(para$shape),
        byrow = TRUE
      )
    }

    logLik <- function(x, y) {
      (dfrechet(y, shape = para$shape, scale = para$scale, log = TRUE))
    }

    new("FLXcomponent",
      parameters = list(shape = para$shape, scale = para$scale), df = para$df, logLik = logLik,
      predict = predict
    )
  }

  z@fit <- function(x, y, w, component, lower = c(1e-10, 1e-10), upper = c(Inf, Inf)) {
    if (!length(component)) {
      start <- c(0.5, 0.1)
    } else {
      start <- unname(unlist(component))
    }

    f <- function(parms) {
      out <- -sum(dfrechet(y, shape = parms[1], scale = parms[2], log = TRUE) * w)
      if (is.infinite(out)) {
        out <- 1e20
      }
      return(out)
    }
    oop <- options(warn = warn)
    on.exit(options(oop))
    # Adapt optimization method to allow for upper and lower bounds
    # Original code: parms <- optim(start, f, method = "Nelder-Mead")$par
    if (any(start < lower)) {
      start[start < lower] <- lower[start < lower]
    }
    if (any(start > upper)) {
      start[start > upper] <- upper[start > upper]
    }

    # parms = frechet.mle(x=y,weights=w,lower = lower,upper=upper,start=start)$coefficients

    parms <- suppressWarnings(nlminb(start, f, lower = lower, upper = upper, control = list(maxit = 1e5)))$par
    z@defineComponent(list(shape = parms[1], scale = parms[2], df = 2))
  }
  z
}


# FLXMCdist1: internal function

prepoc.y.pos.1 <- function(x) {
  if (ncol(x) > 1) {
    stop("for the inverse gaussian family y must be univariate")
  }
  if (any(x < 0)) {
    stop("values must be >= 0")
  }
  x
}