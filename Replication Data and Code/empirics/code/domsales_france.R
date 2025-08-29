# SETUP -------------------------------------------------------------------
#Load libraries
setwd("/Users/fatimamahmood/Documents/Files/EGEI/Thesis/Replication Data and Code")
library(tidyverse)
library(here)
library(xtable) # Output latex tables
library(latex2exp) # Add latex to legends of figures
library(ggpubr) # Combine multiple figures to one
library(distributionsrd)
library(dplyr)
library(purrr)
library(ggplot2)
library(scales)

#Load source code specific to this analysis
source(here("empirics","code","toolbox.R"))
source(here("empirics","code","melitz.R"))
source(here("empirics","code","toolbox_melitz.R"))
source(here("empirics","code","nmad.R"))
source(here("empirics","code","FLXMC_xtra.R"))
source(here("empirics","code","axillary_code.R"))

#Styles
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#c9bd22", "#0072B2", "#D55E00", "#CC79A7")
linetypes = c("solid", "22", "42", "44", "13", "1343", "73", "2262", "12223242","F282", "F4448444", "224282F2", "F1")

#Variable initialization
maxprior = 5 # Maximum number of distribution components
distselect_main = c("lnorm","pareto","doubleparetolognormal",
                    "leftparetolognormal","rightparetolognormal",
                    "invpareto_lnorm_pareto","lnorm_pareto","invpareto_lnorm") # Distributions used in main analysis
sigma = 4 #Elasticity of substitution

# Code based on available data -------------------------------------------------------------
# DATA --------------------------------------------------------------------
load(here("empirics","input","data","orbis_france.Rdata"))

data <- new.env()
load(here("empirics", "input", "data", "orbis_france.Rdata"), envir = data)
domsales2012=x$domsales
domsales2012 = domsales2012[domsales2012>0]
domsales2012 = na.omit(domsales2012)
domsales2012 = domsales2012/mean(domsales2012)

xy = data.frame(var="domsales",ind="Total",varind="domsales.Total",year=2012,value=domsales2012)

#################################################################################################################
# Getting the parameters, AIC/BIC and KS values for all distributions 
#################################################################################################################

## existing distributions
dists <- c("lnorm", "pareto", "invpareto", "gamma", "weibull", "frechet", "exp", "leftparetolognormal",
           "doubleparetolognormal", "rightparetolognormal")
# Run for 1 component (change 1:1 to 1:n if you want more components)
lapply(dists, function(y) {
  lapply(1:1, function(j) {
    process_distribution_fit(y, j, xy, combdist.mle_update)
  })
})

# Supported 1-parameter distributions
one_param_dists <- c("unif1", "rayleigh", "chisq", 
                     "loglogistic", "normal", "t")

# Loop over each 1-parameter distribution 
lapply(one_param_dists, function(y) {
  lapply(1:1, function(j) {
    process_distribution_fit(y, j, xy, combdist.mle_1param)
  })
})

# Supported 2-parameter distributions
two_param_dists <- c("power", "fisk", "lomax", "invlomax")

# Loop over each 2-parameter distribution 
lapply(two_param_dists, function(y) {
  lapply(1:1, function(j) {
    process_distribution_fit(y, j, xy, combdist.mle_2param)
  })
})

# New and supported 3-parameter distributions
three_param_dists <- c("b1", "b2", "br12", "br3", "ib1", "gg", "ug")

# Loop over each 3-parameter distribution 
lapply(three_param_dists, function(y) {
  lapply(1:1, function(j) {
    process_distribution_fit(y, j, xy, combdist.mle_3param)
  })
})

# New supported 4 and 5 distributions
new_dists <- c("gb", "gb1", "gb2", "beta4")

# Loop over 4 and 5 distribution 
lapply(new_dists, function(y) {
  lapply(1:1, function(j) {
    process_distribution_fit(y, j, xy, combdist.mle_4_5_params)
  })
})

#################################################################################################################
# Outputting the results for new distributions into a table for visual inspection
#################################################################################################################

# List of all distribution fitting functions
distributions <- list(
  unif1.mle,
  rayleigh.mle,
  t1.mle,
  loglogis.mle,
  power.mle,
  fisk.mle,
  lomax.mle,
  invlomax.mle,
  ib1.mle,
  ug.mle,
  b1.mle,
  gg.mle,
  b2.mle,
  br12.mle,
  br3.mle,
  gb1.mle,
  gb2.mle,
  beta4.mle,
  gb.mle
)

# Names for each distribution
dist_names <- c(
  "unif1",
  "rayleigh",
  "t1",
  "loglogis",
  "power",
  "fisk",
  "lomax",
  "invlomax",
  "ib1",
  "ug",
  "b1",
  "gg",
  "b2",
  "br12",
  "br3",
  "gb1",
  "gb2",
  "beta4",
  "gb"
)

# Function to extract coefficients
extract_coefs <- function(fit) {
  if (is.list(fit)) {
    if ("par" %in% names(fit)) {
      return(as.numeric(fit$par))
    } else if ("sigma" %in% names(fit)) {
      return(as.numeric(fit$sigma))
    } else if ("beta" %in% names(fit)) {
      return(as.numeric(fit$beta))
    }
  }
  return(NA)
}

# Fit all distributions and store coefficients
coef_list <- list()

for (i in seq_along(distributions)) {
  tryCatch({
    fit <- distributions[[i]](domsales2012)
    coefs <- extract_coefs(fit)
    coef_list[[dist_names[i]]] <- coefs
  }, error = function(e) {
    coef_list[[dist_names[i]]] <- NA
  })
}

# Find maximum number of parameters
max_params <- max(sapply(coef_list, function(x) {
  if (any(is.na(x))) return(0)
  length(x)
}))

# Create results matrix
results <- matrix(NA, nrow = length(dist_names), ncol = max_params)
rownames(results) <- dist_names

# Fill the matrix
for (i in seq_along(coef_list)) {
  if (!all(is.na(coef_list[[i]])) && !is.null(coef_list[[i]])) {
    if (length(coef_list[[i]]) > 0) {
      results[i, 1:length(coef_list[[i]])] <- coef_list[[i]]
    }
  }
}

# Convert to data frame
results_df <- as.data.frame(results)
colnames(results_df) <- paste0("Param", 1:max_params)

# Add distribution names
results_df$Distribution <- rownames(results_df)

# Reorder columns
results_df <- results_df[, c("Distribution", paste0("Param", 1:max_params))]

# Custom formatting function
format_coefficients <- function(x) {
  ifelse(is.na(x), NA, 
         format(round(x, 6), nsmall = 6, scientific = FALSE))
}

# Apply to all parameter columns
for (col in grep("Param", names(results_df))) {
  results_df[[col]] <- format_coefficients(results_df[[col]])
}

#################################################################################################################
# Ranking all distributions based on log-likelihoods, AIC, and BIC
#################################################################################################################

# Define function to identify unreasonable distributions
identify_unreasonable_distributions <- function(results_df) {
  # Thresholds for each distribution's parameters
  reasonable_thresholds <- list(
    "unif1" = c(1000, NA, NA, NA, NA),
    "rayleigh" = c(1000, NA, NA, NA, NA),
    "t" = c(NA, NA, NA, NA, NA),
    "loglogistic" = c(NA, NA, NA, NA, NA),
    "power" = c(2, 2, NA, NA, NA),
    "fisk" = c(1000, 2, NA, NA, NA),
    "lomax" = c(1, 2, NA, NA, NA),
    "invlomax" = c(1e6, 1e6, NA, NA, NA),
    "ib1" = c(500, 5, 5, NA, NA),
    "ug" = c(1e3, 5, 1e3, NA, NA),
    "b1" = c(1000, 5, 500, NA, NA),
    "gg" = c(1, 5, 1, NA, NA),
    "b2" = c(1, 5, 5, NA, NA),
    "br12" = c(5, 1, 5, NA, NA),
    "br3" = c(5, 1, 5, NA, NA),
    "gb1" = c(1, 1, 5, 5, NA),
    "gb2" = c(1, 1, 20, 10, NA),
    "beta4" = c(1, 1, 5, 5, NA),
    "gb" = c(1, 1, 1, 5, 5)
  )
  
  # Initialize vector to store unreasonable distribution names
  unreasonable_dists <- character(0)
  
  # Check each distribution in results_df
  for (dist in rownames(results_df)) {
    # Skip if distribution not in our threshold list
    if (!dist %in% names(reasonable_thresholds)) {
      unreasonable_dists <- c(unreasonable_dists, dist)
      next
    }
    
    params <- as.numeric(results_df[dist, grep("Param", colnames(results_df))])
    thresholds <- reasonable_thresholds[[dist]]
    
    # Check all parameters against thresholds
    for (i in seq_along(params)) {
      if (!is.na(params[i]) && !is.na(thresholds[i]) && abs(params[i]) > thresholds[i]) {
        unreasonable_dists <- c(unreasonable_dists, dist)
        break
      }
    }
  }
  
  return(unique(unreasonable_dists))
}

# Identify unreasonable distributions from results_df
unreasonable_dists <- identify_unreasonable_distributions(results_df)

#Load all saved distribution fits with robust error handling
all_dist_files <- list.files(here("empirics", "output", "fit_cities"),
                             pattern = "\\.Rdata$", full.names = TRUE)

all_dists <- list()
for (file in all_dist_files) {
  load(file)
  dist_name <- gsub("\\.Rdata$", "", basename(file))
  dist_obj <- get(dist_name)
  
  # Skip if the object is empty or doesn't have required components
  if (length(dist_obj) == 0 || 
      !all(c("dist", "np", "n", "loglik", "convergence", "ks_statistic","ks_pvalue") %in% names(dist_obj))) {
    warning(paste("Skipping", dist_name, "- missing required components"))
    next
  }
  
  # Skip if the fit didn't converge (convergence is FALSE)
  if (!dist_obj$convergence[1]) {
    warning(paste("Skipping", dist_name, "- did not converge"))
    next
  }
  
  # Create data frame only if all components have length > 0 and converged
  if (length(dist_obj$dist) > 0 && 
      length(dist_obj$np) > 0 && 
      length(dist_obj$n) > 0 && 
      length(dist_obj$loglik) > 0 &&
      dist_obj$convergence[1]) {
    all_dists[[dist_name]] <- data.frame(
      dist = dist_obj$dist[1],
      np = dist_obj$np[1],
      n = dist_obj$n[1],
      loglik = dist_obj$loglik[1],
      convergence = dist_obj$convergence[1],
      ks_statistic = dist_obj$ks_statistic[1],
      ks_pvalue = dist_obj$ks_pvalue[1],
      stringsAsFactors = FALSE
    )
  } else {
    warning(paste("Skipping", dist_name, "- empty components or did not converge"))
  }
}

# Check if we loaded any distributions
if (length(all_dists) == 0) {
  stop("No valid distributions were loaded")
}

#Combine and calculate AIC/BIC
loglike_all <- do.call(rbind, all_dists) %>%
  mutate(
    NLL = -loglik,
    AIC = 2 * np - 2 * loglik,
    BIC = log(n) * np - 2 * loglik
  )

#Filter out unreasonable distributions
loglike_filtered <- loglike_all %>%
  filter(!dist %in% unreasonable_dists) %>%
  arrange(NLL)

#Rank and annotate BIC differences
loglike_ranked <- loglike_filtered %>%
  group_by(dist) %>%
  mutate(
    id = ifelse(rank(BIC) == 1, 1, 0),
    BIC.base = mean(ifelse(id == 1, BIC, NA), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    ks_statistic=ks_statistic, 
    ks_pvalue=ks_pvalue,
    BIC.diff = abs(BIC.base - BIC),
    BIC.sig = as.numeric(cut(BIC.diff, c(0, 2, 6, 10, Inf))),
    NLL = rank(NLL),
    AIC = rank(AIC),
    BIC = case_when(
      BIC.sig == 4 ~ paste0(rank(BIC), "{\\tiny$^{+++}$}"),
      BIC.sig == 3 ~ paste0(rank(BIC), "{\\tiny$^{++}$}"),
      BIC.sig == 2 ~ paste0(rank(BIC), "{\\tiny$^{+}$}"),
      TRUE ~ paste0(rank(BIC))
    )
  ) %>%
  arrange(NLL) %>%
  select(dist, par = np, loglik, NLL, AIC, BIC, ks_statistic, ks_pvalue)

#Save results
save(loglike_ranked, file = here("empirics", "output", "eval_cities", "loglike_ranked.Rdata"))

# Latex Table
create_summary_table <- function(loglike_ranked) {
  # Merge the ranked results with parameter estimates
  summary_table <- loglike_ranked %>%
    left_join(
      results_df %>% 
        mutate(dist = Distribution) %>% 
        select(-Distribution),
      by = "dist"
    ) %>%
    arrange(NLL) %>%
    select(dist, par, loglik, NLL, AIC, BIC, ks_statistic, ks_pvalue)
  
  # Format for LaTeX output
  print(
    xtable(summary_table, 
           caption = "Goodness-of-fit measures for different distributions",
           label = "tab:dist_fits"),
    include.rownames = FALSE,
    sanitize.text.function = function(x) x,
    caption.placement = "top"
  )
}

# Call the function
create_summary_table(loglike_ranked)

#################################################################################################################
# Plotting Distributions 
#################################################################################################################
output_dir <- here::here("empirics", "output", "distribution_fits")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- helpers ----------
si_lab <- function() scales::label_number(scale_cut = scales::cut_si(" "))

# numeric quantile from a CDF (fallback)
q_from_cdf <- function(u, cdf_fun, lower, upper, tol = 1e-10, maxit = 100) {
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  sapply(u, function(p) uniroot(function(z) cdf_fun(z) - p,
                                interval = c(lower, upper),
                                tol = tol, maxiter = maxit)$root)
}

# generic plot maker
.make_plots_ccdf_qq <- function(x, cdf_fun, q_fun, title_prefix,
                                trim = c(0.001, 0.9995),
                                n_grid = 400, n_pp = 2000) {
  x <- x[is.finite(x) & x > 0]
  ec <- ecdf(x)
  
  # CCDF grid
  xg <- exp(seq(log(quantile(x, trim[1])), log(quantile(x, trim[2])), length.out = n_grid))
  emp_ccdf <- pmax(1 - ec(xg), .Machine$double.eps)
  fit_ccdf <- pmax(1 - cdf_fun(xg), .Machine$double.eps)
  
  # COLORS + LEGEND (added)
  p_ccdf <- ggplot(data.frame(x=xg, emp=emp_ccdf, fit=fit_ccdf), aes(x)) +
    geom_point(aes(y=emp, colour="Empirical"), size=0.7, alpha=0.5) +
    geom_line(aes(y=fit, colour="Fitted"), linewidth=0.7) +
    scale_color_manual(values=c("Empirical"="#000000","Fitted"="#E69F00")) +
    scale_x_log10(labels = si_lab()) +
    scale_y_log10(labels = si_lab()) +
    labs(title=paste0(title_prefix, " CCDF (log–log)"),
         x="Value (log scale)", y="P(Y > y)", colour="") +
    theme(legend.position="bottom")
  
  # Q–Q
  pp    <- ppoints(n_pp)
  q_emp <- quantile(x, probs = pp, names = FALSE)
  q_fit <- q_fun(pp)
  
  # Add legend by mapping points + a 45° reference line as a geom_line with linetype legend
  rng <- range(c(q_emp, q_fit), finite = TRUE)
  line_df <- data.frame(Fitted = rng, Empirical = rng)
  
  p_qq <- ggplot(data.frame(Fitted=q_fit, Empirical=q_emp), aes(Fitted, Empirical)) +
    geom_point(aes(colour="Empirical vs Fitted"), alpha=0.35, size=0.6) +
    geom_line(data=line_df, aes(Fitted, Empirical, linetype="45° line"), colour="grey40") +
    scale_color_manual(values=c("Empirical vs Fitted"="#E69F00")) +
    scale_linetype_manual(values=c("45° line"="dashed")) +
    scale_x_log10(labels = si_lab()) +
    scale_y_log10(labels = si_lab()) +
    labs(title=paste0(title_prefix, " Q–Q (log–log)"),
         x="Fitted quantiles", y="Empirical quantiles",
         colour="", linetype="") +
    theme(legend.position="bottom")
  
  list(p_ccdf = p_ccdf, p_qq = p_qq)
}

# ---------- GB2 ----------
# Needs: pgb2(), gb2.mle(); (qgb2 closed-form provided here)
qgb2 <- function(u, a, b, p, q) {
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  w <- qbeta(u, shape1 = p, shape2 = q)
  b * (w/(1 - w))^(1/a)
}

gb2_plots <- function(x, start = NULL, params = NULL, save_dir = NULL, prefix = "gb2") {
  if (!exists("pgb2") || !exists("gb2.mle")) {
    stop("Please source pgb2() and gb2.mle() before calling.")
  }
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  
  if (is.null(params)) {
    if (is.null(start)) start <- c(a=1, b=median(x), p=2, q=2)
    fit <- gb2.mle(x,
                   start = start,
                   lower = c(a=1e-5, b=1e-8, p=1e-8, q=1e-8),
                   upper = c(a=Inf,  b=Inf,  p=Inf,  q=Inf))
    if (is.null(fit$par)) stop("GB2 fit failed (no parameters).")
    params <- setNames(as.numeric(fit$par), c("a","b","p","q"))
  }
  
  P <- function(t) pgb2(t, a=params["a"], b=params["b"], p=params["p"], q=params["q"])
  Q <- function(pp) qgb2(pp, a=params["a"], b=params["b"], p=params["p"], q=params["q"])
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "GB2")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  
  invisible(list(params = params, p_ccdf = plots$p_ccdf, p_qq = plots$p_qq))
}

# ---------- DPLN ----------
dpln_plots <- function(x, params = NULL, start = NULL, save_dir = NULL, prefix = "dpln") {
  # detect available CDF/Q names
  cdf_name <- if (exists("pdpln")) "pdpln" else if (exists("pdoubleparetolognormal")) "pdoubleparetolognormal" else NA
  q_name   <- if (exists("qdpln")) "qdpln" else if (exists("qdoubleparetolognormal")) "qdoubleparetolognormal" else NA
  if (is.na(cdf_name)) stop("Please source a DPLN CDF: pdpln() or pdoubleparetolognormal().")
  
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  
  # extract params if not given (fall back to mle)
  if (is.null(params)) {
    if (!exists("dpln.mle")) stop("params not provided and dpln.mle() not available.")
    if (is.null(start)) start <- c(meanlog = log(median(x)), sdlog = 1, shape1 = 1, shape2 = 1)
    fit <- dpln.mle(x, start = start)
    if (is.null(fit$par)) stop("DPLN fit failed (no parameters).")
    params <- setNames(as.numeric(fit$par), names(fit$par))
  }
  
  # ensure the names we expect (your example: meanlog, sdlog, shape1, shape2)
  need <- c("meanlog","sdlog","shape1","shape2")
  if (!all(need %in% names(params))) {
    stop("DPLN params must be named: meanlog, sdlog, shape1, shape2.")
  }
  
  # closures
  P <- function(t) do.call(get(cdf_name),
                           c(list(t),
                             as.list(params)[c("meanlog","sdlog","shape1","shape2")]))
  if (!is.na(q_name)) {
    Q <- function(pp) do.call(get(q_name),
                              c(list(pp),
                                as.list(params)[c("meanlog","sdlog","shape1","shape2")]))
  } else {
    xmin <- min(x)
    xmax <- quantile(x, 0.9999) * 100
    Q <- function(pp) q_from_cdf(pp, cdf_fun = P, lower = xmin, upper = xmax)
  }
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "DPLN")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  
  invisible(list(params = params, p_ccdf = plots$p_ccdf, p_qq = plots$p_qq))
}

# ---------- BR3 (a, b, p) — robust grid inversion for Q ----------
# Requires: pbr3() CDF
br3_plots <- function(x, params = NULL, start = NULL, save_dir = NULL, prefix = "br3") {
  if (!exists("pbr3")) stop("Please source pbr3() before calling.")
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  
  if (is.null(params)) {
    if (!exists("br3.mle")) stop("params not provided and br3.mle() not available.")
    if (is.null(start)) start <- c(a=1, b=median(x), p=1)
    fit <- br3.mle(x, start = start)
    if (is.null(fit$par)) stop("BR3 fit failed (no parameters).")
    params <- setNames(as.numeric(fit$par), names(fit$par))
  }
  need <- c("a","b","p")
  if (!all(need %in% names(params))) stop("BR3 params must be named: a, b, p.")
  
  # CDF closure
  P <- function(t) pbr3(t, a = params["a"], b = params["b"], p = params["p"])
  
  # --- robust quantile via grid inversion (no uniroot) ---
  Q <- function(pp) {
    pp <- pmin(pmax(pp, .Machine$double.eps), 1 - .Machine$double.eps)
    
    xmin <- min(x)
    xmax <- quantile(x, 0.9999) * 1e3  # generous upper range for tails
    # log grid over [xmin/2, xmax]
    lo <- max(xmin / 2, .Machine$double.eps)
    xg <- exp(seq(log(lo), log(xmax), length.out = 5000L))
    Fg <- P(xg)
    
    # enforce monotonicity & bounds
    Fg[!is.finite(Fg)] <- 0
    Fg <- pmin(pmax(Fg, 0), 1)
    Fg <- cummax(Fg)  # fix tiny numerical non-monotonicity
    
    # ensure pp inside [min(Fg), max(Fg)]
    loF <- max(min(Fg), .Machine$double.eps)
    hiF <- min(max(Fg), 1 - .Machine$double.eps)
    pp2 <- pmin(pmax(pp, loF), hiF)
    
    approx(x = Fg, y = xg, xout = pp2, ties = "ordered")$y
  }
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "BR3")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# ---------- BR12 (a, b, q) — robust grid inversion for Q ----------
# Requires: pbr12() CDF
br12_plots <- function(x, params = NULL, start = NULL, save_dir = NULL, prefix = "br12") {
  if (!exists("pbr12")) stop("Please source pbr12() before calling.")
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  
  if (is.null(params)) {
    if (!exists("br12.mle")) stop("params not provided and br12.mle() not available.")
    if (is.null(start)) start <- c(a=1, b=median(x), q=1)
    fit <- br12.mle(x, start = start)
    if (is.null(fit$par)) stop("BR12 fit failed (no parameters).")
    params <- setNames(as.numeric(fit$par), names(fit$par))
  }
  need <- c("a","b","q")
  if (!all(need %in% names(params))) stop("BR12 params must be named: a, b, q.")
  
  # CDF closure
  P <- function(t) pbr12(t, a = params["a"], b = params["b"], q = params["q"])
  
  # --- robust quantile via grid inversion (no uniroot) ---
  Q <- function(pp) {
    pp <- pmin(pmax(pp, .Machine$double.eps), 1 - .Machine$double.eps)
    
    xmin <- min(x)
    xmax <- quantile(x, 0.9999) * 1e3
    lo <- max(xmin / 2, .Machine$double.eps)
    xg <- exp(seq(log(lo), log(xmax), length.out = 5000L))
    Fg <- P(xg)
    
    Fg[!is.finite(Fg)] <- 0
    Fg <- pmin(pmax(Fg, 0), 1)
    Fg <- cummax(Fg)
    
    loF <- max(min(Fg), .Machine$double.eps)
    hiF <- min(max(Fg), 1 - .Machine$double.eps)
    pp2 <- pmin(pmax(pp, loF), hiF)
    
    approx(x = Fg, y = xg, xout = pp2, ties = "ordered")$y
  }
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "BR12")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# ---------- Fisk (beta, theta) ----------
# Requires: pfisk() CDF
fisk_plots <- function(x, params = NULL, start = NULL, save_dir = NULL, prefix = "fisk") {
  if (!exists("pfisk")) stop("Please source pfisk() before calling.")
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  
  if (is.null(params)) {
    if (!exists("fisk.mle")) stop("params not provided and fisk.mle() not available.")
    if (is.null(start)) start <- c(beta = median(x), theta = 1)
    fit <- fisk.mle(x, start = start)
    if (is.null(fit$par)) stop("Fisk fit failed (no parameters).")
    params <- setNames(as.numeric(fit$par), names(fit$par))
  }
  need <- c("beta","theta")
  if (!all(need %in% names(params))) stop("Fisk params must be named: beta, theta.")
  
  P <- function(t) pfisk(t, beta = params["beta"], theta = params["theta"])
  xmin <- min(x); xmax <- quantile(x, 0.9999) * 100
  Q <- function(pp) q_from_cdf(pp, cdf_fun = P, lower = xmin, upper = xmax)
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "Fisk")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# ---- Pareto (Type I) ----
# Expects: ppareto(), qpareto() exist with signature (..., xmin=, k=)
# params must be named c(xmin=..., k=...)
pareto_plots <- function(x, params, save_dir = NULL, prefix = "pareto") {
  if (!exists("ppareto") || !exists("qpareto")) stop("ppareto()/qpareto() must exist.")
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  if (!all(c("xmin","k") %in% names(params))) stop("Pareto params must be: xmin, k.")
  
  P <- function(t) ppareto(t, xmin = params["xmin"], k = params["k"])
  Q <- function(pp) qpareto(pp, xmin = params["xmin"], k = params["k"])
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "Pareto")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# ---- Lognormal ----
# Uses base R plnorm/qlnorm; params must be c(meanlog=..., sdlog=...)
lognormal_plots <- function(x, params, save_dir = NULL, prefix = "lognormal") {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  if (!all(c("meanlog","sdlog") %in% names(params))) stop("Lognormal params: meanlog, sdlog.")
  
  P <- function(t) plnorm(t, meanlog = params["meanlog"], sdlog = params["sdlog"])
  Q <- function(pp) qlnorm(pp, meanlog = params["meanlog"], sdlog = params["sdlog"])
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "Lognormal")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# ---- Weibull ----
# Uses base R pweibull/qweibull; params must be c(shape=..., scale=...)
weibull_plots <- function(x, params, save_dir = NULL, prefix = "weibull") {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  if (!all(c("shape","scale") %in% names(params))) stop("Weibull params: shape, scale.")
  
  P <- function(t) pweibull(t, shape = params["shape"], scale = params["scale"])
  Q <- function(pp) qweibull(pp, shape = params["shape"], scale = params["scale"])
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "Weibull")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# ---- Gamma ----
# Uses base R pgamma/qgamma; params must be c(shape=..., rate=...)
gamma_plots <- function(x, params, save_dir = NULL, prefix = "gamma") {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) stop("Not enough positive observations in 'x'.")
  if (!all(c("shape","rate") %in% names(params))) stop("Gamma params: shape, rate.")
  
  P <- function(t) pgamma(t, shape = params["shape"], rate = params["rate"])
  Q <- function(pp) qgamma(pp, shape = params["shape"], rate = params["rate"])
  
  plots <- .make_plots_ccdf_qq(x, cdf_fun = P, q_fun = Q, title_prefix = "Gamma")
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0(prefix, "_ccdf.png")), plots$p_ccdf, width=7, height=5, dpi=300)
    ggsave(file.path(save_dir, paste0(prefix, "_qq.png")),   plots$p_qq,   width=7, height=5, dpi=300)
  }
  invisible(plots)
}

# 1) GB2:
gb2_res <- gb2_plots(domsales2012, save_dir = output_dir, prefix = "gb2_fr2012")

# 2) DPLN with your coefficient object:
params_dpln <- unlist(doubleparetolognormal1$coefficients[[1]][[1]])
dpln_res <- dpln_plots(domsales2012, params = params_dpln, save_dir = output_dir, prefix = "dpln_fr2012")

# 3) BR3 (a,b,p)
params_br3  <- unlist(br31$coefficients[[1]][[1]])   # a, b, p
br3_res     <- br3_plots(domsales2012, params = params_br3, save_dir = output_dir, prefix = "br3_fr2012")
print(br3_res$p_ccdf); print(br3_res$p_qq)

# 4) BR12 (a,b,q)
params_br12 <- unlist(br121$coefficients[[1]][[1]])  # a, b, q
br12_res    <- br12_plots(domsales2012, params = params_br12, save_dir = output_dir, prefix = "br12_fr2012")
print(br12_res$p_ccdf); print(br12_res$p_qq)

# 5) Fisk (beta, theta)
params_fisk <- unlist(fisk1$coefficients[[1]][[1]])  # beta, theta
fisk_res    <- fisk_plots(domsales2012, params = params_fisk, save_dir = output_dir, prefix = "fisk_fr2012")

# 6) Pareto
coef_pareto   <- pareto1$coefficients[[1]][[1]]             # names: xmin, k
params_pareto <- c(xmin = as.numeric(coef_pareto["xmin"]),
                   k    = as.numeric(coef_pareto["k"]))
pareto_res <- pareto_plots(domsales2012, params = params_pareto,
                           save_dir = output_dir, prefix = "pareto_fr2012")

# 7) Lognormal
coef_lnorm   <- lnorm1$coefficients[[1]][[1]]               # 2x1, rows: meanlog, sdlog
params_lnorm <- c(meanlog = as.numeric(coef_lnorm["meanlog", 1]),
                  sdlog   = as.numeric(coef_lnorm["sdlog",   1]))
lognorm_res <- lognormal_plots(domsales2012, params = params_lnorm,
                               save_dir = output_dir, prefix = "lognormal_fr2012")
# 8) Gamma
coef_gamma   <- gamma1$coefficients[[1]][[1]]               # 2x1, rows: shape, rate
params_gamma <- c(shape = as.numeric(coef_gamma["shape", 1]),
                  rate  = as.numeric(coef_gamma["rate",  1]))
gamma_res <- gamma_plots(domsales2012, params = params_gamma,
                         save_dir = output_dir, prefix = "gamma_fr2012")

# 9) Weibull
coef_weib   <- weibull1$coefficients[[1]][[1]]              # 2x1, rows: shape, scale
params_weib <- c(shape = as.numeric(coef_weib["shape", 1]),
                 scale = as.numeric(coef_weib["scale", 1]))
weibull_res <- weibull_plots(domsales2012, params = params_weib,
                             save_dir = output_dir, prefix = "weibull_fr2012")
