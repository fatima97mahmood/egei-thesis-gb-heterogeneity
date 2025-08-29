# gb-heterogeneity-replication

Replication code for master’s thesis on firm-size heterogeneity, comparing Generalized Beta (GB/GB2) and Double Pareto Lognormal (DPLN) distributions using French firm sales (Orbis 2012).

---

# Replication README — Thesis Support Code

**Title:** Parameterising firm-level heterogeneity in quantitative trade models: Narrowing it down from the Generalized Beta distribution  
**Program:** EGEI — Erasmus Mundus Master’s Degree in Economics of Globalization and European Integration  
**Author:** Fatima Mahmood  
**Supervisor:** Ruben Dewitte  
**Submission:** September 2025

## Overview
This repository contains R code to reproduce the empirical analysis comparing the Generalized Beta (GB) family (GB, GB1, GB2, Beta4, BR12/BR3, etc.) with benchmarks such as Lognormal, Pareto, Weibull, Gamma, and the Double-Pareto Lognormal (DPLN). The main application uses French firm domestic sales (2012), normalized by the sample mean. The pipeline fits distributions, computes log-likelihood, NLL, AIC, BIC, KS, and produces CCDF and Q–Q diagnostics.

## Data availability & provenance
- Primary data: Orbis Europe (update 05/06/2020), unconsolidated accounts, France 2012, domestic sales.  
- Confidentiality: Not redistributed. Users must obtain access from Bureau van Dijk/Orbis via their institution.  
- Filters: Strictly positive domestic sales; missing values removed; normalized by sample mean.  
- Expected file: `empirics/input/data/orbis_france.Rdata` (contains a vector of domestic sales; see *Setup gotcha* below).

If you don’t have Orbis, you can still run the code to replicate model fitting on your own numeric vector, and generate figures/tables from saved fit objects.


## Software requirements
- R ≥ 4.0  
- Packages: `tidyverse`, `here`, `xtable`, `latex2exp`, `ggpubr`, `distributionsrd`, `dplyr`, `purrr`, `ggplot2`, `scales`, `MASS`, `flexmix`  
- Suggested reproducibility: add `renv` lockfile.

### Quick setup
```r
install.packages(c(
  "tidyverse","here","xtable","latex2exp","ggpubr","distributionsrd",
  "dplyr","purrr","ggplot2","scales","MASS","flexmix"
))

```
## Instructions to Replicators

### Project root and paths
Open **ReplicationFiles.Rproj** (so `here::here()` resolves to the repo root).  
Remove or comment the absolute `setwd(...)` lines in scripts — `here()` already points correctly.

---

### Place data
Put `orbis_france.Rdata` under: empirics/input/data/orbis_france.Rdata


---

### Run the main script
```r
source(here::here("empirics","code","domsales_france.R"))
```

## This will

- Load data, keep positive values, remove NA, normalize by mean.  
- Fit distributions via `combdist.mle_update`, `combdist.mle_1param/2param/3param/4_5_params`.  
- Save per-distribution `.Rdata` fit objects in `empirics/output/domsales/`.  
- Aggregate, rank (NLL/AIC/BIC), perform KS tests, and save `empirics/output/eval_cities/loglike_ranked.Rdata`.  
- Write a LaTeX summary table (via `xtable`).  
- Produce CCDF and Q–Q figures under `empirics/output/distribution_fits/`.  

---

## Key script(s)

- **`empirics/code/domsales_france.R`**  
  Main driver, loads helpers, data, fits, saves results, plots.  

- **`empirics/code/axillary_code.R`**  
  Implements additional distributions (`unif1`, `rayleigh`, `loglogistic`, `power`, `fisk`, `lomax`, `invlomax`, `IB1`, `UG`, `B1`, `GG`, `B2`, `BR12`, `BR3`, `GB1`, `GB2`, `Beta4`, `GB`) with PDF/CDF/MLE, wrappers, KS tests, plotting utilities.  

- **Other modules (`toolbox*.R`, `melitz.R`, `nmad.R`, `FLXMC_xtra.R`)**  
  Supporting functions for mixtures, Frechet, Melitz-style utilities.  

---

## Expected outputs

- **Fit objects:** `empirics/output/domsales/<dist><components>.Rdata`  
- **Ranked results:** `empirics/output/eval_cities/loglike_ranked.Rdata`  
- **Figures:** `empirics/output/distribution_fits/*_{ccdf,qq}.png`  

---

## Mapping of programs → outputs

| Type    | Script              | Output                                         | Notes |
|---------|---------------------|------------------------------------------------|-------|
| Fits    | `domsales_france.R` | `empirics/output/domsales/*.Rdata`             | One file per distribution |
| Summary | `domsales_france.R` | `empirics/output/eval_cities/loglike_ranked.Rdata` | Includes NLL/AIC/BIC/KS |
| LaTeX   | `domsales_france.R` | console output (copy/paste to `.tex`)          | Goodness-of-fit table |
| Plots   | `axillary_code.R`   | `empirics/output/distribution_fits/*.png`      | CCDF + Q–Q plots |

---

## Reproducibility notes

- Add `set.seed(123)` at the top of `domsales_france.R` for stable results.  
- Optimisers use `optim(..., method="L-BFGS-B")` with sensible bounds.  
- KS tests rely on analytic CDFs where available, otherwise numeric integration.  
- GB-family and BR variants may be sensitive to starting values; defaults provided.  

---

## Setup gotchas

In `domsales_france.R` the script loads to an environment but later references `x$domsales`.  
Unless `orbis_france.Rdata` creates an object `x`, you should replace with:

```r
load(here("empirics","input","data","orbis_france.Rdata"))
domsales2012 <- domsales   # if the object is named domsales
# or
domsales2012 <- data$domsales
```

## How to run on your own data

Replace the data block with any positive numeric vector:

```r
y <- YOUR_VECTOR
y <- y[y > 0]
y <- na.omit(y)
y <- y / mean(y)
xy <- data.frame(var="y", ind="Total", varind="y.Total", year=2012, value=y)
```

Then run the same fitting and plotting functions.

## Citation
@thesis{Mahmood2025GB,
  title={Parameterising firm-level heterogeneity in quantitative trade models: Narrowing it down from the Generalized Beta distribution},
  author={Mahmood, Fatima},
  school={EGEI},
  year={2025},
  month={September}
}

