#' Normalized Absolute Deviation
#'
#' Calculates the Normalized Absolute Deviation between the empirical moments and the moments of the provided distribution. Corresponds to the Kolmogorov-Smirnov test statistic for the zeroth moment.
#' @param x data vector
#' @param r moment parameter
#' @param dist character vector containing distribution
#' @param prior named list of priors, defaults to 1
#' @param coeff named list of coefficients
#' @param stat character vector indicating which statistic should be calculated: none (NULL), the maximum deviation "max" or the sum of deviations "sum". Defaults to NULL.
#' @param ... Additional arguments can be passed to the parametric moment call.
#'
#' @examples
#'
#' x <- rlnorm(1e2, meanlog = -0.5, sdlog = 0.5)
#' nmad_test(x = x, r = 0, dist = "lnorm", coeff = c(meanlog = -0.5, sdlog = 0.5))
#' nmad_test(x = x, r = 0, dist = "lnorm", coeff = c(meanlog = -0.5, sdlog = 0.5), stat = "max")
#' nmad_test(x = x, r = 0, dist = "lnorm", coeff = c(meanlog = -0.5, sdlog = 0.5), stat = "sum")
#' @export

nmad_test <- function(x, r = 0, dist, prior = 1, coeff, lowertrunc=NULL, uppertrunc=NULL, stat = c("NULL", "max", "sum")) {
  stat <- match.arg(stat)

  empirical.CDF <- mempirical(data = x, truncation = x, r = r, lower.tail = FALSE)
  distr.CDF <- mcombdist(r = r, truncation = x, dist = dist, prior = prior, coeff = coeff, lower.tail = FALSE, lowertrunc=lowertrunc, uppertrunc=uppertrunc)

  if (any(is.na(distr.CDF))) {
    return(NA)
  } else {
    nmad <- abs(empirical.CDF - distr.CDF) / max(empirical.CDF)

    if (stat == "max") {
      nmad <- max(nmad, na.rm = TRUE)
    } else if (stat == "sum") {
      nmad <- sum(nmad, na.rm = TRUE)
    }
  }

  return(nmad)
}
