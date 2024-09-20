#' @title filter_site_years
#'
#' @description Subset locations by a minimum number of occurrances
#'
#' @param x data.frame
#' @param site column containing trial location
#' @param year column containing year
#' @param min_times numeric, minimum number of years a site must be present. Default 3.
#'
#' @return a data.frame, a subset of x.
#'
#' @export

filter_site_years <- function(x, site, year, min_times = 3) {

  # Combine site and year, then get unique occurrences
  common_sites <- unique(x[, c(site, year)])

  # Count the number of occurrences for each site
  site_counts <- table(common_sites[[site]])

  # Filter sites with at least min_times occurrences
  frequent_sites <- names(site_counts[site_counts >= min_times])

  # Subset the original data frame based on filtered sites
  filtered_data <- x[x[[site]] %in% frequent_sites, ]

  return(filtered_data)
}
