#' @title get_home_site
#'
#' @description extract home location information for each genotype
#'
#' @param data data.frame from call to `id_home()`
#' @param geno column in `data` with genotype
#' @param site column in `data` with location environment
#'
#' @return a `data.frame` with columns `geno`, `site`, and `home_years` - the number
#' of years a genotype appeared at a home site.
#'
#' @export

get_home_site <- function(data, geno, site) {

  # Check if 'is_home' column exists
  if (!'is_home' %in% names(data)) {
    stop("'is_home' column not found in the data")
  }

  # Filter rows where 'is_home' is TRUE
  out <- data[data$is_home == TRUE, c(geno, site)]

  # Count the number of occurrences of each genotype
  geno_counts <- table(out[[geno]])

  # Get unique genotype-site pairs
  out <- unique(out)

  # Add 'home_years' as the count of years each genotype appeared
  out$home_years <- geno_counts[out[[geno]]]

  return(out)
}
