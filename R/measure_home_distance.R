#' @title measure_home_distance
#'
#' @description calculate the distance between home sites for among genotypes
#'
#' @param data data.frame from a call to `id_home()`
#' @param locations a SpatialPointsDataFrame (recommended) or a data.frame with columns `site`, `lat`, and `long`.
#' @param geno genotype column in `data` and in `locations`
#' @param site spatial environment column in `data`
#' @param lat latitude/y coordinate column in `locations`
#' @param long longitude/x coordinate column in `locations`
#' @param great_circle are coordinates on an ellipse?
#'
#' @return a distance matrix. row/column names are genotypes. Values are distances
#' among between home sites of the genotype pairs.
#'
#' @note If lat/long are in degrees, then set `great_circle` to TRUE. If lat/long
#' are in linear units, set false. Can also pass a SpatialPointsDataFrame, which will
#' override this the `great_circle` setting and coordinates will be derived from the object.
#'
#' @seealso package `"sp"`
#'
#'
#' @import sp
#' @export

measure_home_distance <- function(data, locations, geno, site, lat = NA, long = NA, great_circle = TRUE) {
  # Load required packages

  # Check if locations is a SpatialPointsDataFrame
  if (inherits(locations, "SpatialPoints")) {
    distances <- spDists(locations)
  } else {
    # Ensure lat and long columns are provided
    if (is.na(lat) || is.na(long)) {
      stop("Please provide 'lat' and 'long' column names if not using SpatialPointsDataFrame")
    }

    # Calculate distances based on lat/long
    distances <- as.matrix(locations[, c(long, lat)])
    distances <- spDists(distances, longlat = great_circle)
  }

  # Set row and column names for distances
  rownames(distances) <- locations[[site]]
  colnames(distances) <- locations[[site]]

  # Subset data to get home sites
  home_data <- get_home_site(data, geno, site)

  # Extract genotype and site names
  genotypes <- as.character(home_data[[geno]])
  home_sites <- setNames(as.character(home_data[[site]]), genotypes)

  # Initialize output matrix with row and column names
  n_genotypes <- length(genotypes)
  distance_matrix <- matrix(NA, nrow = n_genotypes, ncol = n_genotypes,
                            dimnames = list(genotypes, genotypes))

  # Fill distance matrix
  for (i in genotypes) {
    home_i <- home_sites[i]
    for (j in genotypes) {
      home_j <- home_sites[j]
      distance_matrix[i, j] <- distances[home_i, home_j]
    }
  }

  return(distance_matrix)
}
