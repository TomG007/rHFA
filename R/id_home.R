#' Scale a Phenotype Column
#'
#' @description This function scales the values of a specified phenotype column in a data frame
#' by standardizing the data (z-scores), and then adds a new column with the scaled values.
#' The new column will have the name `rel_<pheno>`, where `<pheno>` is the name of the original phenotype.
#'
#' @param x A data frame containing the data to scale.
#' @param pheno A string specifying the name of the column in `x` to scale (the phenotype).
#'
#' @return A data frame with an additional column containing the scaled values of the specified phenotype.
#' The new column will be named `rel_<pheno>`.

.scale_pheno <- function(x, pheno) {
  # Create a new column name by prefixing 'rel_' to the name of the phenotype column
  colname <- paste0('rel_', pheno)

  # Scale the values of the specified phenotype column using get()
  x[, (colname) := scale(get(pheno))]

  # Return the updated data frame
  return(x)
}


#' @title Identify Top Phenotype Site (Home Site)
#'
#' @description Identifies the level of `x[, site]` where the mean relative value of
#' `x[, pheno]` is highest. Uses shrinkage estimates of mean relative `x[, pheno]`
#' via `lme4` if `blup = TRUE` and if any site occurs in at least two years.
#' Otherwise, calculates simple means.
#'
#' @param x A data.frame or matrix representing a single unit (e.g., a genotype). It is coerced to a data.frame.
#' @param site A character string specifying the column of the grouping variable (e.g., field site locations).
#' @param pheno A character string specifying the column of the response variable (e.g., yield, biomass).
#' @param blup Logical; if `TRUE`, uses shrinkage estimates from a linear mixed model.
#' @param verbose Logical; if `TRUE`, prints messages to the console about the process.
#'
#' @return The input data frame with an additional logical column `"is_home"`, indicating
#' whether each observation corresponds to the identified home site.

.id_top_pheno <- function(x, site, pheno, blup = TRUE, verbose = TRUE) {

  # Coerce to data.frame and ensure factors are properly handled
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[, site] <- factor(x[, site])
  x[, site] <- droplevels(x[, site])

  # Check if BLUP is appropriate based on the number of site occurrences
  if (blup) {
    site_years <- tapply(x[, site], x[, site], length)
    max_site_years <- max(site_years)

    # Use means if insufficient site years for BLUP
    if (max_site_years < 3) {
      if (verbose) message("Insufficient site data for BLUP. Using mean values.")
      blup <- FALSE
    }
  }

  # Calculate the effect size (site-level mean) based on BLUP or regular means
  if (blup) {
    if (verbose) message("Using linear mixed-effects model (BLUP) to identify the home site.")

    # Build formula and fit the model
    model_formula <- as.formula(paste(pheno, "~ (1 |", site, ")"))
    model <- lmer(model_formula, data = x, control = lmerControl(calc.derivs = FALSE))

    site_effects <- coef(model)[[site]][, 1]
    names(site_effects) <- rownames(coef(model)[[site]])

  } else {
    # Use simple means across the sites
    site_effects <- tapply(x[, pheno], x[, site], mean, simplify = TRUE)
  }

  # Handle case where all site_effects are NA
  if (all(is.na(site_effects))) {
    if (verbose) message("All site effects are NA. Unable to identify a home site.")
    x$is_home <- FALSE  # or another appropriate value
  } else {
    # Identify the site with the maximum mean
    max_site <- names(site_effects)[which.max(site_effects)]
    x$is_home <- x[, site] == max_site
  }

  return(x)
}



#' @title Identify the Home Site Based on the Highest Relative Phenotype Value
#'
#' @description This function identifies the home site (location) where a given
#' variety performs best relative to other varieties across years. It can use
#' shrinkage estimates (BLUP) if required.
#'
#' @param df A data frame containing performance data by site, year, variety, etc.
#' @param site A character string indicating the column name for spatial units (e.g., field sites).
#' @param year A character string indicating the column name for temporal units (e.g., year).
#' @param geno A character string indicating the column name for genetic units (e.g., genotype).
#' @param pheno A character string indicating the column name for phenotype (e.g., yield, biomass).
#' @param blup A logical indicating whether to use shrinkage estimates (default is TRUE).
#' @param verbose A logical indicating whether to print messages (default is TRUE).
#'
#' @return A data frame with two new columns:
#' - `'rel_<pheno>'`: The relative phenotype value of each variety within each site-year.
#' - `'is_home'`: Logical column indicating whether a site is a variety's home site.
#'
#' @import lme4
#' @importFrom data.table rbindlist setDT
#' @importFrom stats coef
#' @export

id_home <- function(df, site, year, geno, pheno, blup = TRUE, verbose = TRUE) {

  # Convert the data.frame to data.table if it's not already
  setDT(df)

  # Create the relative phenotype column name
  rel_colname <- paste0('rel_', pheno)

  # Make site-year vector for splitting the data
  df[, site_year := paste0(get(site), "_", get(year))]

  # Center and scale performance within site-year using lapply
  df_list <- split(df, by = "site_year")
  df_list <- lapply(df_list, function(x) .scale_pheno(x, pheno))

  # Merge the data frames back together
  df <- rbindlist(df_list)

  # Find the highest relative phenotype for each genotype
  geno_list <- split(df, by = geno)
  geno_list <- lapply(geno_list, function(x) .id_top_pheno(x, site, rel_colname, blup, verbose))

  # Merge the data frames back together
  df <- rbindlist(geno_list)

  return(df)
}
