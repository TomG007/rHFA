#' @title Scale a Phenotype Column
#'
#' @description This function scales the values of a specified phenotype column in a data frame
#' by standardizing the data (z-scores), and then adds a new column with the scaled values.
#' The new column will have the name `rel_<pheno>`, where `<pheno>` is the name of the original phenotype.
#'
#' @param x A data frame containing the data to scale.
#' @param pheno A string specifying the name of the column in `x` to scale (the phenotype).
#'
#' @return A data frame with an additional column containing the scaled values of the specified phenotype.
.scale_pheno <- function(x, pheno) {
  # Create a new column name by prefixing 'rel_' to the name of the phenotype column
  colname <- paste0('rel_', pheno)

  # Scale the values of the specified phenotype column and assign them to the new column
  x[, (colname) := scale(get(pheno))]

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

  # Remove rows with missing observations in the phenotype column
  x <- x[!is.na(get(pheno)), ]

  # Ensure site column is a factor
  x[, (site) := droplevels(factor(get(site)))]

  # Check if BLUP is appropriate
  if (blup) {
    site_years <- x[, .N, by = get(site)]$N
    max_site_years <- max(site_years)

    if (max_site_years < 2) {
      if (verbose) message("Cannot use BLUP due to insufficient observations. Using mean values.")
      blup <- FALSE
    }
  }

  # Calculate the effect size (site-level mean) using BLUP or regular means
  if (blup) {
    if (verbose) message("Using linear mixed-effects model (BLUP) to identify the home site.")

    # Fit the linear mixed-effects model
    model_formula <- as.formula(paste0(pheno, " ~ (1 | ", site, ")"))
    model <- lme4::lmer(model_formula, data = x, control = lme4::lmerControl(calc.derivs = FALSE))

    # Extract random effects for the site
    site_effects <- coef(model)[[site]][, 1]
  } else {
    # Calculate the mean of the phenotype within each site
    site_effects <- tapply(get(pheno), get(site), mean, simplify = TRUE)
  }

  # Handle the case where site_effects are NA
  if (all(is.na(site_effects))) {
    if (verbose) message("All site effects are NA. Unable to identify a home site.")
    x[, is_home := FALSE]
  } else {
    # Identify the site with the maximum mean
    max_site <- names(which.max(site_effects))
    x[, is_home := get(site) == max_site]
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
id_home <- function(df, site, year, geno, pheno, blup = TRUE, verbose = FALSE) {

  print("this is the new one...")


  # Convert the data.frame to data.table if it's not already
  setDT(df)

  # Create the relative phenotype column name
  rel_colname <- paste0('rel_', pheno)

  # Center and scale performance within site-year using lapply
  site_year <- paste(df[[site]], df[[year]], sep = '_')
  df_list <- split(df, site_year)
  df_list <- lapply(df_list, .scale_pheno, pheno)
  df <- rbindlist(df_list)

  # Find the highest relative phenotype for each genotype
  geno_list <- split(df, df[[geno]])
  geno_list <- lapply(geno_list, .id_top_pheno, site, rel_colname, blup, verbose)
  df <- rbindlist(geno_list)

  return(df)
}

