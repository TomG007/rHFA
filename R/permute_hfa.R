#' @title Calculate Residuals of a Linear Model
#'
#' @description This function calculates the residuals of a linear model using base R functions.
#' These residuals can be used in subsequent regressions as partial effects.
#'
#' @param ff A formula describing the linear model.
#' @param data A data frame containing the variables in the formula.
#'
#' @return A vector of residuals with length equal to the number of rows in the data.
#'

.quick_resid <- function(ff, data) {

  # Create the model matrix from the formula and data
  model_matrix <- model.matrix(ff, data)

  # Perform QR decomposition on the model matrix
  qr_decomp <- qr(model_matrix)

  # Extract the response variable (Y) from the data frame using the formula
  response_vector <- model.frame(ff, data)[, 1]

  # Calculate the residuals using QR decomposition
  residuals <- qr.resid(qr_decomp, response_vector)

  return(residuals)
}


#' @title Generate Permutation Sets
#'
#' @description Generates structured permutation sets based on site-year combinations.
#'
#' @param x A data frame containing the data to permute.
#' @param site A character string identifying the column with spatial environment information.
#' @param year A character string identifying the column with temporal environment information.
#' @param times An integer indicating the number of permutation sets to generate.
#' @param seed An optional integer for setting the seed for reproducibility (default is NULL).
#'
#' @return A data frame of permutation sets, where row numbers are used to reorder `x`. The first column contains the original data order.
#'

.generate_sets <- function(x, site, year, times, seed = NULL) {

  # Create a control variable by combining site and year information
  control <- as.factor(paste(x[, site], x[, year], sep = '_'))

  # Get the number of rows in the data frame
  N <- nrow(x)

  # Set the seed for reproducibility, if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate permutation sets using the `permute` package
  permuted_sets <- shuffleSet(N, times, control = how(blocks = control))

  # Include the original data order as the first row
  original_data <- seq_len(N)
  all_sets <- rbind(original_data, permuted_sets)

  # Transpose and convert to data frame
  result <- as.data.frame(t(all_sets))

  return(result)
}


#' @title Calculate Home Field Advantage (HFA)
#'
#' @description This function calculates the home field advantage based on the provided formula.
#' It uses efficient matrix operations for performance.
#'
#' @param x A data frame containing the data to calculate HFA.
#' @param ff A formula to use (e.g., part_pheno ~ geno + is_home).
#' @param part_pheno A character string representing the phenotype column name.
#' @param geno A character string representing the genotype column name.
#' @param year A character string representing the year column name.
#' @param LAPACK A logical indicating whether to use LAPACK for QR decomposition (default is TRUE).
#'
#' @return A matrix of the home field advantage coefficients.
#'
.calculate_hfa <- function(x, ff, part_pheno, geno, year, LAPACK = TRUE) {

  # Create the model matrix using the formula and data
  model_matrix <- model.matrix(ff, data = droplevels(x))

  # Perform QR decomposition using LAPACK if specified
  qr_decomp <- qr(model_matrix, LAPACK = LAPACK)

  # Extract the response variable (phenotype) from the data
  response_vector <- x[, part_pheno, drop = FALSE]

  # Calculate the coefficients using the QR decomposition
  home_coef <- as.matrix(qr.coef(qr_decomp, response_vector))

  # Identify rows corresponding to 'is_homeTRUE' in the coefficients
  is_home <- grepl('is_homeTRUE', rownames(home_coef))
  home_coef <- home_coef[is_home, , drop = FALSE]

  # Clean up the row names to remove unnecessary parts
  rownames(home_coef) <- gsub(':is_homeTRUE', '', rownames(home_coef))
  rownames(home_coef) <- gsub(geno, '', rownames(home_coef))
  rownames(home_coef) <- gsub(year, '', rownames(home_coef))
  rownames(home_coef) <- gsub('TRUE', '', rownames(home_coef))

  return(home_coef)
}


#' @title Perform Two-Tailed Test on Permutation Values
#'
#' @description This function performs a two-tailed test on a matrix of permutation values.
#' Each column in the matrix represents the result of a different permutation, and the first column contains the original data.
#'
#' @param x A matrix where each column represents a permutation, and the first column is the original data.
#'
#' @return A vector of p-values for each row in the matrix.
#'
.two_tailed <- function(x) {

  # Calculate the proportion of permutations where the value is greater than or equal to the original data
  alpha <- rowSums(sweep(x, 1, x[, 1], `>=`)) / ncol(x)

  # Calculate two-tailed p-values
  p_val <- apply(cbind(alpha, 1 - alpha), 1, min) * 2

  return(p_val)
}


#' @title Calculate Confidence Intervals for Permutation Results
#'
#' @description This function calculates the median and 90% confidence intervals
#' (5th and 95th percentiles) for the difference between observed and permutation data.
#'
#' @param x A matrix of permutation results, where rows represent values and columns represent permutations.
#' The first column contains the observed data.
#'
#' @return A matrix with columns for the median, 5th percentile (q05), and 95th percentile (q95) of (observed - permutation).
#'

.calculate_intervals <- function(x) {

  # Calculate the difference between observed and permutations (observed - permutation)
  difference <- sweep(x, 1, x[, 1], '-') * -1

  # Compute the median, 5th percentile (q05), and 95th percentile (q95) across rows
  result <- cbind(
    median = apply(difference, 1, median, na.rm = TRUE),
    p05 = apply(difference, 1, quantile, 0.05, na.rm = TRUE),
    p95 = apply(difference, 1, quantile, 0.95, na.rm = TRUE)
  )

  return(result)
}


#' @title Test the Magnitude and Significance of Home Field Advantage
#'
#' @description Test the magnitude and significance of the home field advantage
#' versus what would be expected by chance using permutations.
#'
#' @param data A data frame containing the data.
#' @param level A parameter at which to apply the home field advantage. Default is "population".
#' @param site The column name containing spatial environmental information.
#' @param year The column name containing temporal environmental information.
#' @param geno The column name containing genotype or variety information.
#' @param pheno The column name containing phenotype or performance information.
#' @param popn The column name differentiating sub-populations of genotypes within the dataset.
#' @param times The number of permutations to generate (default is 99).
#' @param blup_home Logical; whether to use shrinkage estimates for home field advantage (default is TRUE).
#' @param parallel Logical; whether to perform computations in parallel (default is TRUE).
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with test results and permutation data.
#'
#' @import Matrix
#' @import permute
#' @import parallel
#' @importFrom stats as.formula lm median model.frame model.matrix na.omit p.adjust quantile setNames
#' @importFrom utils type.convert
#' @export

permute_hfa <- function(data,
                        level = c('population', 'genotype', 'year', 'site'),
                        site = NA,
                        year = NA,
                        geno = NA,
                        pheno = NA,
                        popn = NA,
                        times = 99,
                        blup_home = TRUE,
                        parallel = TRUE,
                        seed = NULL) {

  # Check for valid column names in the data
  if (!all(c(site, year, geno, pheno) %in% colnames(data))) {
    stop("One or more of the provided column names do not exist in the data frame.")
  }

  # Set up parallel processing
  ncpu <- ifelse(parallel, detectCores() - 1, 1)

  # Check platform compatibility for parallelism
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(ncpu)
  } else {
    cl <- NULL
  }

  # Match the level argument
  level <- match.arg(level)

  # Set LAPACK option based on the level
  LAPACK <- !(level %in% c('site', 'year'))

  # Create new column names for partial and relative phenotypes
  part_pheno <- paste0('part_', pheno)
  rel_pheno <- paste0('rel_', pheno)

  # Construct the formula for HFA effect
  ff <- switch(level,
               'population' = as.formula(paste("~", geno, "+ is_home")),
               'genotype' = as.formula(paste("~", geno, "+", geno, ":is_home")),
               'year' = as.formula(paste("~", geno, "+ year:is_home")),
               'site' = as.formula(paste("~", geno, "+ site:is_home")))

  # Select and format data into list of data frames for each population


  #
  #  begin new code
  #

  if (!is.na(popn)) {
    data_subset <- na.omit(data[, c(site, year, geno, pheno, popn)])
  } else {
    data_subset <- na.omit(data[, c(site, year, geno, pheno)])
  }

  print(str(data_subset))  # Check structure after subsetting

  setDT(data_subset)  # Convert to data.table

  #
  #  end new code, uncomment following line
  #


#  data_subset <- na.omit(data[, c(site, year, geno, pheno, popn)])

  if (!is.na(popn)) {
    data_list <- split(data_subset, data_subset[, popn])
  } else {
    data_list <- list(data_subset)
  }

  # Identify if home site needs to be determined
  id_home_site <- !any(grepl('is_home', colnames(data)))

  if (id_home_site) {
    if (.Platform$OS.type == "windows") {
      data_list <- parLapply(cl, data_list, function(x) {
        id_home(x, site, year, geno, pheno, blup = blup_home, verbose = FALSE)
      })
    } else {
      data_list <- mclapply(data_list, function(x) {
        id_home(x, site, year, geno, pheno, blup = blup_home, verbose = FALSE)
      }, mc.cores = ncpu)
    }
  }

  # Calculate partial and relative yields
  if (.Platform$OS.type == "windows") {
    data_list <- parLapply(cl, data_list, function(x) {
      x <- id_home(x, site, year, geno, pheno, blup_home, verbose = FALSE)
      x[, part_pheno] <- .quick_resid(as.formula(paste(pheno, '~', year, '*', site)), data = x)
      x[, part_pheno] <- scale(x[, part_pheno], scale = FALSE)
      return(x)
    })
  } else {
    data_list <- mclapply(data_list, function(x) {
      x <- id_home(x, site, year, geno, pheno, blup_home, verbose = FALSE)
      x[, part_pheno] <- .quick_resid(as.formula(paste(pheno, '~', year, '*', site)), data = x)
      x[, part_pheno] <- scale(x[, part_pheno], scale = FALSE)
      return(x)
    }, mc.cores = ncpu)
  }

  # Permute HFA within each population
  results <- lapply(data_list, function(x) {
    sets <- .generate_sets(x, site, year, times, seed)

    coef_permute <- if (.Platform$OS.type == "windows") {
      parLapply(cl, sets, function(ss) {
        x[, c(rel_pheno, part_pheno)] <- x[ss, ]
        x <- split(x, x[, geno])
        x <- lapply(x, .id_top_pheno, site, rel_pheno, blup = blup_home, verbose = FALSE)
        x <- do.call(rbind, x)
        .calculate_hfa(x, ff, part_pheno, geno, year, LAPACK)
      })
    } else {
      mclapply(sets, function(ss) {
        x[, c(rel_pheno, part_pheno)] <- x[ss, ]
        x <- split(x, x[, geno])
        x <- lapply(x, .id_top_pheno, site, rel_pheno, blup = blup_home, verbose = FALSE)
        x <- do.call(rbind, x)
        .calculate_hfa(x, ff, part_pheno, geno, year, LAPACK)
      }, mc.cores = ncpu)
    }

    coef_permute <- do.call(cbind, coef_permute)
    colnames(coef_permute) <- c('observed', paste0('perm', 1:(ncol(coef_permute) - 1)))

    test <- cbind(
      intervals = .calculate_intervals(coef_permute),
      p_val = .two_tailed(coef_permute)
    )

    list(results = test, perms = coef_permute)
  })

  if (.Platform$OS.type == "windows") {
    stopCluster(cl)
  }

  if (level == 'population') {
    test_results <- do.call(rbind, lapply(results, `[[`, "results"))
    perms <- do.call(rbind, lapply(results, `[[`, "perms"))
  } else {
    lvl <- switch(level, 'genotype' = geno, 'year' = year, 'site' = site)
    test_results <- do.call(rbind, lapply(results, function(res) {
      data.frame(popn = names(results), lvl = rownames(res$results), res$results)
    }))
    perms <- lapply(results, `[[`, "perms")
  }

  return(list(home_field = test_results, perms = perms))
}
