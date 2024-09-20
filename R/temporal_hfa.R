#' @title Calculate Home Field Advantage (HFA) for Each Year
#'
#' @description This function calculates the home field advantage (HFA) coefficients for each year based on the specified formula.
#'
#' @param x A data frame containing the data.
#' @param ff A formula to use in the linear model, which should include the term <year>_num:is_home.
#' @param year A character string representing the year column in the formula.
#'
#' @details The formula should contain a term such as <year>_num:is_home and specify a response variable (e.g., phenotype).
#'
#' @return A data frame containing the HFA coefficients, standard errors, p-values, and adjusted p-values for each year.
#'

.calculate_temporal_hfa <- function(x, ff, year) {

  # Extract the response variable from the formula
  response_var <- as.character(ff)[2]

  # Fit the linear model using the specified formula
  model <- lm(ff, data = x)

  # Extract the summary and coefficients
  model_summary <- summary(model)
  home_coef <- model_summary$coefficients

  # Filter rows with the ':is_homeTRUE' pattern
  home_coef <- home_coef[grepl(':is_homeTRUE', rownames(home_coef)), ]

  # Clean up the row names to remove ':is_homeTRUE' and the year from the coefficients
  rownames(home_coef) <- gsub(':is_homeTRUE', '', rownames(home_coef))
  rownames(home_coef) <- gsub(year, '', rownames(home_coef))

  # Convert the row names (years) to a column and adjust p-values
  home_coef_df <- data.frame(
    year = rownames(home_coef),
    year_num = type.convert(rownames(home_coef)),
    home_coef,
    p.adj = p.adjust(home_coef[, 'Pr(>|t|)']),
    stringsAsFactors = FALSE
  )

  # Clean up column names for clarity
  colnames(home_coef_df) <- gsub('year', year, colnames(home_coef_df))
  colnames(home_coef_df) <- gsub('Pr...t..', 'p.value', colnames(home_coef_df))
  colnames(home_coef_df) <- gsub('Std..Error', 'Std.Error', colnames(home_coef_df))

  # Return the formatted data frame with coefficients and relevant statistics
  return(home_coef_df)
}

