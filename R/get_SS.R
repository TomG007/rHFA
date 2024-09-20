#' @title get_ss
#'
#' @description calculate the sum of squares for a model and present a summary of
#' the predictors, their sum of squares, percentage variance explained, F-values,
#' and p-values in a data frame. It relies on the Anova function from the car package.
#'
#' @param model A data frame containing the model data.
#'
#' @return A data frame with the statistical measures
#'
#' @import car
#' @export

get_ss <- function(model) {

  # Check if model is valid
  if (is.null(model)) {
    stop("The input model is NULL. Please provide a valid model.")
  }

  # Perform Anova
  a <- Anova(model)

  # Check if Anova has results
  if (is.null(a) || nrow(a) == 0) {
    stop("No valid ANOVA results. Check your model.")
  }

  # Create output data frame
  out <- data.frame(
    Predictor = rownames(a),
    SumSq = a$`Sum Sq`,
    PercentVariance = round(a$`Sum Sq` / sum(a$`Sum Sq`) * 100, 2),
    F_value = round(a$`F value`, 4),
    p_value = signif(a$`Pr(>F)`, 3)
  )

# ----- Original dataframe column names --------
#  # Create output data frame
#  out <- data.frame(
#    PREDICTOR = rownames(a),
#    SUMSQ = a$`Sum Sq`,
#    pVar = round(a$`Sum Sq` / sum(a$`Sum Sq`) * 100, 2),
#    F_val = round(a$`F value`, 4),
#    p_val = signif(a$`Pr(>F)`, 3)
#  )

  return(out)
}
