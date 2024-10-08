% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permute_hfa.R
\name{.quick_resid}
\alias{.quick_resid}
\title{Calculate Residuals of a Linear Model}
\usage{
.quick_resid(ff, data)
}
\arguments{
\item{ff}{A formula describing the linear model.}

\item{data}{A data frame containing the variables in the formula.}
}
\value{
A vector of residuals with length equal to the number of rows in the data.
}
\description{
This function calculates the residuals of a linear model using base R functions.
These residuals can be used in subsequent regressions as partial effects.
}
