% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_site_years.R
\name{filter_site_years}
\alias{filter_site_years}
\title{filter_site_years}
\usage{
filter_site_years(x, site, year, min_times = 3)
}
\arguments{
\item{x}{data.frame}

\item{site}{column containing trial location}

\item{year}{column containing year}

\item{min_times}{numeric, minimum number of years a site must be present. Default 3.}
}
\value{
a data.frame, a subset of x.
}
\description{
Subset locations by a minimum number of occurrances
}
