% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measure_home_distance.R
\name{measure_home_distance}
\alias{measure_home_distance}
\title{measure_home_distance}
\usage{
measure_home_distance(
  data,
  locations,
  geno,
  site,
  lat = NA,
  long = NA,
  great_circle = TRUE
)
}
\arguments{
\item{data}{data.frame from a call to `id_home()`}

\item{locations}{a SpatialPointsDataFrame (recommended) or a data.frame with columns `site`, `lat`, and `long`.}

\item{geno}{genotype column in `data` and in `locations`}

\item{site}{spatial environment column in `data`}

\item{lat}{latitude/y coordinate column in `locations`}

\item{long}{longitude/x coordinate column in `locations`}

\item{great_circle}{are coordinates on an ellipse?}
}
\value{
a distance matrix. row/column names are genotypes. Values are distances
among between home sites of the genotype pairs.
}
\description{
calculate the distance between home sites for among genotypes
}
\note{
If lat/long are in degrees, then set `great_circle` to TRUE. If lat/long
are in linear units, set false. Can also pass a SpatialPointsDataFrame, which will
override this the `great_circle` setting and coordinates will be derived from the object.
}
\seealso{
package `"sp"`
}
