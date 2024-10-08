% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadBoundaries_utils.R
\name{create_sector_df}
\alias{create_sector_df}
\title{Create sector data frame from circle data}
\usage{
create_sector_df(df)
}
\arguments{
\item{df}{A data frame containing circle data with segment IDs, x and y
locations.}
}
\value{
A data frame with the coordinates for each sector.
}
\description{
This function takes a data frame of circle data with segment IDs and
locations, and returns a data frame of sectors with their scaled coordinates.
}
