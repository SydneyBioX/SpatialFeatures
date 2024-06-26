% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities_MCC.R
\name{get_neighbors}
\alias{get_neighbors}
\title{Get Neighbors for a Given Cell}
\usage{
get_neighbors(cell_id, df, threshold = 2)
}
\arguments{
\item{cell_id}{The identifier of the cell to find neighbors for.}

\item{df}{A data frame with cell information including 'segment_id' and 'sample_id'.}

\item{threshold}{Distance threshold to consider cells as neighbors (default is 2).}
}
\value{
A vector of segment IDs representing the neighbors of the specified cell.
}
\description{
Get Neighbors for a Given Cell
}
