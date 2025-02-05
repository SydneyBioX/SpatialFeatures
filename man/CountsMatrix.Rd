% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EntropyMatrix_utils.R
\name{CountsMatrix}
\alias{CountsMatrix}
\title{Extract Counts Matrix from MoleculeExperiment based on Assay Type}
\usage{
CountsMatrix(me, assayName, nCores = 1, ...)
}
\arguments{
\item{me}{A MoleculeExperiment object.}

\item{assayName}{A character string indicating the assay type. Supported
values include "sub-sector", "sub-concentric", "sub-combo",
"super-concentric", and "super-combo".}

\item{nCores}{Number of cores for parallel processing (default 1)}

\item{...}{arguments passing to \code{MoleculeExperiment::countMolecules}}
}
\value{
A counts matrix corresponding to the specified assay type.
}
\description{
This function retrieves a counts matrix from a MoleculeExperiment object
based on the given assay type.
}
\examples{
# Assuming `data_obj` is your Molecule Experiment object
# cm = CountsMatrix(data_obj, assayName = "sub-sector")
}
