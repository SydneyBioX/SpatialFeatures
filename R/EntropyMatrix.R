#' Calculate the Entropy of a Counts Matrix
#'
#' @description 
#' This function computes the entropy of a counts matrix from a Molecule Experiment object based on the given assay type.
#'
#' @param me A Molecule Experiment object.
#' @param assayName A character string specifying the assay type. Supported values include 
#' "sub-sector", "sub-concentric", "sub-combo", "super-concentric", and "super-combo".
#'
#' @return A matrix representing the entropy values corresponding to the given assay type.
#' @export
#' @examples
#' \dontrun{
#' # Assuming `data_obj` is your Molecule Experiment object
#' em = EntropyMatrix(data_obj, assayName = "sub-sector")
#' }
EntropyMatrix <- function(me, assayName, ...) {
  # Ensure assayName is valid
  if (!assayName %in% c("sub_sector", "sub_concentric", 'sub_combo', "super_sector", "super_concentric", "super_combo")) {
    stop("Invalid assayName provided!")
  }
  # Assuming 'me' is globally available or should pass it as an argument
  countsmatrix = MoleculeExperiment::countMolecules(me, moleculesAssay = "detected", boundariesAssay = assayName, matrixOnly = TRUE, ...)
  
  # Extract the counts matrix
  counts_matrix <- CountsMatrix(assayName, countsmatrix)
  entropy_matrix <- matrix_entropy(counts_matrix)
  return(entropy_matrix)
}