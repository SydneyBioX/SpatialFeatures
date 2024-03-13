#' Extract Counts Matrix from Molecule Experiment based on Assay Type
#'
#' @description 
#' This function retrieves a counts matrix from a Molecule Experiment object based on the given assay type.
#'
#' @param me A Molecule Experiment object.
#' @param assayName A character string indicating the assay type. Supported values include 
#' "sub-sector", "sub-concentric", "sub-combo", "super-concentric", and "super-combo".
#'
#' @return A counts matrix corresponding to the specified assay type.
#' @export
#' @importFrom MoleculeExperiment countMolecules
#' @examples
#' \dontrun{
#' # Assuming `data_obj` is your Molecule Experiment object
#' cm = CountsMatrix(data_obj, assayName = "sub-sector")
#' }

CountsMatrix <- function(me, assayName, ...) {
  # Assuming 'me' is globally available or should pass it as an argument
  counts_matrix = MoleculeExperiment::countMolecules(me, moleculesAssay = "detected", boundariesAssay = assayName, matrixOnly = TRUE, ...)

  # Determine which method to use for counts matrix transformation
  if (assayName == "sub_sector") {
    return(counts_matrix)
  } else if (assayName == "sub_concentric") {
    return(annuli_counts(counts_matrix))
  } else if (assayName == "super_sector") {
    return(counts_matrix)
  } else if (assayName == "super_concentric") {
    # Convert to annuli counts and then remove inner polygons
    out_concentric_super <- annuli_counts(counts_matrix)
    return(delete_inner(out_concentric_super))
  } else {
    stop("Unknown assay type!")
  }
}