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
CountsMatrix <- function(assayName, counts_matrix) {
  # Ensure assayName is valid
  if (!assayName %in% c("sub_sector", "sub_concentric", 'sub_combo', "super_sector", "super_concentric", "super_combo")) {
    stop("Invalid assayName provided!")
  }

  # Determine which method to use for counts matrix transformation
  if (assayName == "sub_sector") {
    return(counts_matrix)
  } else if (assayName == "sub_concentric") {
    return(annuli_counts(counts_matrix))
  } else if (assayName == "sub_combo") {
    return(combo_counts(counts_matrix))
  } else if (assayName == "super_sector") {
    return(counts_matrix)
  } else if (assayName == "super_concentric") {
    # Convert to annuli counts and then remove inner polygons
    out_concentric_super <- annuli_counts(counts_matrix)
    return(delete_inner(out_concentric_super))
  } else if (assayName == "super_combo") {
    # Convert to combo counts and then remove inner polygons
    out_combo_super <- combo_counts(counts_matrix)
    return(delete_inner(out_combo_super))
  } else {
    stop("Unknown assay type!")
  }
}