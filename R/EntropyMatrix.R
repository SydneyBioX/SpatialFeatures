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

EntropyMatrix <- function(me, assayNames, ...) {
  # Ensure assayNames are valid
  if (!all(assayNames %in% c("sub_sector", "sub_concentric", 'sub_combo', "super_sector", "super_concentric", "super_combo"))) {
    stop("Invalid assayName(s) provided!")
  }
  
  # Generate a list of entropy matrices
  entropy_matrices <- lapply(assayNames, function(assay) {
    counts_matrix <- CountsMatrix(me, assay, ...)
    matrix_entropy(counts_matrix)
  })
  
  # Name the list based on the assay names
  names(entropy_matrices) <- assayNames
  
  return(entropy_matrices)
}