#' Calculate the Entropy of a Counts Matrix
#'
#' @description
#' This function computes the entropy of a counts matrix from a Molecule Experiment object based on the given assay type.
#'
#' @param me A Molecule Experiment object.
#' @param featureTypes A character string specifying the feature type. Supported values include
#' "sub_sector", "sub_concentric", "super_sector", and "super_concentric".
#' @param nCores Number of cores
#' @param ... arguments passing to CountsMatrix
#'
#' @return A matrix representing the entropy values corresponding to the given assay type.
#' @export
#' @examples
#' data(example_me)
#' me <- loadBoundaries(me)
#' ent <- EntropyMatrix(me, c("sub_sector", "sub_concentric", "super_sector", "super_concentric"), nCores = 1)
#' lapply(ent, head, n = 4)
EntropyMatrix <- function(me, featureTypes, nCores = 1, ...) {
  # Ensure featureTypes are valid
  # if (!all(featureTypes %in% c("sub_sector", "sub_concentric", 'sub_combo', "super_sector", "super_concentric", "super_combo"))) {
  #   stop("Invalid assayName(s) provided!")
  # }

  # Generate a list of entropy matrices
  entropy_matrices <- lapply(featureTypes, function(assay) {
    counts_matrix <- CountsMatrix(me = me, assayName = assay, nCores = nCores, ...)
    matrix_entropy(counts_matrix, nCores = nCores)
  })

  # Name the list based on the assay names
  names(entropy_matrices) <- featureTypes

  return(entropy_matrices)
}
