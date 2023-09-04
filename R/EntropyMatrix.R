#' Compute Entropy Matrix for Specified Assay Name
#'
#' This function calculates the entropy matrix for the specified assay name.
#' 
#' @param me A data object that contains the input data.
#' @param assayName A character string specifying the assay name. It must be one 
#'   of "sub-sector", "sub-concentric", "sub-combo", "super-concentric", or 
#'   "super-combo".
#'
#' @return A matrix containing the entropy data.
#' @export
EntropyMatrix <- function(me, assayName) {
  # Ensure assayName is valid
  if (!assayName %in% c("sub-sector", "sub-concentric", 'sub-combo', 
                        "super-concentric", "super-combo")) {
    stop("Invalid assayName provided!")
  }
  # Extract the counts matrix
  counts_matrix <- CountsMatrix(me, assayName)
  entropy_matrix <- matrix_entropy(counts_matrix)
  return(entropy_matrix)
}