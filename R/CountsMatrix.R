#' Compute Counts Matrix for Specified Assay Name
#'
#' This function calculates the counts matrix for the specified assay name.
#' 
#' @param me A data object that contains the input data.
#' @param assayName A character string specifying the assay name. It must be one 
#'   of "sub-sector", "sub-concentric", "sub-combo", "super-concentric", or 
#'   "super-combo".
#'
#' @return A matrix containing the counts data.
#' @export
CountsMatrix <- function(me, assayName) {
  # me = loadBoundariesToME(me)
  # Ensure assayName is valid
  if (!assayName %in% c("sub-sector", "sub-concentric", 'sub-combo', 
                        "super-concentric", "super-combo")) {
    stop("Invalid assayName provided!")
  }
  
  # Assuming 'me' is globally available or you should pass it as an argument
  spe = countMolecules(me, moleculesAssay = "detected", boundariesAssay = assayName, buffer = 0, matrixOnly = FALSE, nCores = 1)
  
  # Extract the counts matrix from the spe object
  counts_matrix <- counts(spe)
  
  # Determine which method to use for counts matrix transformation
  if (assayName == "sub-sector") {
    return(counts_matrix)
  } else if (assayName == "sub-concentric") {
    return(annuli_counts(counts_matrix))
  } else if (assayName == "sub-combo") {
    return(combo_counts(counts_matrix))
  } else if (assayName == "super-concentric") {
    # Convert to annuli counts and then remove inner polygons
    out_concentric_super <- annuli_counts(counts_matrix)
    return(delete_inner(out_concentric_super))
  } else if (assayName == "super-combo") {
    # Convert to combo counts and then remove inner polygons
    out_combo_super <- combo_counts(counts_matrix)
    return(delete_inner(out_combo_super))
  } else {
    stop("Unknown assay type!")
  }
}