#' Comprehensive Function to Compute MCC Matrix from SpatialExperiment
#'
#' This function wraps several processing steps including extracting cell boundaries,
#' counting molecules, calculating neighbor sums, ordering and sparsifying the matrix,
#' and computing the MCC matrix.
#'
#' @param me A SpatialExperiment object.
#' @return A sparse matrix of Matthews correlation coefficients.
#' @export
#' @examples
#' # Assuming 'small_me' is a SpatialExperiment object
#' # not run
#' # data(small_me, package = "MoleculeExperiment")
#' # mcc_matrix <- MCCMatrix(small_me)
MCCMatrix <- function(me) {
  # Step 1: Extract boundaries and calculate neighbors
  results <- SpatialFeatures::extract_data_and_neighbors(me)
  # Step 2: Count molecules in the SpatialExperiment object
  cell <- SpatialFeatures::count_molecules_from_SE(me)
  mat <- assay(cell, "counts")
  # Step 3: Calculate neighbor sums
  combined_mat <- calculate_neighbour_sums(mat, results$cellneighbours)
  # Step 4: Order and convert to sparse matrix
  mat_sparse_ordered <- SpatialFeatures::order_and_sparse_matrix(combined_mat)
  # Step 5: Compute the MCC matrix
  mcc_matrix <- SpatialFeatures::compute_mcc_matrix(mat_sparse_ordered)
  return(mcc_matrix)
}
