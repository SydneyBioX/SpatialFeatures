#' Compute the MCC matrix for spatial cell data.
#' 
#' @param me A SpatialExperiment object.
#' @return A sparse matrix with MCC calculations.
#' @export
MCCMatrix <- function(me) {
  # Step 1: Extract boundaries and calculate neighbors
  results <- extract_data_and_neighbors(me)
  
  # Step 2: Count molecules in the SpatialExperiment object
  cell <- count_molecules_from_SE(me)
  mat <- counts(cell)
  
  # Step 3: Calculate neighbor sums
  combined_mat <- calculate_neighbour_sums(mat, results$cellneighbours)
  
  # Step 4: Order the matrix and convert to sparse format
  mat_sparse_ordered <- order_and_sparse_matrix(combined_mat)
  
  # Step 5: Compute the MCC for all pairs of genes for a given cell and its neighbors
  # Assuming 'mat_sparse_ordered' is the matrix you want to use
  mcc_results <- mclapply(seq(1, ncol(mat_sparse_ordered), by=2), function(i) {
    cell_col <- mat_sparse_ordered[, i, drop=FALSE]
    neighbor_col <- mat_sparse_ordered[, i+1, drop=FALSE]
    calculate_mcc_for_cell(cell_col, neighbor_col)
  }, mc.cores = 10) # This uses parallel processing
  
  mcc_df <- do.call(cbind, mcc_results)
  
  colnames(mcc_df) <- colnames(mat_sparse_ordered)[seq(1, ncol(mat_sparse_ordered), by=2)]
  
  return(mcc_df)
}