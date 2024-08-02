#' Convert cumulative counts in concentric polygons to annuli polygon counts
#'
#' @param mat A matrix of cumulative counts.
#'
#' @importFrom terra t diff
#' @return A matrix of annuli polygon counts.
annuli_counts = function(mat) {
  fac = sub("_[0-9]+$", "", colnames(mat))
  # tmat = terra::t(mat)
  tmat = terra::t(mat)
  tmat_split = split.data.frame(tmat, fac)
  # tmat_split_diffs = lapply(tmat_split, terra::diff)
  tmat_split_diffs = lapply(tmat_split, terra::diff)
  diffs = do.call(rbind, tmat_split_diffs)
  cts_new = terra::t(diffs)
  cts_new[cts_new < 0] <- 0
  return(cts_new)
}

#' Extract Counts Matrix from Molecule Experiment based on Assay Type
#'
#' @description
#' This function retrieves a counts matrix from a Molecule Experiment object based on the given assay type.
#'
#' @param me A Molecule Experiment object.
#' @param assayName A character string indicating the assay type. Supported values include
#' "sub-sector", "sub-concentric", "sub-combo", "super-concentric", and "super-combo".
#' @param nCores Number of cores for parallel processing (default 1)
#' @param ... arguments passing to MoleculeExperiment::countMolecules
#'
#' @return A counts matrix corresponding to the specified assay type.
#' @importFrom MoleculeExperiment countMolecules
#' @examples
#' \dontrun{
#' # Assuming `data_obj` is your Molecule Experiment object
#' cm = CountsMatrix(data_obj, assayName = "sub-sector")
#' }
CountsMatrix <- function(me, assayName, nCores = 1, ...) {
  counts_matrix = MoleculeExperiment::countMolecules(me, moleculesAssay = "detected", boundariesAssay = assayName, matrixOnly = TRUE, nCores = nCores, ...)

  if (assayName %in% c("subsector", "supersector")) {
    return(counts_matrix)
  }

  if (assayName %in% c("subconcentric", "superconcentric")) {
    return(annuli_counts(counts_matrix))
  }

  stop("Unknown assay type!")
}


#' Extract unique segments from matrix column names
#'
#' @param mat A matrix to process.
#' @return A character vector of unique segments.
get_cells <- function(mat) {
  # Extract the unique cell identifiers before the last underscore
  cells <- unique(gsub("(.*)_[^_]+$", "\\1", colnames(mat)))
  return(cells)
}

#' Extract a submatrix corresponding to a segment from a matrix
#'
#' @param mat A matrix to process.
#' @param cell A character string indicating the segment to extract.
#' @return A submatrix of the provided matrix.
extract_cell_segments <- function(mat, cell) {
  # Create a regex pattern to match all segments of the cell
  pattern <- paste0("^", cell, "_\\d+$")
  matching_cols <- grepl(pattern, colnames(mat))
  cell_segments_mat <- mat[, matching_cols, drop = FALSE]
  return(cell_segments_mat)
}

#' Calculate entropy from a probability vector
#'
#' @param p A numeric vector of probabilities.
#' @return The calculated entropy value.
calculate_entropy <- function(p) {
  p <- p[p > 0] # remove zeroes
  return(-sum(p * log(p)))
}

#' Calculate prop for entropy
#'
#' @param counts a vector of integers
#' @return the proportion
countsprop = function(counts) {
  total <- sum(counts)
  if (total == 0) {
    return(rep(0, length(counts)))
  } else {
    return(counts / total)
  }
}

#' Compute entropy for a given matrix
#'
#' @param cell_segments_mat A matrix to compute entropy for.
#' @return A vector of entropy values.
compute_cell_entropy <- function(cell_segments_mat) {
  # Convert to proportions
  proportions <- apply(cell_segments_mat, 1, countsprop)

  # Compute entropy for each gene
  gene_entropies <- apply(proportions, 2, calculate_entropy)
  return(gene_entropies)
}

#' Compute entropy for the entire matrix
#' @param mat A matrix to compute entropy for.
#' @param nCores Number of cores
#' @return A data frame of entropy values for the matrix.
matrix_entropy <- function(mat, nCores = 1) {
  mat <- as.matrix(mat)  # Convert dgCMatrix to a regular matrix if necessary
  cells <- get_cells(mat)

  results <- parallel::mclapply(cells, function(cell) {
    cell_segments_mat <- extract_cell_segments(mat, cell)
    compute_cell_entropy(cell_segments_mat)
  }, mc.cores = nCores)

  df_entropy <- as.data.frame(do.call(cbind, results))
  rownames(df_entropy) <- rownames(mat)
  colnames(df_entropy) <- cells

  return(df_entropy)
}
