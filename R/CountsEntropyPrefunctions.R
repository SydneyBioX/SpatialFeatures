#' Convert Cumulative Counts in Concentric Polygons to Counts in Annuli Polygons
#' 
#' Converts cumulative counts of concentric polygons to actual counts in annuli polygons.
#' 
#' @param mat Matrix containing cumulative counts.
#' 
#' @return A matrix with updated counts.
annuli_counts <- function(mat) {
  # Find unique cell prefixes
  cell_prefixes <- unique(sub("_[0-9]+$", "", colnames(mat)))
  
  for(prefix in cell_prefixes) {
    # Identify columns associated with the current cell prefix
    cols <- grep(paste0("^", prefix, "_"), colnames(mat))
    
    # If there's only one column for this prefix, skip it
    if(length(cols) < 2) next
    
    # For every column except the first, subtract the previous column
    for(i in length(cols):2) {
      mat[,cols[i]] <- mat[,cols[i]] - mat[,cols[i-1]]
    }
  }
  return(mat)
}

#' Convert Cumulative Counts in Concentric Polygons to Counts in Combo Polygons
#' 
#' Converts cumulative counts of concentric polygons to actual counts in combo polygons.
#' 
#' @param mat Matrix containing cumulative counts.
#' 
#' @return A matrix with updated counts.
combo_counts <- function(mat) {
  # Find unique cell prefixes
  cell_prefixes <- unique(gsub("(_[0-9]+)_[0-9]+$", "\\1", colnames(mat)))
  for(prefix in cell_prefixes) {
    # Identify columns associated with the current cell prefix
    cols <- grep(paste0("^", prefix, "_"), colnames(mat))
    
    # If there's only one column for this prefix, skip it
    if(length(cols) < 2) next
    
    # For every column except the first, subtract the previous column
    for(i in length(cols):2) {
      mat[,cols[i]] <- mat[,cols[i]] - mat[,cols[i-1]]
    }
  }
  return(mat)
}

#' Delete Innermost Columns
#' 
#' Deletes columns ending with '0', which represent the innermost polygon.
#' 
#' @param m Matrix containing counts.
#' 
#' @return A matrix without innermost polygon columns.
delete_inner <- function(m) {
  # Identify columns that don't end with "_0"
  cols_to_keep <- !grepl("_0$", colnames(m))
  
  # Subset the matrix to keep only these columns
  m_sub <- m[, cols_to_keep]
  return(m_sub)
}

#' Calculate Entropy
#' 
#' Calculates the entropy for a given probability distribution.
#' 
#' @param p Numeric vector representing a probability distribution.
#' 
#' @return The entropy of the distribution.
entropy <- function(p) {
  p <- p[p > 0] # remove zeroes
  -sum(p * log(p))
}

#' Get Unique Segments
#' 
#' Extracts unique segments from the column names of a matrix.
#' 
#' @param mat Matrix from which to extract segment names.
#' 
#' @return A character vector of unique segments.
get_segments <- function(mat) {
  unique(gsub("_.+$", "", colnames(mat)))
}

#' Extract Sub Matrix
#' 
#' Extracts a sub-matrix corresponding to a specific segment from a larger matrix.
#' 
#' @param mat Matrix from which to extract the sub-matrix.
#' @param segment Character string representing the segment of interest.
#' 
#' @return A matrix corresponding to the specified segment.
extract_sub_mat <- function(mat, segment) {
  matching_cols <- grepl(segment, colnames(mat))
  mat[, matching_cols, drop = FALSE]
}

#' Compute Entropy for Matrix
#' 
#' Computes entropy for each row of a matrix.
#' 
#' @param mat Matrix for which to compute entropy.
#' 
#' @return A data frame of entropy values for each row of the matrix.
compute_entropy <- function(mat) {
  # Convert to proportions
  row_sums <- rowSums(mat)
  proportions <- t(apply(mat, 1, function(row) {
    if (sum(row) == 0) {
      return(rep(0, length(row)))
    } else {
      return(row / sum(row))
    }
  }))
  apply(proportions, 1, entropy)
}
