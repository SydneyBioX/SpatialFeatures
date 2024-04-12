# ==============================================================================
# helper functions for MCC main functions
# ==============================================================================

#' Extract Cell Boundaries from SpatialExperiment Object
#'
#' @param me A SpatialExperiment object containing spatial and molecular data.
#' @return A data frame with cell boundaries extracted from the SpatialExperiment object.
#' @export
extract_boundaries <- function(me) {
  df_boundary <- data.frame(MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE))
  return(df_boundary)
}

#' Determine if Cells are Neighbors
#'
#' @param df1 A data frame with x and y coordinates of the first set of cells.
#' @param df2 A data frame with x and y coordinates of the second set of cells.
#' @param threshold A numeric value specifying the distance threshold to consider cells as neighbors.
#' @return TRUE if any cell in df1 is a neighbor of any cell in df2 within the specified threshold, otherwise FALSE.
#' @export
are_neighbors <- function(df1, df2, threshold) {
  nearest_pairs <- nn2(data = df1[, c("x_location", "y_location")], 
                       query = df2[, c("x_location", "y_location")], 
                       k = 1, 
                       searchtype = "standard")
  
  if(any(nearest_pairs$nn.dists < threshold)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Get Neighbors for a Given Cell
#'
#' @param cell_id The identifier of the cell to find neighbors for.
#' @param df A data frame with cell information including 'segment_id' and 'sample_id'.
#' @param threshold Distance threshold to consider cells as neighbors (default is 2).
#' @return A vector of segment IDs representing the neighbors of the specified cell.
#' @export
get_neighbors <- function(cell_id, df, threshold = 2) {
  cell_data <- df %>% filter(segment_id == cell_id)
  potential_neighbors <- df %>% filter(sample_id == cell_data$sample_id[1] & segment_id != cell_id)
  neighbor_ids <- unique(potential_neighbors$segment_id)
  
  true_neighbors <- sapply(neighbor_ids, function(id) {
    neighbor_data <- df %>% filter(segment_id == id)
    are_neighbors(cell_data, neighbor_data, threshold)
  })
  
  neighbor_ids[true_neighbors]
}

#' Get All Neighbors in a Dataset
#'
#' @param df A data frame with cell information, including 'segment_id' and 'sample_id'.
#' @param threshold Distance threshold to consider cells as neighbors (default is 2).
#' @return A data frame listing pairs of cells and their neighbors.
#' @export
get_all_neighbors <- function(df, threshold = 2) {
  all_cells <- unique(df$segment_id)
  results <- list()
  
  for(cell_id in all_cells) {
    cell_sample_id <- df %>% filter(segment_id == cell_id) %>% select(sample_id) %>% distinct() %>% pull(sample_id)
    neighbors <- get_neighbors(cell_id, df, threshold)
    
    for(neighbor in neighbors) {
      neighbor_sample_id <- df %>% filter(segment_id == neighbor) %>% select(sample_id) %>% distinct() %>% pull(sample_id)
      results[[length(results) + 1]] <- data.frame(cell = paste(cell_sample_id, cell_id, sep = "."), 
                                                   neighbor = paste(neighbor_sample_id, neighbor, sep = "."))
    }
  }
  
  do.call(rbind, results)
}

#' Compute Matthews Correlation Coefficient
#'
#' @param table A 2x2 contingency table.
#' @return The Matthews correlation coefficient calculated from the contingency table.
#' @export
compute_mcc <- function(table) {
  a <- table[1, 1]
  b <- table[1, 2]
  c <- table[2, 1]
  d <- table[2, 2]
  
  denominator <- sqrt((a + b) * (a + c) * (b + d) * (c + d))
  if (denominator == 0 || is.na(denominator)) {
    return(0)
  } else {
    mcc <- (a * d - b * c) / denominator
    return(mcc)
  }
}

#' Extract Data and Neighbors from SpatialExperiment
#'
#' @param me A SpatialExperiment object.
#' @return A list containing a data frame of boundaries and a data frame of cell neighbors.
#' @export
extract_data_and_neighbors <- function(me) {
  df <- extract_boundaries(me)
  cellneighbours <- get_all_neighbors(df)
  list(df = df, cellneighbours = cellneighbours)
}

#' Count Molecules in SpatialExperiment
#'
#' @param me A SpatialExperiment object.
#' @return An object with molecule counts for each cell.
#' @export
count_molecules_from_SE <- function(me) {
  countMolecules(
    me,
    moleculesAssay = "detected",
    boundariesAssay = "cell",
    buffer = 0,
    matrixOnly = FALSE,
    nCores = 1
  )
}

#' Calculate Sum of Gene Counts from Neighbors
#'
#' @param mat A matrix of gene counts.
#' @param cellneighbours A data frame or list with neighbor information.
#' @return A matrix with original gene counts and additional columns for neighbor sums.
#' @export
calculate_neighbour_sums <- function(mat, cellneighbours) {
  library(dplyr)
  # Function to get the neighbor IDs for a given cell ID
  get_cell_neighbours <- function(cell_id) {
    neighbours <- cellneighbours %>% 
      filter(cell == cell_id) %>%
      pull(neighbor)
    return(neighbours)
  }
  
  # Function to sum up gene counts from neighbors
  get_neighbour_counts <- function(cell_id) {
    neighbours <- get_cell_neighbours(cell_id)
    # Check if there are no valid neighbours
    if (length(neighbours) == 0 || all(!neighbours %in% colnames(mat))) {
      return(matrix(0, nrow = nrow(mat), ncol = 1))
    }
    
    # Subset the matrix, ensuring it does not drop dimensions
    neighbour_mat <- mat[, neighbours, drop = FALSE]
    
    # Force neighbour_mat to be a matrix in case it's not
    # This is crucial when there's only one neighbour resulting in a vector instead of a matrix
    if (is.vector(neighbour_mat)) {
      neighbour_mat <- matrix(neighbour_mat, nrow = nrow(mat))
    }
    
    # Additional check to ensure neighbour_mat is never a vector
    if (length(dim(neighbour_mat)) < 2) {
      neighbour_mat <- matrix(neighbour_mat, nrow = nrow(mat), ncol = 1)
    }
    
    # Calculate gene counts from neighbours
    gene_counts <- rowSums(neighbour_mat, na.rm = TRUE)
    
    return(matrix(gene_counts, nrow = nrow(mat), ncol = 1))
  }
  
  # Initialize an empty matrix to store neighbor sums
  neighbour_sums <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  colnames(neighbour_sums) <- paste0(colnames(mat), "_n")
  
  # Iterate over each cell and calculate neighbor sums
  for (cell in colnames(mat)) {
    neighbour_sums[, paste0(cell, "_n")] <- get_neighbour_counts(cell)
  }
  
  # Combine the original matrix with the neighbor sums
  combined_mat <- cbind(mat, neighbour_sums)
  return(combined_mat)
}

#' Order and Convert Matrix to Sparse Format
#'
#' @param mat A matrix to be ordered and converted.
#' @return A sparse matrix in ordered format.
#' @export
order_and_sparse_matrix <- function(mat) {
  base_names <- gsub("_n$", "", colnames(mat))
  ordered_indices <- as.vector(rbind(match(base_names, colnames(mat)), 
                                     match(paste0(base_names, "_n"), colnames(mat))))
  combined_mat <- mat[, ordered_indices]
  mat_sparse_ordered <- as(combined_mat, "dgCMatrix")
  mat_sparse_ordered
}

#' Compute MCC Matrix for Spatial Cell Data
#'
#' @param mat A matrix with gene counts for cells and their neighbors.
#' @return A matrix with Matthews correlation coefficients for gene pairs across cells.
#' @export
compute_mcc_matrix <- function(mat) {
  gene_names <- rownames(mat)
  cell_names <- colnames(mat)[seq(1, ncol(mat), by = 2)]
  
  mcc_matrix <- matrix(0, nrow = choose(length(gene_names), 2), ncol = length(cell_names))
  rownames(mcc_matrix) <- combn(gene_names, 2, function(x) paste(x, collapse = "_"))
  colnames(mcc_matrix) <- cell_names
  
  for (cell_index in seq_along(cell_names)) {
    cell <- cell_names[cell_index]
    neighbor <- paste0(cell, "_n")
    
    mcc_values <- combn(gene_names, 2, function(gene_pair) {
      gene1_index <- which(rownames(mat) == gene_pair[1])
      gene2_index <- which(rownames(mat) == gene_pair[2])
      
      # Construct the table based on expression counts or presence
      table <- matrix(c(mat[gene1_index, cell_index * 2 - 1], mat[gene1_index, cell_index * 2],
                        mat[gene2_index, cell_index * 2 - 1], mat[gene2_index, cell_index * 2]),
                      nrow = 2, byrow = TRUE)
      
      mcc <- compute_mcc(table)
      return(mcc)
    })
    
    mcc_matrix[, cell_index] <- mcc_values
  }
  
  return(mcc_matrix)
}