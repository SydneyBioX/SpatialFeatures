# ==============================================================================
# helper functions for MCC main functions
# ==============================================================================

#' Extract boundaries from a SpatialExperiment object.
#' 
#' @param me A SpatialExperiment object.
#' @return A dataframe containing boundaries.
#' @export
extract_boundaries <- function(me) {
  df_boundary <- data.frame(boundaries(me, assayName = "cell", flatten = TRUE))
  return(df_boundary)
}

#' Check if two cells are neighbors based on their boundary points.
#' 
#' @param df1 Dataframe of the first cell's boundaries.
#' @param df2 Dataframe of the second cell's boundaries.
#' @param threshold Distance threshold to determine if cells are neighbors.
#' @return Logical indicating if the cells are neighbors.
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

#' Get neighbors for a given cell in a dataframe.
#' 
#' @param cell_id ID of the cell to find neighbors for.
#' @param df Dataframe containing cell data.
#' @param threshold Distance threshold to determine if cells are neighbors.
#' @return A vector containing IDs of the neighboring cells.
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

#' Get neighbors for all cells and save to dataframe.
#' 
#' @param df Dataframe containing cell data.
#' @param threshold Distance threshold to determine if cells are neighbors.
#' @return A dataframe containing cells and their respective neighbors.
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

#' Compute Matthews Correlation Coefficient (MCC) for a given 2x2 table.
#' 
#' @param table A 2x2 matrix or table.
#' @return The computed MCC value.
#' @export
compute_mcc <- function(table) {
  a <- table[1, 1]
  b <- table[1, 2]
  c <- table[2, 1]
  d <- table[2, 2]
  
  denominator <- sqrt((a+b)*(a+c)*(b+d)*(c+d))
  if (denominator == 0) {
    return(0)
  } else {
    return((a*d - b*c) / denominator)
  }
}

#' Extract boundaries and calculate neighbors for a SpatialExperiment object.
#' 
#' @param me A SpatialExperiment object.
#' @return A list containing boundary data and neighbors.
#' @export
extract_data_and_neighbors <- function(me) {
  df <- extract_boundaries(me)
  cellneighbours <- get_all_neighbors(df)
  list(df = df, cellneighbours = cellneighbours)
}

#' Count molecules in a SpatialExperiment object.
#' 
#' @param me A SpatialExperiment object.
#' @return A cell object with counted molecules.
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

#' Calculate neighbor sums for a given matrix.
#' 
#' @param mat A matrix of counted molecules.
#' @param cellneighbours A dataframe containing cells and their respective neighbors.
#' @return A matrix with neighbor sums.
#' @export
calculate_neighbour_sums <- function(mat, cellneighbours) {
  # Extract neighbors from cellneighbours for each cell
  get_cell_neighbours <- function(cell_id) {
    neighbours <- cellneighbours %>% 
      filter(cell == cell_id) %>%
      pull(neighbor)
    return(neighbours)
  }
  
  # Sum up gene counts from neighbors
  get_neighbour_counts <- function(cell_id) {
    neighbours <- get_cell_neighbours(cell_id)
    if(length(neighbours) == 0 || all(!neighbours %in% colnames(mat))) {
      return(matrix(0, nrow=nrow(mat), ncol=1))
    }
    neighbour_mat <- mat[, neighbours, drop=FALSE]
    gene_counts <- rowSums(neighbour_mat, na.rm = TRUE)
    return(matrix(gene_counts, nrow=length(gene_counts), ncol=1))
  }
  
  neighbour_sums <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  colnames(neighbour_sums) <- paste0(colnames(mat), "_n")
  for(cell in colnames(mat)) {
    neighbour_sums[, paste0(cell, "_n")] <- get_neighbour_counts(cell)
  }
  
  combined_mat <- cbind(mat, neighbour_sums)
  combined_mat
}

#' Order the matrix and convert to sparse format.
#' 
#' @param mat A matrix to be ordered and converted.
#' @return A sparse ordered matrix.
#' @export
order_and_sparse_matrix <- function(mat) {
  base_names <- gsub("_n$", "", colnames(mat))
  ordered_indices <- as.vector(rbind(match(base_names, colnames(mat)), 
                                     match(paste0(base_names, "_n"), colnames(mat))))
  combined_mat <- mat[, ordered_indices]
  mat_sparse_ordered <- as(combined_mat, "dgCMatrix")
  mat_sparse_ordered
}

#' Calculate MCC for all pairs of genes for a given cell and its neighbors.
#' 
#' @param cell_col Column of a matrix representing a cell.
#' @param neighbor_col Column of a matrix representing a cell's neighbor.
#' @return A vector containing MCC values for gene pairs.
#' @export
calculate_mcc_for_cell <- function(cell_col, neighbor_col) {
  # Generating all ordered pairs of genes
  gene_pairs <- expand.grid(rownames(mat_sparse), rownames(mat_sparse))
  
  # Removing pairs where both genes are the same
  gene_pairs <- gene_pairs[gene_pairs$Var1 != gene_pairs$Var2, ]
  
  sapply(1:nrow(gene_pairs), function(i) {
    gene1 <- gene_pairs$Var1[i]
    gene2 <- gene_pairs$Var2[i]
    
    table <- matrix(0, 2, 2)
    table[1, 1] <- cell_col[gene1, ]
    table[1, 2] <- neighbor_col[gene1, ]
    table[2, 1] <- cell_col[gene2, ]
    table[2, 2] <- neighbor_col[gene2, ]
    
    compute_mcc(table)
  }) %>% 
    set_names(paste("mcc", gene_pairs$Var1, gene_pairs$Var2, sep="_"))
}
