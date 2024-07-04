# ==============================================================================
# helper functions for Entropy main functions
# ==============================================================================

# Required Libraries

#' Extract boundaries from a given object
#'
#' @param me The object to extract boundaries from.
#' @return A data frame of extracted boundaries.
extract_boundaries <- function(me) {
  df_boundary <- data.frame(MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE))
  return(df_boundary)
}

#' Extract boundaries and calculate centroids
#'
#' @param me The object to process.
#' @return A list containing a data frame of boundaries and a data frame of centroids.
#' @export
extract_boundaries_and_centroids <- function(me) {
  df_boundary <- extract_boundaries(me)
  df_circle <- calculate_centroids(df_boundary)
  return(list(df_boundary = df_boundary, df_circle = df_circle))
}


#' Generate scale factors for all values of k
#'
#' @param k A numeric value.
#' @return A vector of scale factors.
generate_scale_factors_all <- function(k) {
  scale_factors <- seq(from = 1/(k+1), to = 1, by = 1/(k+1))
  return(scale_factors)
}

#' Compute entropy for the entire matrix
#' @export
#' @param mat A matrix to compute entropy for.
#' @return A data frame of entropy values for the matrix.
matrix_entropy <- function(mat, nCores = 1) {
  mat <- as.matrix(mat)  # Convert dgCMatrix to a regular matrix if necessary
  cells <- get_cells(mat)

  results <- mclapply(cells, function(cell) {
    cell_segments_mat <- extract_cell_segments(mat, cell)
    compute_cell_entropy(cell_segments_mat)
  }, mc.cores = nCores)

  df_entropy <- as.data.frame(do.call(cbind, results))
  rownames(df_entropy) <- rownames(mat)
  colnames(df_entropy) <- cells

  return(df_entropy)
}





#' Split Molecule Experiment Object into Regions
#'
#' This function splits a given Molecule Experiment (ME) object into
#' 12 different rectangular regions. The regions are determined based on
#' the central coordinates and extents of the ME object. The function
#' returns a dataframe that associates each segment with a region.
#'
#' @param me A Molecule Experiment (ME) object.
#'
#' @return A dataframe with three columns: region (numeric region identifier
#'         ranging from 1 to 12), sample_id (unique identifier for each sample),
#'         and segment_ids (comma-separated segment IDs that belong to each
#'         region-sample combination). If there are no segments for a particular
#'         region-sample combination, segment_ids will be an empty string.
#'
#' @examples
#' \dontrun{
#' # Assuming `me_obj` is your Molecule Experiment object
#' regions_df <- split_rectangle(me_obj)
#' }
#'
#' @export
split_rectangle <- function(me) {
  # Compute central coordinates
  central_coords <- compute_central_coordinates(me)

  # Extract extent
  xmin = extent(me, assayName = "detected")[1]
  xmax = extent(me, assayName = "detected")[2]
  ymin = extent(me, assayName = "detected")[3]
  ymax = extent(me, assayName = "detected")[4]

  # Determine midpoints
  xmid = xmin + (xmax - xmin) / 2
  ymid1 = ymin + (ymax - ymin) / 3  # One-third the height for top row boundary
  ymid2 = ymin + 2 * (ymax - ymin) / 3  # Two-thirds the height for bottom row boundary

  # Determine which region each segment belongs to
  central_coords$region <-
    case_when(
      central_coords$x_central < xmid & central_coords$y_central > 2*ymid2 ~ 1,
      central_coords$x_central < xmid & central_coords$y_central > ymid2 ~ 2,
      central_coords$x_central < xmid & central_coords$y_central > ymid1 ~ 3,
      central_coords$x_central < xmid ~ 4,
      central_coords$x_central >= xmid & central_coords$y_central > 2*ymid2 ~ 5,
      central_coords$x_central >= xmid & central_coords$y_central > ymid2 ~ 6,
      central_coords$x_central >= xmid & central_coords$y_central > ymid1 ~ 7,
      central_coords$x_central >= xmid ~ 8,
      central_coords$x_central < xmid & central_coords$y_central <= ymid1 ~ 9,
      central_coords$x_central >= xmid & central_coords$y_central <= ymid1 ~ 10,
      central_coords$x_central < xmid & central_coords$y_central <= ymid2 ~ 11,
      central_coords$x_central >= xmid & central_coords$y_central <= ymid2 ~ 12
    )

  # Group by region and summarize
  result <- central_coords %>%
    group_by(region) %>%
    summarise(
      sample_ids = paste(unique(sample_id), collapse = ", "),
      segment_ids = paste(unique(segment_id), collapse = ", "),
      .groups = "drop"
    )

  # Create a reference dataframe for all regions
  all_regions <- data.frame(region = 1:12)

  # Left join with the reference dataframe to ensure all regions are included
  result <- left_join(all_regions, result, by = "region")

  # If segments or sample_ids column is NA, replace with empty string
  result$segment_ids[is.na(result$segment_ids)] <- ""
  result$sample_ids[is.na(result$sample_ids)] <- ""

  return(result)
}
