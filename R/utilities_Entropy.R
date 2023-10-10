# ==============================================================================
# helper functions for Entropy main functions
# ==============================================================================

# Required Libraries
#' @importFrom dplyr %>% group_by mutate ungroup arrange
#' @importFrom tidyr group_split
#' @importFrom stringr str_extract
#' @importFrom stats entropy

#' Extract boundaries from a given object
#'
#' @param me The object to extract boundaries from.
#' @return A data frame of extracted boundaries.
extract_boundaries <- function(me) {
  df_boundary <- data.frame(MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE))
  return(df_boundary)
}

#' Calculate centroids from boundaries
#'
#' @param df_boundary A data frame of boundaries.
#' @return A data frame of centroids.
calculate_centroids <- function(df_boundary) {
  df_circle <- df_boundary %>%
    group_by(segment_id) %>%
    mutate(x_central = mean(head(x_location, -1)), 
           y_central = mean(head(y_location, -1))) %>%
    ungroup()
  return(df_circle)
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

#' Create sectors from a given data frame
#'
#' @param df A data frame to create sectors from.
#' @return A data frame of sectors.
create_sectors <- function(df) {
  df <- df[order(df$angle),]
  sectors <- list()
  
  for (i in 1:nrow(df)) {
    p1 <- df[i,]
    p2 <- if (i == nrow(df)) df[1,] else df[i+1,]
    
    sector <- data.frame(
      x_location_sector = c(p1$x_location, p2$x_location, p1$x_central, p1$x_location),
      y_location_sector = c(p1$y_location, p2$y_location, p1$y_central, p1$y_location),
      segment_id = p1$segment_id,
      sample_id = p1$sample_id,
      sector_id = paste0(p1$segment_id, "_", sprintf("%02d", i))
    )
    sectors[[i]] <- sector
  }
  
  do.call(rbind, sectors)
}

#' Generate scale factors for all values of k
#'
#' @param k A numeric value.
#' @return A vector of scale factors.
generate_scale_factors_all <- function(k) {
  scale_factors <- seq(from = 1/(k+1), to = 1, by = 1/(k+1))
  return(scale_factors)
}

#' Modify sectors based on a given data frame
#'
#' @param df A data frame to process.
#' @return A modified data frame of sectors.
create_sectors_modified <- function(df) {
  df <- df[order(df$angle),]
  sectors <- list()
  
  concentric_order <- as.numeric(str_extract(df$concentric_id[1], "(?<=_)[0-9]+$"))
  
  for (i in 1:nrow(df)) {
    p1 <- df[i,]
    p2 <- if (i == nrow(df)) df[1,] else df[i+1,]
    
    sector <- data.frame(
      x_location_sector = c(p1$x_location, p2$x_location, p1$x_central, p1$x_location),
      y_location_sector = c(p1$y_location, p2$y_location, p1$y_central, p1$y_location),
      segment_id = p1$segment_id,
      sample_id = p1$sample_id,
      combo_id = paste0(p1$segment_id, "_", sprintf("%02d", i), "_", concentric_order)
    )
    sectors[[i]] <- sector
  }
  
  do.call(rbind, sectors)
}

#' Create sectors for scaled values
#'
#' @param df A data frame to process.
#' @param scale A numeric scale factor.
#' @param k A numeric value, defaulting to 8.
#' @return A data frame of sectors for scaled values.
create_sectors_for_scaled <- function(df, scale, k = 8) {
  # Initial transformations
  df <- df %>%
    mutate(x_location = x_scaled, y_location = y_scaled) %>%
    group_by(segment_id) %>%
    mutate(angle = atan2(y_scaled - y_central, x_scaled - x_central) + pi)
  
  # Apply the create_sectors_modified function on each group
  sectors_list <- lapply(split(df, df$segment_id), create_sectors_modified)
  
  # Combine results, arrange and add concentric_id
  df_sectors <- do.call(rbind, sectors_list) %>%
    arrange(sample_id) %>%
    mutate(concentric_id = paste0(segment_id, '_', scale*(k+1)))
  
  return(df_sectors)
}

#' Generate scale factors for all outside values of ks
#'
#' @param ks A numeric value.
#' @return A vector of scale factors for outside values.
generate_scale_factors_all_outside <- function(ks) {
  scale_factors <- c(seq(from = 1, to = 2, by = 1/(ks+1)))
  return(scale_factors)
}

#' Create a scaled data frame
#'
#' @param s A numeric scaling factor.
#' @param df_circle A data frame to scale.
#' @param k A numeric value, defaulting to 8.
#' @return A scaled data frame.
create_scaled_df <- function(s, df_circle, k = 8) {
  scaling_factor <- s * (k + 1)
  
  df_circle %>%
    group_by(segment_id) %>%
    mutate(
      x_diff = x_location - x_central,
      y_diff = y_location - y_central,
      x_scaled = x_central + s * x_diff,
      y_scaled = y_central + s * y_diff,
      scale = s,
      concentric_id = paste0(segment_id, "_", scaling_factor)
    ) %>%
    select(-x_diff, -y_diff) # Remove temporary columns
}

#' Create sector data frame from circle data
#'
#' This function takes a data frame of circle data with segment IDs and locations, 
#' and returns a data frame of sectors with their scaled coordinates.
#'
#' @param df A data frame containing circle data with segment IDs, x and y locations.
#' @return A data frame with the coordinates for each sector.
create_sector_df <- function(df) {

  # Distinct boundary points
  df <- df %>%
    distinct(segment_id, x_location, y_location, .keep_all = TRUE) %>%
    group_by(segment_id) %>%
    mutate(
      x_diff = x_location - x_central,
      y_diff = y_location - y_central,
      x_scaled = x_central + 2 * x_diff,
      y_scaled = y_central + 2 * y_diff
    ) %>%
    ungroup()

  # Create sectors
  sectors <- list()
  df_grouped <- split(df, df$segment_id)

  for (segment in names(df_grouped)) {
    current_segment <- df_grouped[[segment]]
    n_rows <- nrow(current_segment)

    for (i in 1:(n_rows-1)) {  # We iterate till n_rows-2, to skip the last point
      p1 <- current_segment[i,]
      p2 <- current_segment[i+1,]

      sector <- data.frame(
        x_location_sector = c(p1$x_location, p2$x_location, p2$x_scaled, p1$x_scaled, p1$x_location),
        y_location_sector = c(p1$y_location, p2$y_location, p2$y_scaled, p1$y_scaled, p1$y_location),
        segment_id = p1$segment_id,
        sample_id = p1$sample_id,
        sector_id = paste0(p1$segment_id, "_", sprintf("%02d", i))
      )
      sectors[[length(sectors) + 1]] <- sector
    }

    # To account for the last boundary point connecting with the first, without repetition
    p1 <- current_segment[n_rows,]
    p2 <- current_segment[1,]

    sector <- data.frame(
      x_location_sector = c(p1$x_location, p2$x_location, p2$x_scaled, p1$x_scaled, p1$x_location),
      y_location_sector = c(p1$y_location, p2$y_location, p2$y_scaled, p1$y_scaled, p1$y_location),
      segment_id = p1$segment_id,
      sample_id = p1$sample_id,
      sector_id = paste0(p1$segment_id, "_", sprintf("%02d", n_rows-1))
    )
    sectors[[length(sectors) + 1]] <- sector
  }

  df_sectors <- do.call(rbind, sectors)
  # Rename columns as per request
  df_sectors <- df_sectors %>%
    dplyr::rename(
      x_section = x_location_sector,
      y_section = y_location_sector,
      area_id = sector_id
    )

  return(df_sectors)
}

#' Modify the area_id format
#'
#' This function takes an area_id string and modifies the last segment 
#' if it's a single-digit number between 1 and 9 by adding a leading zero.
#'
#' @param id A string representing the area_id.
#'
#' @return A string with the modified area_id format.
modify_area_id <- function(id) {
  # Use a regular expression to check if the pattern matches, and replace it if so
  modified_id <- sub("_(\\d)$", "_0\\1", id)
  return(modified_id)
}

#' Convert cumulative counts in concentric polygons to annuli polygon counts
#'
#' @param mat A matrix of cumulative counts.
#' @return A matrix of annuli polygon counts.
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

#' Convert cumulative counts in concentric polygons to combo polygon counts
#'
#' @param mat A matrix of cumulative counts.
#' @return A matrix of combo polygon counts.
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

#' Delete columns ending with '0' from a matrix
#'
#' @param m A matrix to process.
#' @return A modified matrix without columns ending in '0'.
delete_inner <- function(m) {
  # Identify columns that don't end with "_0"
  cols_to_keep <- !grepl("_09$", colnames(m))
  
  # Subset the matrix to keep only these columns
  m_sub <- m[, cols_to_keep]
  return(m_sub)
}

#' Calculate entropy from a probability vector
#'
#' @param p A numeric vector of probabilities.
#' @return The calculated entropy value.
calculate_entropy <- function(p) {
  p <- p[p > 0] # remove zeroes
  return(-sum(p * log(p)))
}

#' Extract unique segments from matrix column names
#'
#' @param mat A matrix to process.
#' @return A character vector of unique segments.
get_segments <- function(mat) {
  segments <- unique(gsub("_.+$", "", colnames(mat)))
  return(segments)
}

#' Extract a submatrix corresponding to a segment from a matrix
#'
#' @param mat A matrix to process.
#' @param segment A character string indicating the segment to extract.
#' @return A submatrix of the provided matrix.
extract_sub_mat <- function(mat, segment) {
  matching_cols <- grepl(segment, colnames(mat))
  sub_mat <- mat[, matching_cols, drop = FALSE]
  return(sub_mat)
}

#' Compute entropy for a given matrix
#'
#' @param mat A matrix to compute entropy for.
#' @return A vector of entropy values.
compute_entropy <- function(mat) {
  # Convert to proportions
  proportions <- t(apply(mat, 1, function(row) {
    total <- sum(row)  # Compute sum of the row
    if (total == 0) {
      return(rep(0, length(row)))
    } else {
      return(row / total)
    }
  }))
  out_mat = apply(proportions, 1, calculate_entropy)
  return(out_mat)
}

#' Compute entropy for the entire matrix
#'
#' @param mat A matrix to compute entropy for.
#' @return A data frame of entropy values for the matrix.
matrix_entropy <- function(mat) {
  segments <- get_segments(mat)
  # Determine the number of cores available for parallel processing
  no_of_cores <- 10
  
  results <- mclapply(segments, function(seg) {
    sub_mat <- extract_sub_mat(mat, seg)
    compute_entropy(sub_mat)
  }, mc.cores = no_of_cores)
  # results <- lapply(segments, function(seg) {
  #   sub_mat <- extract_sub_mat(mat, seg)
  #   compute_entropy(sub_mat)
  # })
  df_entropy <- as.data.frame(do.call(cbind, results))
  rownames(df_entropy) <- rownames(mat)
  colnames(df_entropy) <- paste0(segments)
  return(df_entropy)
}

#' Compute Central Coordinates of Boundaries
#'
#' @description
#' This function computes the central coordinates for each segment in the given 
#' Molecule Experiment (ME) object. The function takes advantage of the `boundaries` 
#' function to extract boundary coordinates and then computes the central coordinate 
#' by averaging the x and y locations.
#'
#' @param me A Molecule Experiment (ME) object that holds boundary data.
#' @param assay_name The name of the assay to extract boundaries from, defaults to "cell".
#' @return A dataframe with segment_id, sample_id, x_central, and y_central columns.
#' @examples
#' \dontrun{
#' # Assuming `me_obj` is your Molecule Experiment object
#' central_coordinates <- compute_central_coordinates(me_obj, assay_name = "cell")
#' }

compute_central_coordinates <- function(me, assay_name = "cell") {
  # Extract boundaries
  df <- MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE)
  
  # Calculate the central points for each segment_id
  result_df <- df %>%
    group_by(segment_id, sample_id) %>%
    summarise(
      x_central = mean(x_location),
      y_central = mean(y_location),
      .groups = "drop"
    )
  return(result_df)
}

#' Split Rectangle into Regions Based on Central Coordinates
#'
#' @description
#' This function determines which region each segment belongs to within a rectangle. 
#' It computes the central coordinates and then based on these central points, segments 
#' the rectangle into 12 regions. The function finally associates each segment to one 
#' of the 12 regions.
#'
#' @param me A Molecule Experiment (ME) object that holds boundary data.
#' @return A dataframe that associates each segment with a region within the rectangle.
#' @examples
#' \dontrun{
#' # Assuming `me_obj` is your Molecule Experiment object
#' segmented_regions <- split_rectangle(me_obj)
#' }
split_dataframe_by_region <- function(df, regions_segments) {
  # Create an empty data.frame with the same columns as df
  empty_df <- data.frame(lapply(df, function(col) numeric(0)))
  
  # Use mclapply to filter the dataframe based on segment_id and sample_id for each row in regions_segments
  df_list <- mclapply(seq_len(nrow(regions_segments)), function(i) {
    ids <- unlist(strsplit(regions_segments$segment_ids[i], ", "))
    samples <- unlist(strsplit(regions_segments$sample_ids[i], ", "))
    
    # Filter the dataframe based on segment_id and sample_id
    filtered_df <- df %>%
      filter(segment_id %in% ids & sample_id %in% samples) 
    
    if (nrow(filtered_df) == 0) {
      return(empty_df)
    } else {
      return(filtered_df)
    }
  }, mc.cores = 10) # Setting the number of cores to 10
  
  return(df_list)
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
  xmin = extent(me)[1]
  xmax = extent(me)[2]
  ymin = extent(me)[3]
  ymax = extent(me)[4]
  
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