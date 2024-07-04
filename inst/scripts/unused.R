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
create_sectors_for_scaled <- function(df, scale, k = 5) {
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
generate_scale_factors_all_outside <- function(k) {
  scale_factors <- c(seq(from = 1, to = 2, by = 1/(k+1)))
  return(scale_factors)
}



#' Create a scaled data frame
#'
#' @param s A numeric scaling factor.
#' @param df_circle A data frame to scale.
#' @param k A numeric value, defaulting to 8.
#' @return A scaled data frame.
create_scaled_df_sub <- function(s, df_circle, k = 5) {
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




#' Create a scaled data frame
#'
#' @param s A numeric scaling factor.
#' @param df_circle A data frame to scale.
#' @param k A numeric value, defaulting to 8.
#' @return A scaled data frame.
create_scaled_df_super <- function(s, df_circle, k = 5) {
  scaling_factor <- s * (k + 1)-k

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



#' Delete columns ending with '0' from a matrix
#'
#' @param m A matrix to process.
#' @return A modified matrix without columns ending in '0'.
delete_inner <- function(m) {
  # Identify columns that don't end with "_0"
  cols_to_keep <- !grepl("_01$", colnames(m))

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
get_cells <- function(mat) {
  # Extract the unique cell identifiers before the last underscore
  cells <- unique(gsub("(.*)_[^_]+$", "\\1", colnames(mat)))
  return(cells)
}



#' Extract a submatrix corresponding to a segment from a matrix
#'
#' @param mat A matrix to process.
#' @param segment A character string indicating the segment to extract.
#' @return A submatrix of the provided matrix.
extract_cell_segments <- function(mat, cell) {
  # Create a regex pattern to match all segments of the cell
  pattern <- paste0("^", cell, "_\\d+$")
  matching_cols <- grepl(pattern, colnames(mat))
  cell_segments_mat <- mat[, matching_cols, drop = FALSE]
  return(cell_segments_mat)
}



#' Compute entropy for a given matrix
#'
#' @param mat A matrix to compute entropy for.
#' @return A vector of entropy values.
compute_cell_entropy <- function(cell_segments_mat) {
  # Convert to proportions
  proportions <- apply(cell_segments_mat, 1, function(counts) {
    total <- sum(counts)
    if (total == 0) {
      return(rep(0, length(counts)))
    } else {
      return(counts / total)
    }
  })

  # Compute entropy for each gene
  gene_entropies <- apply(proportions, 2, calculate_entropy)
  return(gene_entropies)
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
split_dataframe_by_region <- function(df, regions_segments, nCores = 1) {
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
  }, mc.cores = nCores) # Setting the number of cores to 10

  return(df_list)
}


#' Load Boundaries to Molecule Experiment Object
#'
#' @description
#' This function takes a Molecule Experiment (ME) object as input and enriches it
#' with boundary data from a variety of assays. The currently supported assays are:
#'
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Sub-combo polygons
#' - Super-concentric polygons
#' - Super-combo polygons
#'
#' The function internally processes each assay, extracts respective boundary information,
#' and updates the original ME object with this new data.
#'
#' @param me A Molecule Experiment (ME) object that needs its boundaries populated based
#'        on the assays mentioned above.
#' @return An enriched Molecule Experiment (ME) object that includes boundary information
#'         for the supported assays.
#' @examples
#' data(small_me, package = "MoleculeExperiment")
#' me_list <- loadBoundariesToME(small_me)

loadBoundariesToME <- function(me, k = 5, r = 6) {

  regions_segments <- split_rectangle(me)

  # # List of assay names to process
  # assay_names <- c("sub_sector", "sub_concentric", "sub_combo",
  #                  "super_sector", "super_concentric", "super_combo")
  assay_names <- c("sub_sector", "sub_concentric",
                   "super_sector", "super_concentric")

  # Generate all the data using the GenerateFeatureData function
  results_all <- SpatialFeatures::GenerateFeatureData(me, k = k)

  df_boundaries <- MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE)

  # Read the molecule data
  df_molecule <- MoleculeExperiment::molecules(me, assayName = "detected", flatten = TRUE)

  process_segment <- function(i) {
    # List to store the results for each assay name
    results_list <- list()
    moleculesMEList_current <- list()

    boundariesMEList_current <- list()
    df_bdy = split_dataframe_by_region(df_boundaries, regions_segments)[[i]]

    boundariesMEList_current <- dataframeToMEList(df = df_bdy, dfType = "boundaries", assayName = "cell",
                                                  sampleCol = "sample_id", factorCol = "segment_id",
                                                  xCol = "x_location", yCol = "y_location",
                                                  keepCols = "essential", scaleFactor = 1)

    for (assayName in assay_names) {

      # Split the dataframe based on region
      df_tmp <- split_dataframe_by_region(results_all[[assayName]], regions_segments)[[i]]

      # Check if the data frame is empty
      if (nrow(df_tmp) == 0) next

      # Modify column names based on the given logic
      df_tmp$segment_id <- df_tmp$area_id
      df_tmp$x_location <- df_tmp$x_section
      df_tmp$y_location <- df_tmp$y_section

      # Filter the location of Molecules data
      if (assayName == "super_sector") {
        coord_ranges <- df_tmp %>%
          dplyr::group_by(sample_id) %>%
          dplyr::summarise(
            xmin = min(x_location, na.rm = TRUE),
            xmax = max(x_location, na.rm = TRUE),
            ymin = min(y_location, na.rm = TRUE),
            ymax = max(y_location, na.rm = TRUE),
            .groups = "drop"
          )

        filtered_molecule <- df_molecule %>%
          dplyr::inner_join(coord_ranges, by = "sample_id") %>%
          dplyr::filter(x_location >= xmin & x_location <= xmax &
                          y_location >= ymin & y_location <= ymax)

        filtered_molecule$x_coords <- filtered_molecule$x_location
        filtered_molecule$y_coords <- filtered_molecule$y_location

        moleculesMEList_current <- dataframeToMEList(filtered_molecule, dfType = "molecules", assayName = "detected",
                                                     sampleCol = "sample_id", factorCol = "features", xCol = "x_coords",
                                                     yCol = "y_coords")
      }

      # Convert to ME List using dataframeToMEList function
      results_list[[assayName]] <- dataframeToMEList(df = df_tmp, dfType = "boundaries", assayName = assayName,
                                                     sampleCol = "sample_id", factorCol = "segment_id",
                                                     xCol = "x_location", yCol = "y_location",
                                                     keepCols = "essential", scaleFactor = 1)
    }

    me_copy <- MoleculeExperiment(
      molecules = moleculesMEList_current,
      boundaries = boundariesMEList_current
    )

    # Set the boundaries for each assay in me_copy
    for (assay in assay_names) {
      boundaries(me_copy, assay) <- results_list[[assay]]
    }

    return(me_copy)
  }

  # Only process when i = 6
  me_copy_r <- process_segment(r)

  return(me_copy_r)
}
