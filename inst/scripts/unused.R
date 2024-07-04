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
