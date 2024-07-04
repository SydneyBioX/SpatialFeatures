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
  xmin = MoleculeExperiment::extent(me, assayName = "detected")[1]
  xmax = MoleculeExperiment::extent(me, assayName = "detected")[2]
  ymin = MoleculeExperiment::extent(me, assayName = "detected")[3]
  ymax = MoleculeExperiment::extent(me, assayName = "detected")[4]

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

#' Split Rectangle into Regions Based on Central Coordinates
#'
#' @description
#' This function determines which region each segment belongs to within a rectangle.
#' It computes the central coordinates and then based on these central points, segments
#' the rectangle into 12 regions. The function finally associates each segment to one
#' of the 12 regions.
#'
#' @param df A dataframe
#' @param regions_segments Regions segments
#' @param nCores Number of cores
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

