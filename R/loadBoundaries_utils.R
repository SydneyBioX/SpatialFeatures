#' Generate Feature Data Based on Assay Type Combining Various Polygon Assays
#'
#' @description
#' This comprehensive function processes the given data (`me`) and returns
#' feature data based on multiple assay types. The supported assays
#' encapsulated within this function include:
#'
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Sub-combo polygons
#' - Super-concentric polygons
#' - Super-combo polygons
#'
#' These assays represent various spatial subdivisions or transformations, each
#' addressing unique feature characteristics. The results from each assay are
#' combined and returned as a list of dataframes.
#'
#' @param me The input data object. This object contains spatial information,
#'        typically regarding segments or regions of interest.
#' @param featureTypes Character vector containing the feature types to
#' calculate.
#' @param k A numeric value indicating the scaling factor for concentric
#' polygons.
#'        Defaults to 5.
#' @return A list containing dataframes for each assay type:
#' - `subsector`: Feature data for sub-sector polygons.
#' - `subconcentric`: Feature data for sub-concentric polygons.
#' - `supersector`: Feature data for super-sector polygons.
#' - `superconcentric`: Feature data for super-concentric polygons.
#' @importFrom dplyr group_by mutate select ungroup rowwise distinct %>% bind_rows
#' @importFrom purrr map_dfr
#' @importFrom parallel mclapply
#' @importFrom rlang .data
#' @importFrom utils head
GenerateFeatureData <- function(me,
                                featureTypes = c("subsector", "subconcentric","supersector", "superconcentric"),
                                k = 5) {

  . <- NULL

  results <- extract_boundaries_and_centroids(me)
  df_circle <- results$df_circle

  featureTypes <- match.arg(featureTypes, several.ok = TRUE)

  if ("subsector" %in% featureTypes) {
    # For subsector
    subsector <- df_circle %>%
      dplyr::distinct(.data$segment_id, .data$x_location,
                      .data$y_location, .keep_all = TRUE) %>%
      dplyr::mutate(angle = atan2(.data$y_location - .data$y_central,
                                  .data$x_location - .data$x_central) + pi) %>%
      base::split(.$segment_id) %>%
      lapply(create_sectors) %>%
      dplyr::bind_rows() %>%
      dplyr::rename(x_section = .data$x_location_sector,
                    y_section = .data$y_location_sector,
                    area_id = .data$sector_id)
  } else {
    subsector <- NULL
  }

  if ("supersector" %in% featureTypes) {
    # For supersector
    supersector <- create_sector_df(df_circle)
  } else {
    supersector <- NULL
  }

  if ("subconcentric" %in% featureTypes) {
    common_scale_factors <- generate_scale_factors_all(k)
    common_scaled_df <- common_scale_factors %>%
      purrr::map_dfr(~create_scaled_df_sub(.x, df_circle, k))
    subconcentric <- common_scaled_df %>%
      dplyr::select(x_section = "x_scaled", y_section = "y_scaled",
                    "segment_id", "sample_id", area_id = "concentric_id")
    # modify the area_id column
    subconcentric <- subconcentric %>%
      dplyr::mutate(area_id = vapply(.data$area_id, modify_area_id, FUN.VALUE = "character")) # switched from sapply to vapply
  } else {
    subconcentric <- NULL
  }

  if ("superconcentric" %in% featureTypes) {
    # For superconcentric
    super_scale_factors <- generate_scale_factors_all_outside(k)
    super_scaled_df <- super_scale_factors %>%
      purrr::map_dfr(~create_scaled_df_super(.x, df_circle, k))
    superconcentric <- super_scaled_df %>%
      dplyr::select(x_section = "x_scaled", y_section = "y_scaled",
                    "segment_id", "sample_id", area_id = "concentric_id")
    # Modify the area_id column for each table
    superconcentric <- superconcentric %>%
      dplyr::mutate(area_id = vapply(.data$area_id, modify_area_id, FUN.VALUE = "character")) # switched from sapply to vapply
  } else {
    superconcentric <- NULL
  }

  featureData <- sapply(featureTypes, get, simplify = FALSE,
                        envir = environment(),
                        # FUN.VALUE = "list"
                        ) # changed from sapply to vapply # changed back

  return(featureData)
}

#' Create sectors from a given data frame
#'
#' @param df A data frame to create sectors from.
#' @return A data frame of sectors.
create_sectors <- function(df) {
  df <- df[order(df$angle),]
  sectors <- list()

  for (i in seq_len(nrow(df))) {
    p1 <- df[i,]
    p2 <- if (i == nrow(df)) df[1,] else df[i+1,]

    sector <- data.frame(
      x_location_sector = c(p1$x_location, p2$x_location,
                            p1$x_central, p1$x_location),
      y_location_sector = c(p1$y_location, p2$y_location,
                            p1$y_central, p1$y_location),
      segment_id = p1$segment_id,
      sample_id = p1$sample_id,
      sector_id = paste0(p1$segment_id, "_", sprintf("%02d", i))
    )
    sectors[[i]] <- sector
  }

  do.call(rbind, sectors)
}


#' Extract boundaries from a given object
#'
#' @param me The object to extract boundaries from.
#' @importFrom MoleculeExperiment boundaries
#' @return A data frame of extracted boundaries.
extract_boundaries <- function(me) {
  df_boundary <- data.frame(
    MoleculeExperiment::boundaries(me,assayName = "cell", flatten = TRUE)
    )
  return(df_boundary)
}

#' Calculate centroids from boundaries
#'
#' @param df_boundary A data frame of boundaries.
#' @return A data frame of centroids.
#' @importFrom dplyr %>% group_by mutate ungroup
calculate_centroids <- function(df_boundary) {
  df_circle <- df_boundary %>%
    dplyr::group_by(.data$segment_id) %>%
    dplyr::mutate(x_central = mean(utils::head(.data$x_location, -1)),
           y_central = mean(utils::head(.data$y_location, -1))) %>%
    dplyr::ungroup()
  return(df_circle)
}

#' Extract boundaries and calculate centroids
#'
#' @param me The object to process.
#' @return A list containing a data frame of boundaries and a data frame of
#' centroids.
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


#' Create a scaled data frame
#'
#' @param s A numeric scaling factor.
#' @param df_circle A data frame to scale.
#' @param k A numeric value, defaulting to 8.
#' @return A scaled data frame.
#' @importFrom dplyr group_by mutate select %>%
create_scaled_df_sub <- function(s, df_circle, k = 5) {
  scaling_factor <- s * (k + 1)

  df_circle %>%
    dplyr::group_by(.data$segment_id) %>%
    dplyr::mutate(
      x_diff = .data$x_location - .data$x_central,
      y_diff = .data$y_location - .data$y_central,
      x_scaled = .data$x_central + s * .data$x_diff,
      y_scaled = .data$y_central + s * .data$y_diff,
      scale = s,
      concentric_id = paste0(.data$segment_id, "_", scaling_factor)
    ) %>%
    dplyr::select(-"x_diff", -"y_diff") # Remove temporary columns
}

#' Create sector data frame from circle data
#'
#' This function takes a data frame of circle data with segment IDs and
#' locations, and returns a data frame of sectors with their scaled coordinates.
#'
#' @param df A data frame containing circle data with segment IDs, x and y
#' locations.
#' @return A data frame with the coordinates for each sector.
#' @importFrom dplyr %>% distinct group_by mutate ungroup
create_sector_df <- function(df) {

  # Distinct boundary points
  df <- df %>%
    dplyr::distinct(.data$segment_id, .data$x_location, .data$y_location,
                    .keep_all = TRUE) %>%
    dplyr::group_by(.data$segment_id) %>%
    dplyr::mutate(
      x_diff = .data$x_location - .data$x_central,
      y_diff = .data$y_location - .data$y_central,
      x_scaled = .data$x_central + 2 * .data$x_diff,
      y_scaled = .data$y_central + 2 * .data$y_diff
    ) %>%
    dplyr::ungroup()

  # Create sectors
  sectors <- list()
  df_grouped <- base::split(df, df$segment_id)

  for (segment in names(df_grouped)) {
    current_segment <- df_grouped[[segment]]
    n_rows <- nrow(current_segment)

    for (i in seq_len(n_rows-1)) { # We iterate till n_rows-2, skip the last point
      p1 <- current_segment[i,]
      p2 <- current_segment[i+1,]

      sector <- data.frame(
        x_location_sector = c(p1$x_location, p2$x_location, p2$x_scaled,
                              p1$x_scaled, p1$x_location),
        y_location_sector = c(p1$y_location, p2$y_location, p2$y_scaled,
                              p1$y_scaled, p1$y_location),
        segment_id = p1$segment_id,
        sample_id = p1$sample_id,
        sector_id = paste0(p1$segment_id, "_", sprintf("%02d", i))
      )
      sectors[[length(sectors) + 1]] <- sector
    }

    # To account for the last boundary point connecting with the first
    # without repetition
    p1 <- current_segment[n_rows,]
    p2 <- current_segment[1,]

    sector <- data.frame(
      x_location_sector = c(p1$x_location, p2$x_location, p2$x_scaled,
                            p1$x_scaled, p1$x_location),
      y_location_sector = c(p1$y_location, p2$y_location, p2$y_scaled,
                            p1$y_scaled, p1$y_location),
      segment_id = p1$segment_id,
      sample_id = p1$sample_id,
      sector_id = base::paste0(p1$segment_id, "_", base::sprintf("%02d", n_rows-1))
    )
    sectors[[length(sectors) + 1]] <- sector
  }

  df_sectors <- do.call(rbind, sectors)
  # Rename columns as per request
  df_sectors <- df_sectors %>%
    dplyr::rename(
      x_section = .data$x_location_sector,
      y_section = .data$y_location_sector,
      area_id = .data$sector_id
    )

  return(df_sectors)
}

#' Generate scale factors for all outside values of ks
#'
#' @param k A numeric value.
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
#' @importFrom dplyr %>% group_by mutate select
create_scaled_df_super <- function(s, df_circle, k = 5) {
  scaling_factor <- s * (k + 1)-k

  df_circle %>%
    dplyr::group_by(.data$segment_id) %>%
    dplyr::mutate(
      x_diff = .data$x_location - .data$x_central,
      y_diff = .data$y_location - .data$y_central,
      x_scaled = .data$x_central + s * .data$x_diff,
      y_scaled = .data$y_central + s * .data$y_diff,
      scale = s,
      concentric_id = paste0(.data$segment_id, "_", scaling_factor)
    ) %>%
    dplyr::select(-"x_diff", -"y_diff") # Remove temporary columns
}

#' Modify the area_id format
#'
#' This function takes an area_id string and modifies the last segment
#' if it's a single-digit number between 1 and 9 by adding a leading zero.
#'
#' @param id A string representing the area_id.
#' @return A string with the modified area_id format.
modify_area_id <- function(id) {
  # Use a regular expression to check if the pattern matches, and replace if so
  modified_id <- sub("_(\\d)$", "_0\\1", id)
  return(modified_id)
}
