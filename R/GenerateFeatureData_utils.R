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
