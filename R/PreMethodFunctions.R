#' Extract Boundaries
#' 
#' Extracts boundaries from the ME object.
#' 
#' @param me An object containing information about the experiment.
#' 
#' @return A data frame containing boundary data.
extract_boundaries <- function(me) {
  df_boundary <- data.frame(MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE))
  return(df_boundary)
}

#' Calculate Centroids
#' 
#' Calculates centroids for given boundaries.
#' 
#' @param df_boundary A data frame containing boundary data.
#' 
#' @return A data frame containing centroids for each segment.
calculate_centroids <- function(df_boundary) {
  df_circle <- df_boundary %>%
    group_by(segment_id) %>%
    mutate(x_central = mean(head(x_location, -1)), 
           y_central = mean(head(y_location, -1))) %>%
    ungroup()
  return(df_circle)
}

#' Extract Boundaries and Calculate Centroids
#' 
#' A combined function to extract boundaries and calculate centroids from the ME object.
#' 
#' @param me An object containing information about the experiment.
#' 
#' @return A list containing boundary and centroid data frames.
extract_boundaries_and_centroids <- function(me) {
  df_boundary <- extract_boundaries(me)
  df_circle <- calculate_centroids(df_boundary)
  return(list(df_boundary = df_boundary, df_circle = df_circle))
}


#' Create Sectors
#' 
#' Creates sectors based on given data.
#' 
#' @param df A data frame containing the required fields to create sectors.
#' 
#' @return A data frame containing sector data.
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

#' Generate Scale Factors (All)
#' 
#' Generates scale factors for given value of k.
#' 
#' @param k Numeric value for scaling.
#' 
#' @return A vector of scale factors.
generate_scale_factors_all <- function(k) {
  scale_factors <- seq(from = 1/(k+1), to = 1, by = 1/(k+1))
  return(scale_factors)
}

#' Create Modified Sectors
#' 
#' Creates sectors with modifications based on given data.
#' 
#' @param df A data frame containing the required fields to create modified sectors.
#' 
#' @return A data frame containing modified sector data.
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

#' Create Sectors for Scaled Data
#' 
#' Creates sectors based on scaled data.
#' 
#' @param df A data frame containing the required fields.
#' @param scale Numeric value for scaling.
#' 
#' @return A data frame containing sector data based on scaled inputs.
create_sectors_for_scaled <- function(df, scale, k = 8) {
  df <- df %>%
    group_by(segment_id) %>%
    mutate(x_location = x_scaled, y_location = y_scaled, 
           angle = atan2(y_scaled - y_central, x_scaled - x_central) + pi) %>%
    ungroup()
  
  df_sectors <- df %>%
    group_by(segment_id) %>%
    group_split() %>%
    lapply(create_sectors_modified) %>%
    do.call(rbind, .) %>%
    arrange(sample_id) %>%
    mutate(concentric_id = paste0(segment_id, '_', scale*(k+1)))
  
  return(df_sectors)
}

#' Generate Scale Factors (All Outside)
#' 
#' Generates scale factors for given value of ks.
#' 
#' @param ks Numeric value for scaling.
#' 
#' @return A vector of scale factors for outer regions.
generate_scale_factors_all_outside <- function(ks) {
  scale_factors <- c(seq(from = 1, to = 2, by = 1/(ks+1)))
  return(scale_factors)
}


#' Generate Scaled Dataframe
#' 
#' Creates a dataframe with scaled polygon coordinates based on a given scale factor.
#' 
#' @param s Numeric scale factor.
#' @param df_circle Dataframe containing circle coordinates and segment IDs.
#' @param k Numeric value for scaling in the main function.
#' 
#' @return A dataframe with scaled coordinates.
create_scaled_df <- function(s, df_circle, k = 8) {
  df_circle %>%
    group_by(segment_id) %>%
    mutate(x_scaled = x_central + s * (x_location - x_central),
           y_scaled = y_central + s * (y_location - y_central),
           scale = s,
           concentric_id = paste0(segment_id, "_", s * (k + 1))) %>%
    ungroup()
}