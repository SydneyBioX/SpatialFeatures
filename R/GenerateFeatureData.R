#' Generate Feature Data Based on Assay Name
#'
#' This function processes the input data according to the specified assay name 
#' and returns the respective feature data.
#' 
#' @param me A data object that contains the input data.
#' @param assayName A character string specifying the assay name. It must be one 
#'   of "sub-sector", "sub-concentric", "sub-combo", "super-concentric", or 
#'   "super-combo".
#'
#' @return A dataframe containing the processed feature data.
#' @export
GenerateFeatureData <- function(me, assayName, k = 8) {
  # Ensure assayName is valid
  if (!assayName %in% c("sub-sector", "sub-concentric", 'sub-combo', 
                        "super-concentric", "super-combo")) {
    stop("Invalid assayName provided!")
  }
  
  # Common function for extracting boundaries and centroids
  results <- extract_boundaries_and_centroids(me)
  
  # Sector-polygon
  if (assayName == "sub-sector") {
    # Extract boundaries and centroids using the function
    df_ctv = results$df_circle
    
    # Remove duplicate boundary points for each segment_id
    df_ctv <- df_ctv %>%
      group_by(segment_id, x_location, y_location) %>%
      #slice(1) %>%
      ungroup()
    
    # Compute angles
    df_ctv <- df_ctv %>%
      rowwise() %>%
      mutate(angle = atan2(y_location - y_central, x_location - x_central) + pi)
    
    # Split the dataframe by segment_id and create sectors for each segment
    df_sectors <- df_ctv %>%
      group_by(segment_id) %>%
      group_split() %>%
      lapply(create_sectors) %>%
      do.call(rbind, .) %>% 
      dplyr::rename(x_section = x_location_sector, y_section = y_location_sector, area_id = sector_id) %>%
      arrange(sample_id)
    return(df_sectors)
  }
  
  # Concentric-polygon
  if (assayName == "sub-concentric") {
    # Extract boundaries and centroids using the function
    df_boundary = results$df_boundary
    df_circle = results$df_circle
    
    # Define the number of polygons
    scale_factors <- generate_scale_factors_all(k)
    
    # Generate a dataframe with the scaled polygons
    df_concentric_initial <- map_dfr(scale_factors, ~create_scaled_df(.x, df_circle, k))
    
    # Select required columns
    df_concentric <- df_concentric_initial %>%
      select(x_scaled, y_scaled, segment_id, sample_id, concentric_id) %>%
      dplyr::rename(x_section = x_scaled, y_section = y_scaled, area_id = concentric_id) %>%
      arrange(sample_id)
    return(df_concentric)
  }
  
  # Combo-polygon
  if (assayName == "sub-combo") {
    # Extract boundaries and centroids using the function
    df_boundary = results$df_boundary
    df_circle = results$df_circle
    
    # Define the number of polygons
    scale_factors <- generate_scale_factors_all(k)
    
    # Generate a dataframe with the scaled polygons
    df_concentric_initial <- map_dfr(scale_factors, ~create_scaled_df(.x, df_circle, k))
    
    # Remove duplicate boundary points for each segment_id
    df_concentric_initial <- df_concentric_initial %>%
      group_by(concentric_id, x_scaled, y_scaled) %>%
      #slice(1) %>%
      ungroup()
    
    # Generate concentric polygons and sectors for each of them
    df_combo <- df_concentric_initial %>%
      group_by(concentric_id) %>%
      group_split() %>%
      lapply(function(df) {
        create_sectors_for_scaled(df, unique(df$scale))
      }) %>%
      do.call(rbind, .) %>%
      select(x_location_sector, y_location_sector, segment_id, sample_id, combo_id) %>%
      dplyr::rename(x_section = x_location_sector, y_section = y_location_sector, area_id = combo_id) %>%
      arrange(sample_id)
    return(df_combo)
  }
  
  # Concentric-polygon-super
  if (assayName == "super-concentric") {
    # Extract boundaries and centroids using the function
    results_super = extract_boundaries_and_centroids(me)
    df_boundary_super = results_super$df_boundary
    df_circle_super = results_super$df_circle
    
    # Define the number of polygons
    scale_factors_outside <- generate_scale_factors_all_outside(k)
    
    # Generate a dataframe with the scaled polygons for outside boundaries
    df_concentric_super_initial <- map_dfr(scale_factors_outside, ~create_scaled_df(.x, df_circle_super, k))
    
    # Select required columns
    df_concentric_super <- df_concentric_super_initial %>%
      select(x_scaled, y_scaled, segment_id, sample_id, concentric_id) %>%
      dplyr::rename(x_section = x_scaled, y_section = y_scaled, area_id = concentric_id) %>%
      arrange(sample_id)
    return(df_concentric_super)
  }
  
  # Combo-polygon-super
  if (assayName == "super-combo") {
    # Extract boundaries and centroids using the function
    results_super = extract_boundaries_and_centroids(me)
    df_boundary_super = results_super$df_boundary
    df_circle_super = results_super$df_circle
    
    # Define the number of polygons
    scale_factors_outside <- generate_scale_factors_all_outside(k)
    
    # Generate a dataframe with the scaled polygons for outside boundaries
    df_concentric_super_initial <- map_dfr(scale_factors_outside, ~create_scaled_df(.x, df_circle_super, k))
    
    # Remove duplicate boundary points for each segment_id
    df_concentric_super_initial <- df_concentric_super_initial %>%
      group_by(concentric_id, x_scaled, y_scaled) %>%
      #slice(1) %>%
      ungroup()
    
    # Generate concentric polygons and sectors for each of them
    df_combo_super <- df_concentric_super_initial %>%
      group_by(concentric_id) %>%
      group_split() %>%
      lapply(function(df) {
        create_sectors_for_scaled(df, unique(df$scale))
      }) %>%
      do.call(rbind, .) %>%
      select(x_location_sector, y_location_sector, segment_id, sample_id, combo_id) %>%
      dplyr::rename(x_section = x_location_sector, y_section = y_location_sector, area_id = combo_id) %>%
      arrange(sample_id)
    return(df_combo_super)
  }
}